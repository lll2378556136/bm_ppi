"""Incremental protein segmentation over completed DPAM chunks.

DPAM domain prediction runs per chunk under ``data/interim/DPAM_runs/full_v1``.
Rather than waiting for all chunks to finish, this module picks up every chunk
whose ``validation.json`` reports ``hard_success == true``, merges their
``*_domains`` tables, restricts the protein FASTA to those chunks' target IDs
(from ``*_struc.list``), and delegates segmentation to the reference
``protein_segmenter`` in ``third_party/pirna_ppi``.

Proteins absent from the merged DPAM domains table (zero-hit or
circuit-breaker-terminated proteins) automatically fall back to InterProScan
domain annotations inside the segmenter. The run is idempotent: re-run it after
more chunks finish to extend coverage.

Run with:

    pixi run python -m bm_ppi.segment_proteins
    pixi run python -m bm_ppi.segment_proteins --chunks 1,2,33 --max-segment-length 600
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from Bio import SeqIO
from loguru import logger
import typer

from bm_ppi.config import INTERIM_DATA_DIR, PROJ_ROOT, RAW_DATA_DIR

app = typer.Typer(add_completion=False)

PIRNA_PPI_DIR = PROJ_ROOT / "third_party" / "pirna_ppi"
DEFAULT_FASTA = RAW_DATA_DIR / "GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa"
DEFAULT_CHUNKS_TSV = INTERIM_DATA_DIR / "DPAM_runs" / "full_v1" / "chunks.tsv"
DEFAULT_IPRSCAN_DIR = INTERIM_DATA_DIR / "interproscan" / "output"
DEFAULT_STRUCTURE_DIR = INTERIM_DATA_DIR / "DPAM_input" / "bm_ppi"
DEFAULT_OUTPUT_DIR = INTERIM_DATA_DIR / "segmented_proteins"


def _import_segmenter():
    """Import the reference segmenter after adding pirna_ppi to sys.path.

    pirna_ppi.config also calls ``logger.remove(0)`` on import; since
    bm_ppi.config already removed that handler, the remove is made tolerant
    and any handlers pirna_ppi adds are dropped again to keep single-sink
    logging.
    """
    if str(PIRNA_PPI_DIR) not in sys.path:
        sys.path.insert(0, str(PIRNA_PPI_DIR))

    original_remove = logger.remove
    handlers_before = set(logger._core.handlers)

    def _tolerant_remove(*args, **kwargs):
        try:
            return original_remove(*args, **kwargs)
        except ValueError:
            return None

    try:
        logger.remove = _tolerant_remove
        from pirna_ppi.protein_segmenter import cli as seg_cli
        from pirna_ppi.protein_segmenter.config import SegmenterConfig
    finally:
        logger.remove = original_remove
        for handler_id in set(logger._core.handlers) - handlers_before:
            original_remove(handler_id)

    # polars infers the protein-ID columns as Int64 for our all-numeric IDs,
    # which would never match the string IDs parsed from FASTA headers.
    # Re-key both domain dicts as strings before they reach run_segmentation.
    def _strkeyed(parser):
        def wrapped(*args, **kwargs):
            return {str(k): v for k, v in parser(*args, **kwargs).items()}
        return wrapped

    seg_cli.parse_dpam_domains = _strkeyed(seg_cli.parse_dpam_domains)
    seg_cli.parse_iprscan = _strkeyed(seg_cli.parse_iprscan)

    return seg_cli.run_segmentation, SegmenterConfig


def read_completed_chunks(chunks_tsv: Path) -> dict[int, dict]:
    """Map chunk_id -> {dataset, input_dir} for chunks with hard_success validation."""
    import csv

    completed = {}
    with open(chunks_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            input_dir = Path(row["input_dir"])
            validation_path = input_dir / "validation.json"
            if not validation_path.is_file():
                continue
            validation = json.loads(validation_path.read_text())
            if validation.get("hard_success"):
                completed[int(row["chunk_id"])] = {
                    "dataset": row["dataset"],
                    "input_dir": input_dir,
                    "n_targets": int(row["n_targets"]),
                }
    return completed


def read_structure_list(path: Path) -> list[str]:
    """Read a DPAM ``*_struc.list`` protein ID list."""
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def stage_inputs(
    chunk_ids: list[int],
    completed: dict[int, dict],
    fasta_path: Path,
    iprscan_dir: Path,
    staging_dir: Path,
) -> tuple[Path, Path, Path, dict[str, int]]:
    """Build the three segmenter input files restricted to selected chunks.

    Returns (subset_fasta, merged_domains, merged_iprscan, protein_to_chunk).
    """
    staging_dir.mkdir(parents=True, exist_ok=True)
    for old in staging_dir.iterdir():
        old.unlink()

    protein_to_chunk: dict[str, int] = {}
    domains_parts: list[str] = []
    for chunk_id in chunk_ids:
        info = completed[chunk_id]
        input_dir, dataset = info["input_dir"], info["dataset"]
        struc_list = input_dir / f"{dataset}_struc.list"
        domains_path = input_dir / f"{dataset}_domains"
        for pid in read_structure_list(struc_list):
            protein_to_chunk[pid] = chunk_id
        # Merge domains tables, keeping the header only once
        lines = domains_path.read_text().splitlines()
        domains_parts.extend(lines if not domains_parts else lines[1:])

    merged_domains = staging_dir / "dpam_domains.tsv"
    merged_domains.write_text("\n".join(domains_parts) + "\n")

    subset_fasta = staging_dir / "proteins.fasta"
    n_written = 0
    with open(subset_fasta, "w") as out:
        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.id in protein_to_chunk:
                SeqIO.write(record, out, "fasta")
                n_written += 1
    if n_written != len(protein_to_chunk):
        raise ValueError(
            f"FASTA subset mismatch: wrote {n_written} records but "
            f"struc.list files list {len(protein_to_chunk)} proteins"
        )

    merged_iprscan = staging_dir / "iprscan.tsv"
    with open(merged_iprscan, "w") as out:
        for chunk_file in sorted(iprscan_dir.glob("chunk_*.tsv")):
            out.write(chunk_file.read_text())

    return subset_fasta, merged_domains, merged_iprscan, protein_to_chunk


def write_provenance(
    output_dir: Path,
    chunk_ids: list[int],
    completed: dict[int, dict],
    protein_to_chunk: dict[str, int],
    detail_path: Path,
) -> None:
    """Write per-chunk and per-protein audit tables for this run."""
    import csv

    # Domain source per protein from the segmenter's detail table
    domain_source: dict[str, str] = {}
    with open(detail_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n_domains = int(row["n_domains"])
            status = row["status"]
            if status == "no_domain" or n_domains == 0:
                domain_source[row["protein_id"]] = "none"
            else:
                domain_source[row["protein_id"]] = "annotated"

    with open(output_dir / "processed_chunks.tsv", "w") as fh:
        fh.write("chunk_id\tdataset\tn_targets\n")
        for chunk_id in chunk_ids:
            info = completed[chunk_id]
            fh.write(f"{chunk_id}\t{info['dataset']}\t{info['n_targets']}\n")

    with open(output_dir / "processed_proteins.tsv", "w") as fh:
        fh.write("protein_id\tchunk_id\tdomain_status\n")
        for pid, chunk_id in sorted(protein_to_chunk.items(), key=lambda x: x[1]):
            fh.write(f"{pid}\t{chunk_id}\t{domain_source.get(pid, 'unprocessed')}\n")


@app.command()
def main(
    chunks: str = typer.Option(
        "auto", help="Comma-separated chunk IDs, or 'auto' for all validated chunks."
    ),
    fasta_path: Path = typer.Option(DEFAULT_FASTA, help="Full proteome FASTA."),
    chunks_tsv: Path = typer.Option(DEFAULT_CHUNKS_TSV, help="DPAM chunk manifest."),
    iprscan_dir: Path = typer.Option(DEFAULT_IPRSCAN_DIR, help="InterProScan chunk output dir."),
    structure_dir: Path = typer.Option(DEFAULT_STRUCTURE_DIR, help="PDB/PAE structure dir."),
    output_dir: Path = typer.Option(DEFAULT_OUTPUT_DIR, help="Segmentation output dir."),
    max_segment_length: int = typer.Option(600, help="Length threshold for segmentation."),
) -> None:
    run_segmentation, SegmenterConfig = _import_segmenter()

    completed = read_completed_chunks(chunks_tsv)
    logger.info(f"Found {len(completed)} validated DPAM chunks: {sorted(completed)}")

    if chunks == "auto":
        chunk_ids = sorted(completed)
    else:
        chunk_ids = sorted(int(c) for c in chunks.split(","))
        invalid = [c for c in chunk_ids if c not in completed]
        if invalid:
            raise typer.BadParameter(
                f"Chunks {invalid} are not validated-complete (available: {sorted(completed)})"
            )
    if not chunk_ids:
        raise typer.BadParameter("No completed DPAM chunks to process.")

    logger.info(f"Segmenting proteins from chunks: {chunk_ids}")
    staging_dir = output_dir / "_staging"
    subset_fasta, merged_domains, merged_iprscan, protein_to_chunk = stage_inputs(
        chunk_ids, completed, fasta_path, iprscan_dir, staging_dir
    )
    logger.info(f"Staged {len(protein_to_chunk)} proteins in {staging_dir}")

    config = SegmenterConfig(
        fasta_path=subset_fasta,
        dpam_domains_path=merged_domains,
        iprscan_path=merged_iprscan,
        structure_dir=structure_dir,
        output_dir=output_dir,
        max_segment_length=max_segment_length,
    )
    run_segmentation(config)

    write_provenance(
        output_dir, chunk_ids, completed, protein_to_chunk,
        output_dir / "segment_detail.tsv",
    )
    logger.success("Incremental segmentation complete.")


if __name__ == "__main__":
    app()
