"""Validate and split the Bombyx mori proteome FASTA for InterProScan.

Reads the working protein set
(``data/raw/GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa``),
checks IDs and alphabet, then writes fixed-size chunk FASTAs plus a
manifest to ``data/interim/interproscan/input/`` so a SLURM array job can
process one chunk per task.

Run with:

    pixi run python -m bm_ppi.prepare_interproscan_input
"""

from __future__ import annotations

from collections import Counter
import csv
from pathlib import Path

from Bio import SeqIO
from loguru import logger
import typer

from bm_ppi.config import INTERIM_DATA_DIR, RAW_DATA_DIR

app = typer.Typer()

# Standard 20 amino acids. Ambiguity codes (X, B, Z) and selenocysteine /
# pyrrolysine (U, O) are tolerated by HMMER-based analyses; '*' and 'J'
# are reported as problematic because some analyses reject them.
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
TOLERATED_AA = STANDARD_AA | set("XBZUO")


@app.command()
def main(
    input_fasta: Path = RAW_DATA_DIR / "GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa",
    output_dir: Path = INTERIM_DATA_DIR / "interproscan" / "input",
    chunk_size: int = 500,
) -> None:
    records = list(SeqIO.parse(str(input_fasta), "fasta"))
    logger.info(f"Read {len(records)} sequences from {input_fasta}")

    # ---- Validation ----
    ids = [rec.id for rec in records]
    duplicated = [seq_id for seq_id, n in Counter(ids).items() if n > 1]
    if duplicated:
        raise typer.BadParameter(
            f"{len(duplicated)} duplicated sequence IDs, e.g. {duplicated[:5]}"
        )

    char_counts: Counter[str] = Counter()
    n_bad = 0
    lengths = []
    for rec in records:
        seq = str(rec.seq).upper()
        lengths.append(len(seq))
        nonstandard = set(seq) - STANDARD_AA
        if nonstandard:
            n_bad += 1
            char_counts.update(nonstandard)
        if not seq:
            raise typer.BadParameter(f"Empty sequence: {rec.id}")

    logger.info(f"Lengths: min={min(lengths)}, max={max(lengths)}, total_residues={sum(lengths)}")
    if char_counts:
        logger.warning(f"{n_bad} sequences contain non-standard residues: {dict(char_counts)}")
        fatal = set(char_counts) - TOLERATED_AA
        if fatal:
            raise typer.BadParameter(
                f"Unsupported characters {sorted(fatal)}; clean the FASTA first."
            )
    else:
        logger.info("All sequences use the standard 20-aa alphabet")

    # ---- Split into chunks ----
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "manifest.tsv"
    chunks = [records[i : i + chunk_size] for i in range(0, len(records), chunk_size)]

    with open(manifest_path, "w", newline="") as manifest_fh:
        writer = csv.writer(manifest_fh, delimiter="\t")
        writer.writerow(["chunk", "file", "n_seqs", "n_residues"])
        for idx, chunk in enumerate(chunks, start=1):
            chunk_name = f"chunk_{idx:04d}"
            chunk_path = output_dir / f"{chunk_name}.faa"
            n_written = SeqIO.write(chunk, str(chunk_path), "fasta")
            writer.writerow(
                [chunk_name, chunk_path.name, n_written, sum(len(r.seq) for r in chunk)]
            )

    logger.success(f"Wrote {len(chunks)} chunk FASTAs to {output_dir}; manifest: {manifest_path}")
    logger.info(
        f"Submit the array job with: sbatch --array=1-{len(chunks)} scripts/run_interproscan.slurm"
    )


if __name__ == "__main__":
    app()
