"""Merge per-chunk InterProScan TSV outputs into one annotation table.

Reads ``data/interim/interproscan/output/chunk_*.tsv`` (produced by
``scripts/run_interproscan.slurm``), normalizes the variable-width rows to
the 15-column InterProScan layout (TSV columns depend on the -goterms /
-iprlookup / -pa flags), and writes a single headered table to
``data/processed/interproscan_domains.tsv``.

Run with:

    pixi run python -m bm_ppi.merge_interproscan
"""

from __future__ import annotations

from collections import Counter
import csv
from pathlib import Path

from loguru import logger
import typer

from bm_ppi.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR

app = typer.Typer()

COLUMNS = [
    "protein_id",
    "seq_md5",
    "seq_length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "stop",
    "evalue",
    "status",
    "run_date",
    "interpro_accession",
    "interpro_description",
    "go_terms",
    "pathways",
]


@app.command()
def main(
    input_dir: Path = INTERIM_DATA_DIR / "interproscan" / "output",
    manifest_path: Path = INTERIM_DATA_DIR / "interproscan" / "input" / "manifest.tsv",
    output_path: Path = PROCESSED_DATA_DIR / "interproscan_domains.tsv",
) -> None:
    chunk_files = sorted(input_dir.glob("chunk_*.tsv"))
    if not chunk_files:
        raise typer.BadParameter(f"No chunk_*.tsv files found in {input_dir}")

    # Report chunks that have not finished yet (missing output files).
    if manifest_path.is_file():
        with open(manifest_path, newline="") as fh:
            expected = {row["chunk"] for row in csv.DictReader(fh, delimiter="\t")}
        done = {p.stem for p in chunk_files}
        missing = sorted(expected - done)
        if missing:
            logger.warning(
                f"{len(missing)}/{len(expected)} chunks have no output yet, e.g. {missing[:5]}"
            )

    rows: list[list[str]] = []
    for chunk_file in chunk_files:
        with open(chunk_file, newline="") as fh:
            for fields in csv.reader(fh, delimiter="\t"):
                if len(fields) < 11:
                    logger.warning(f"Skipping malformed line in {chunk_file.name}: {fields[:3]}")
                    continue
                rows.append(fields[: len(COLUMNS)] + [""] * (len(COLUMNS) - len(fields)))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(COLUMNS)
        writer.writerows(rows)

    proteins_with_hits = {row[0] for row in rows}
    per_analysis = Counter(row[3] for row in rows)
    unique_ipr = {row[11] for row in rows if len(row) > 11 and row[11] != "-"}

    logger.success(f"Merged {len(chunk_files)} chunks -> {output_path}")
    logger.info(
        f"{len(rows)} annotation lines covering {len(proteins_with_hits)} proteins; "
        f"{len(unique_ipr)} unique InterPro entries"
    )
    for analysis, n in per_analysis.most_common():
        logger.info(f"  {analysis}: {n} hits")


if __name__ == "__main__":
    app()
