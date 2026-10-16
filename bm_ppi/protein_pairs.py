"""Expand Stage 1 CC-compatible GeneID pairs into segmented protein pairs.

Only the proteins present in the four segmenter FASTA files are eligible.
The Stage 1 table is expected to contain one row per unordered protein pair;
its rows are streamed without retaining the expanded pairs in memory.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import gzip
from itertools import product
from pathlib import Path
import re
from typing import Annotated

from Bio import SeqIO
from loguru import logger
import typer

from bm_ppi.config import INTERIM_DATA_DIR

app = typer.Typer(add_completion=False)

DEFAULT_PAIRS = (
    INTERIM_DATA_DIR
    / "uniprot"
    / "cc_spatial_compatibility"
    / "stage1"
    / "stage1_spatially_compatible_pairs.tsv.gz"
)
DEFAULT_INPUT_DIR = INTERIM_DATA_DIR / "segmented_proteins"
DEFAULT_OUTPUT_DIR = INTERIM_DATA_DIR / "protein_pairs_stage1"
INPUT_FASTAS = (
    "segmented.fasta",
    "unsegmented.fasta",
    "no_domain.fasta",
    "oversized_segments.fasta",
)
OUTPUT_FILES = (
    "passed_sequences.fasta",
    "failed_sequences.fasta",
    "protein_pairs.fasta.gz",
    "summary.tsv",
    "processed_chunks.tsv",
)
SEGMENT_ID = re.compile(r"^(.+)-seg([1-9][0-9]*)$")
REQUIRED_COLUMNS = {"protein_A", "protein_B", "pair_class"}


@dataclass(frozen=True)
class Segment:
    segment_id: str
    sequence: str


def load_segments(
    input_dir: Path, max_length: int
) -> tuple[
    dict[str, list[Segment]], list[tuple[Segment, str]], list[tuple[Segment, str]], Counter
]:
    """Validate all four FASTAs and index only length-qualified segments."""
    if not input_dir.is_dir():
        raise ValueError(f"Segmentation directory does not exist: {input_dir}")
    indexed: dict[str, list[Segment]] = defaultdict(list)
    passed: list[tuple[Segment, str]] = []
    failed: list[tuple[Segment, str]] = []
    seen: set[str] = set()
    all_parents: set[str] = set()

    for filename in INPUT_FASTAS:
        path = input_dir / filename
        if not path.is_file():
            raise ValueError(f"Missing required segmentation FASTA: {path}")
        count = 0
        with path.open("r", encoding="utf-8") as handle:
            first_line = handle.readline()
            if first_line and not first_line.startswith(">"):
                raise ValueError(f"Invalid FASTA header in {path}")
            handle.seek(0)
            for record in SeqIO.parse(handle, "fasta"):
                count += 1
                match = SEGMENT_ID.fullmatch(record.id)
                if match is None:
                    raise ValueError(f"Invalid segment ID in {path}: {record.id!r}")
                if record.id in seen:
                    raise ValueError(f"Duplicate segment ID: {record.id}")
                seen.add(record.id)
                parent = match.group(1)
                all_parents.add(parent)
                indexed.setdefault(parent, [])
                sequence = str(record.seq).upper()
                if not sequence or re.fullmatch(r"[A-Z*]+", sequence) is None:
                    raise ValueError(f"Empty or invalid sequence for {record.id} in {path}")
                segment = Segment(record.id, sequence)
                if len(sequence) <= max_length:
                    indexed[parent].append(segment)
                    passed.append((segment, filename))
                else:
                    failed.append((segment, filename))
        logger.info(f"{filename}: {count} segments")

    if not seen:
        raise ValueError("No segments found in the four FASTA files")
    stats = Counter(
        segments_total=len(seen),
        segments_passed=len(passed),
        segments_failed=len(failed),
        parents_total=len(all_parents),
        parents_passed=sum(bool(segments) for segments in indexed.values()),
        parents_failed_only=sum(not segments for segments in indexed.values()),
    )
    return indexed, passed, failed, stats


def check_pairs_header(pairs_path: Path) -> None:
    """Check the compressed Stage 1 table before creating output files."""
    if not pairs_path.is_file():
        raise ValueError(f"Stage 1 pair table does not exist: {pairs_path}")
    with gzip.open(pairs_path, "rt", encoding="utf-8", newline="") as handle:
        header = csv.reader(handle, delimiter="\t")
        fields = next(header, [])
        if len(fields) != len(set(fields)) or not REQUIRED_COLUMNS.issubset(fields):
            raise ValueError(
                f"Invalid Stage 1 columns: need {sorted(REQUIRED_COLUMNS)} without duplicates"
            )


def write_fasta(records: list[tuple[Segment, str]], path: Path) -> None:
    """Write a filtered FASTA with the reference script's source annotation."""
    with path.open("x", encoding="utf-8") as handle:
        for segment, source in records:
            handle.write(f">{segment.segment_id} source={source}\n{segment.sequence}\n")


def expand_pairs(
    pairs_path: Path, indexed: dict[str, list[Segment]], output_path: Path, stats: Counter
) -> None:
    """Stream the validated parent table into a compressed fragment-pair FASTA."""
    with (
        gzip.open(pairs_path, "rt", encoding="utf-8", newline="") as source,
        gzip.open(output_path, "xt", encoding="utf-8", compresslevel=6) as output,
    ):
        reader = csv.DictReader(source, delimiter="\t")
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(row[field] is None for field in REQUIRED_COLUMNS):
                raise ValueError(f"Malformed Stage 1 row at line {line_number}")
            a, b = row["protein_A"].strip(), row["protein_B"].strip()
            if not a or not b or a == b or row["pair_class"].strip() != "Same locality":
                raise ValueError(f"Invalid Stage 1 pair at line {line_number}: {a!r}, {b!r}")
            stats["pairs_read"] += 1
            segments_a, segments_b = indexed.get(a), indexed.get(b)
            if a not in indexed or b not in indexed:
                stats["pairs_missing_parent"] += 1
                continue
            if not segments_a or not segments_b:
                stats["pairs_length_filtered"] += 1
                continue
            stats["pairs_expandable"] += 1
            for segment_a, segment_b in product(segments_a, segments_b):
                output.write(
                    f">{segment_a.segment_id}_{segment_b.segment_id}\n"
                    f"{segment_a.sequence}:{segment_b.sequence}\n"
                )
                stats["fragment_pairs_written"] += 1


def run(
    pairs_path: Path,
    input_dir: Path,
    output_dir: Path,
    max_length: int = 600,
    overwrite: bool = False,
) -> Counter:
    """Generate filtered FASTAs, compressed fragment pairs, and a count summary."""
    if max_length < 1:
        raise ValueError("max_length must be positive")
    indexed, passed, failed, stats = load_segments(input_dir, max_length)
    check_pairs_header(pairs_path)
    if output_dir.resolve() == input_dir.resolve():
        raise ValueError("Output directory must differ from segmentation input directory")
    existing = [name for name in OUTPUT_FILES if (output_dir / name).exists()]
    if existing and not overwrite:
        raise ValueError(f"Output files already exist in {output_dir}; use --overwrite")
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in existing:
        (output_dir / name).unlink()
    write_fasta(passed, output_dir / OUTPUT_FILES[0])
    write_fasta(failed, output_dir / OUTPUT_FILES[1])
    expand_pairs(pairs_path, indexed, output_dir / OUTPUT_FILES[2], stats)

    # Copy the segmenter's chunk manifest so each run records its coverage
    chunks_path = input_dir / "processed_chunks.tsv"
    n_chunks = 0
    if chunks_path.is_file():
        lines = chunks_path.read_text().splitlines()
        (output_dir / OUTPUT_FILES[4]).write_text("\n".join(lines) + "\n")
        n_chunks = max(len(lines) - 1, 0)
    else:
        logger.warning(f"No processed_chunks.tsv in {input_dir}; skipping chunk provenance")

    with (output_dir / OUTPUT_FILES[3]).open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["category", "n_entries"])
        writer.writerows(
            (key, stats[key])
            for key in (
                "segments_total",
                "segments_passed",
                "segments_failed",
                "parents_total",
                "parents_passed",
                "parents_failed_only",
                "pairs_read",
                "pairs_expandable",
                "pairs_missing_parent",
                "pairs_length_filtered",
                "fragment_pairs_written",
            )
        )
        writer.writerow(["chunks_covered", n_chunks])
        writer.writerow(["stage1_table", pairs_path.name])
    logger.info(f"Stage 1 rows: {stats['pairs_read']}; expandable: {stats['pairs_expandable']}")
    logger.info(f"Fragment pairs written: {stats['fragment_pairs_written']}")
    logger.warning(
        "Partial coverage: only proteins present in the current segmented FASTAs were used"
    )
    return stats


@app.command()
def main(
    pairs_path: Annotated[Path, typer.Option(help="Stage 1 gzipped CC pair TSV.")] = DEFAULT_PAIRS,
    input_dir: Annotated[
        Path, typer.Option(help="Directory with four segmented FASTAs.")
    ] = DEFAULT_INPUT_DIR,
    output_dir: Annotated[
        Path, typer.Option(help="New directory for pair output.")
    ] = DEFAULT_OUTPUT_DIR,
    max_length: Annotated[int, typer.Option(help="Maximum segment length in amino acids.")] = 600,
    overwrite: Annotated[
        bool, typer.Option(help="Replace existing output files in output_dir.")
    ] = False,
) -> None:
    """Generate fragment combinations only for Stage 1 Same locality pairs."""
    try:
        run(pairs_path, input_dir, output_dir, max_length, overwrite)
    except (OSError, ValueError) as exc:
        logger.error(str(exc))
        raise typer.Exit(code=1) from exc


if __name__ == "__main__":
    app()
