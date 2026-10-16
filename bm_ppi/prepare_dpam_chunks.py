"""Create deterministic, isolated DPAM chunks from a canonical input manifest."""

from __future__ import annotations

from dataclasses import dataclass, field
import heapq
import json
import os
from pathlib import Path
import shutil
from typing import Annotated, Any

from loguru import logger
import polars as pl
from tqdm import tqdm
import typer

from bm_ppi.config import INTERIM_DATA_DIR

app = typer.Typer(add_completion=False)

DEFAULT_MANIFEST_PATH = INTERIM_DATA_DIR / "DPAM_input" / "bm_ppi_input_manifest.tsv"
DEFAULT_RUN_ROOT = INTERIM_DATA_DIR / "DPAM_runs" / "full_v1"
CHUNK_COLUMNS = [
    "array_index",
    "chunk_id",
    "dataset",
    "input_dir",
    "structure_list",
    "n_targets",
    "total_aa",
    "max_length",
    "total_weight",
]
TARGET_COLUMNS = [
    "array_index",
    "chunk_id",
    "dataset",
    "fasta_order",
    "protein_id",
    "length",
    "weight",
]


@dataclass(frozen=True, slots=True)
class Target:
    """Canonical DPAM target and its estimated processing weight."""

    fasta_order: int
    protein_id: str
    length: int
    weight: int
    pdb_path: Path
    json_path: Path


@dataclass(slots=True)
class Chunk:
    """One balanced collection of DPAM targets."""

    array_index: int
    targets: list[Target] = field(default_factory=list)
    total_weight: int = 0
    total_aa: int = 0

    @property
    def chunk_id(self) -> int:
        """Return a human-facing one-based chunk number."""
        return self.array_index + 1

    @property
    def dataset(self) -> str:
        """Return the unique DPAM dataset name for this chunk."""
        return f"bm_ppi_c{self.chunk_id:04d}"

    @property
    def max_length(self) -> int:
        """Return the maximum target length in this chunk."""
        return max((target.length for target in self.targets), default=0)

    def add(self, target: Target) -> None:
        """Add one target and update aggregate load statistics."""
        self.targets.append(target)
        self.total_weight += target.weight
        self.total_aa += target.length


def target_weight(length: int, baseline_length: int = 500) -> int:
    """Estimate DPAM cost using fixed overhead plus a quadratic length term."""
    return baseline_length**2 + length**2


def load_targets(manifest_path: Path) -> list[Target]:
    """Load prepared targets and validate canonical input paths."""
    dataframe = pl.read_csv(manifest_path, separator="\t")
    required = {
        "fasta_order",
        "protein_id",
        "length",
        "status",
        "dpam_pdb",
        "dpam_json",
    }
    missing = required - set(dataframe.columns)
    if missing:
        raise ValueError(f"Input manifest is missing columns: {sorted(missing)}")

    prepared = dataframe.filter(pl.col("status") == "prepared").sort("fasta_order")
    if prepared.is_empty():
        raise ValueError(f"No prepared targets found in manifest: {manifest_path}")
    duplicate_ids = (
        prepared.group_by("protein_id")
        .len()
        .filter(pl.col("len") > 1)
        .get_column("protein_id")
        .to_list()
    )
    if duplicate_ids:
        raise ValueError(f"Duplicate prepared Protein IDs: {duplicate_ids[:5]}")

    targets: list[Target] = []
    for row in prepared.iter_rows(named=True):
        pdb_path = Path(row["dpam_pdb"])
        json_path = Path(row["dpam_json"])
        if not pdb_path.is_file() or not json_path.is_file():
            raise ValueError(
                f"Canonical DPAM files missing for {row['protein_id']}: "
                f"PDB={pdb_path.is_file()}, JSON={json_path.is_file()}"
            )
        length = int(row["length"])
        targets.append(
            Target(
                fasta_order=int(row["fasta_order"]),
                protein_id=str(row["protein_id"]),
                length=length,
                weight=target_weight(length),
                pdb_path=pdb_path,
                json_path=json_path,
            )
        )
    return targets


def assign_chunks(targets: list[Target], n_chunks: int) -> list[Chunk]:
    """Assign targets by largest-processing-time scheduling."""
    if n_chunks < 1:
        raise ValueError("n_chunks must be at least 1")
    if n_chunks > len(targets):
        raise ValueError("n_chunks cannot exceed the number of prepared targets")

    chunks = [Chunk(array_index=index) for index in range(n_chunks)]
    heap = [(0, index) for index in range(n_chunks)]
    heapq.heapify(heap)

    ordered = sorted(
        targets,
        key=lambda target: (-target.weight, target.fasta_order, target.protein_id),
    )
    for target in ordered:
        _, index = heapq.heappop(heap)
        chunks[index].add(target)
        heapq.heappush(heap, (chunks[index].total_weight, index))
    return chunks


def create_hardlink(source: Path, destination: Path) -> None:
    """Create and verify a same-filesystem hardlink."""
    os.link(source, destination)
    source_stat = source.stat()
    destination_stat = destination.stat()
    if (
        source_stat.st_dev != destination_stat.st_dev
        or source_stat.st_ino != destination_stat.st_ino
    ):
        raise OSError(f"Hardlink verification failed: {source} -> {destination}")


def prepare_chunk(
    chunk: Chunk,
    staging_root: Path,
    final_root: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Create one isolated chunk root and return summary/detail rows."""
    chunk_name = f"chunk_{chunk.chunk_id:04d}"
    staging_chunk_root = staging_root / chunk_name
    dataset_dir = staging_chunk_root / chunk.dataset
    dataset_dir.mkdir(parents=True)

    sorted_targets = sorted(chunk.targets, key=lambda target: target.fasta_order)
    for target in sorted_targets:
        create_hardlink(
            target.pdb_path,
            dataset_dir / f"{target.protein_id}.pdb",
        )
        create_hardlink(
            target.json_path,
            dataset_dir / f"{target.protein_id}.json",
        )

    structure_list_name = f"{chunk.dataset}_struc.list"
    (staging_chunk_root / structure_list_name).write_text(
        "".join(f"{target.protein_id}\n" for target in sorted_targets)
    )

    final_chunk_root = final_root / chunk_name
    summary = {
        "array_index": chunk.array_index,
        "chunk_id": chunk.chunk_id,
        "dataset": chunk.dataset,
        "input_dir": str(final_chunk_root.resolve()),
        "structure_list": str(
            (final_chunk_root / structure_list_name).resolve()
        ),
        "n_targets": len(sorted_targets),
        "total_aa": chunk.total_aa,
        "max_length": chunk.max_length,
        "total_weight": chunk.total_weight,
    }
    details = [
        {
            "array_index": chunk.array_index,
            "chunk_id": chunk.chunk_id,
            "dataset": chunk.dataset,
            "fasta_order": target.fasta_order,
            "protein_id": target.protein_id,
            "length": target.length,
            "weight": target.weight,
        }
        for target in sorted_targets
    ]
    return summary, details


def validate_assignments(chunks: list[Chunk], targets: list[Target]) -> None:
    """Verify that chunk assignments form a disjoint complete partition."""
    expected = {target.protein_id for target in targets}
    assigned = [target.protein_id for chunk in chunks for target in chunk.targets]
    if len(assigned) != len(set(assigned)):
        counts = pl.Series("protein_id", assigned).value_counts()
        duplicates = counts.filter(pl.col("count") > 1)["protein_id"].to_list()
        raise ValueError(f"Protein IDs assigned to multiple chunks: {duplicates[:5]}")
    assigned_set = set(assigned)
    if assigned_set != expected:
        missing = sorted(expected - assigned_set)
        unexpected = sorted(assigned_set - expected)
        raise ValueError(
            f"Chunk partition mismatch: missing={missing[:5]}, "
            f"unexpected={unexpected[:5]}"
        )


def ensure_run_root_available(run_root: Path, force: bool) -> None:
    """Refuse implicit overwrite of a chunk run root."""
    if run_root.exists() and not force:
        raise FileExistsError(f"DPAM run root already exists: {run_root}")
    if run_root.exists():
        shutil.rmtree(run_root)


@app.command()
def main(
    manifest_path: Annotated[
        Path,
        typer.Option(
            help="Canonical DPAM input manifest produced by prepare_dpam_input."
        ),
    ] = DEFAULT_MANIFEST_PATH,
    run_root: Annotated[
        Path,
        typer.Option(help="Versioned root for all isolated DPAM chunk runs."),
    ] = DEFAULT_RUN_ROOT,
    n_chunks: Annotated[
        int,
        typer.Option(min=1, help="Number of balanced chunks."),
    ] = 64,
    force: Annotated[
        bool,
        typer.Option(
            "--force",
            help="Replace an existing run root and all of its chunk data.",
        ),
    ] = False,
) -> None:
    """Build deterministic DPAM chunk datasets without running DPAM."""
    try:
        if not manifest_path.is_file():
            raise ValueError(f"Input manifest not found: {manifest_path}")
        staging_root = run_root.parent / f".{run_root.name}.staging"
        if staging_root.exists():
            raise FileExistsError(
                f"Staging run root already exists; inspect before retrying: "
                f"{staging_root}"
            )
        ensure_run_root_available(run_root, force)
        staging_root.mkdir(parents=True)

        targets = load_targets(manifest_path)
        chunks = assign_chunks(targets, n_chunks)
        validate_assignments(chunks, targets)

        chunk_rows: list[dict[str, Any]] = []
        target_rows: list[dict[str, Any]] = []
        for chunk in tqdm(chunks, desc="Creating DPAM chunks"):
            summary, details = prepare_chunk(chunk, staging_root, run_root)
            chunk_rows.append(summary)
            target_rows.extend(details)

        pl.DataFrame(chunk_rows).select(CHUNK_COLUMNS).write_csv(
            staging_root / "chunks.tsv",
            separator="\t",
        )
        pl.DataFrame(target_rows).select(TARGET_COLUMNS).write_csv(
            staging_root / "chunk_targets.tsv",
            separator="\t",
        )
        run_metadata = {
            "source_manifest": str(manifest_path.resolve()),
            "n_targets": len(targets),
            "n_chunks": n_chunks,
            "weight_formula": "500^2 + sequence_length^2",
            "total_aa": sum(target.length for target in targets),
            "total_weight": sum(target.weight for target in targets),
        }
        (staging_root / "run_manifest.json").write_text(
            json.dumps(run_metadata, indent=2) + "\n"
        )
        staging_root.rename(run_root)

        counts = [len(chunk.targets) for chunk in chunks]
        amino_acids = [chunk.total_aa for chunk in chunks]
        weights = [chunk.total_weight for chunk in chunks]
        logger.info(
            f"Chunk targets min/max={min(counts)}/{max(counts)}; "
            f"amino acids min/max={min(amino_acids)}/{max(amino_acids)}; "
            f"weight max/min={max(weights) / min(weights):.4f}"
        )
        logger.success(
            f"Prepared {len(targets)} targets in {n_chunks} isolated chunks at "
            f"{run_root}"
        )
    except (FileExistsError, OSError, ValueError) as error:
        logger.error(str(error))
        raise typer.Exit(code=1) from error


if __name__ == "__main__":
    app()
