"""Prepare validated DPAM inputs from ColabFold rank-001 outputs."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path
import re
import shutil
from typing import Annotated, Any

from Bio import SeqIO
from loguru import logger
import polars as pl
from tqdm import tqdm
import typer

from bm_ppi.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR, RAW_DATA_DIR

app = typer.Typer(add_completion=False)

BATCH_IDS = ("001", "002", "003", "004")
JOBNAME_PATTERN = re.compile(r"[^A-Za-z0-9._-]")
DATASET_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
MANIFEST_COLUMNS = [
    "fasta_order",
    "protein_id",
    "colabfold_jobname",
    "length",
    "status",
    "reason",
    "source_output_dir",
    "source_pdb",
    "source_json",
    "dpam_pdb",
    "dpam_json",
]


@dataclass(frozen=True, slots=True)
class ReferenceProtein:
    """Reference FASTA record and its deterministic ColabFold jobname."""

    fasta_order: int
    protein_id: str
    jobname: str
    length: int


@dataclass(frozen=True, slots=True)
class Rank001Pair:
    """Matched rank-001 ColabFold structure and scores files."""

    pdb_path: Path
    scores_path: Path | None


@dataclass(frozen=True, slots=True)
class PreparedScores:
    """Validated scores converted to the JSON schema required by DPAM."""

    data: list[dict[str, Any]]
    sequence_length: int


def default_output_dirs() -> list[Path]:
    """Return main and five-model ColabFold output directories."""
    main = [PROCESSED_DATA_DIR / f"colabfold_output_gcf_{batch}" for batch in BATCH_IDS]
    five_model = [
        PROCESSED_DATA_DIR / f"colabfold_output_gcf_{batch}_5m"
        for batch in BATCH_IDS
    ]
    return [*main, *five_model]


DEFAULT_OUTPUT_DIRS = default_output_dirs()
DEFAULT_REFERENCE_FASTA = (
    RAW_DATA_DIR / "GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa"
)
DEFAULT_DPAM_INPUT_DIR = INTERIM_DATA_DIR / "DPAM_input"


def sanitize_colabfold_jobname(description: str) -> str:
    """Reconstruct the filename-safe query name used by ColabFold."""
    return JOBNAME_PATTERN.sub("_", description)


def load_reference_proteins(fasta_path: Path) -> list[ReferenceProtein]:
    """Load FASTA records and validate ID/jobname uniqueness."""
    records: list[ReferenceProtein] = []
    seen_ids: set[str] = set()
    seen_jobnames: set[str] = set()

    for index, record in enumerate(SeqIO.parse(str(fasta_path), "fasta")):
        protein_id = str(record.id)
        jobname = sanitize_colabfold_jobname(record.description)
        if protein_id in seen_ids:
            raise ValueError(f"Duplicate FASTA protein ID: {protein_id}")
        if jobname in seen_jobnames:
            raise ValueError(f"Duplicate sanitized ColabFold jobname: {jobname}")
        seen_ids.add(protein_id)
        seen_jobnames.add(jobname)
        records.append(
            ReferenceProtein(
                fasta_order=index,
                protein_id=protein_id,
                jobname=jobname,
                length=len(record.seq),
            )
        )

    if not records:
        raise ValueError(f"No protein records found in FASTA: {fasta_path}")
    return records


def collect_selected_jobnames(a3m_dirs: list[Path]) -> set[str]:
    """Collect a unique jobname subset from one or more A3M directories."""
    selected: set[str] = set()
    for a3m_dir in a3m_dirs:
        if not a3m_dir.is_dir():
            raise ValueError(f"A3M directory not found: {a3m_dir}")
        for path in sorted(a3m_dir.glob("*.a3m")):
            if path.stem in selected:
                raise ValueError(f"Duplicate A3M jobname across inputs: {path.stem}")
            selected.add(path.stem)
    if not selected:
        raise ValueError("No A3M files found in the selected A3M directories")
    return selected


def rank001_quality(scores_path: Path | None) -> tuple[float, float]:
    """Return a deterministic quality key for a rank-001 candidate."""
    if scores_path is None:
        return float("-inf"), float("-inf")
    try:
        with scores_path.open() as handle:
            scores = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return float("-inf"), float("-inf")

    plddt = scores.get("plddt")
    mean_plddt = (
        sum(plddt) / len(plddt)
        if isinstance(plddt, list) and plddt
        else float("-inf")
    )
    ptm = scores.get("ptm")
    return mean_plddt, ptm if isinstance(ptm, (int, float)) else float("-inf")


def collect_model_candidates(model_dir: Path) -> dict[str, list[Rank001Pair]]:
    """Collect rank-001 candidates from one five-model subdirectory."""
    scores_by_key: dict[tuple[str, str], Path] = {}
    for scores_path in sorted(model_dir.glob("*_scores_rank_001_*.json")):
        jobname, suffix = scores_path.stem.split("_scores_rank_001_", maxsplit=1)
        scores_by_key[(jobname, suffix)] = scores_path

    candidates: dict[str, list[Rank001Pair]] = {}
    for pdb_path in sorted(model_dir.glob("*_unrelaxed_rank_001_*.pdb")):
        jobname, suffix = pdb_path.stem.split("_unrelaxed_rank_001_", maxsplit=1)
        candidates.setdefault(jobname, []).append(
            Rank001Pair(
                pdb_path=pdb_path,
                scores_path=scores_by_key.get((jobname, suffix)),
            )
        )
    return candidates


def collect_rank001_pairs(output_dirs: list[Path]) -> dict[str, Rank001Pair]:
    """Collect one unambiguous rank-001 PDB/scores pair per jobname."""
    pairs: dict[str, Rank001Pair] = {}
    for output_dir in output_dirs:
        if not output_dir.is_dir():
            raise ValueError(f"ColabFold output directory not found: {output_dir}")

        scores_by_jobname: dict[str, Path] = {}
        for scores_path in sorted(output_dir.glob("*_scores_rank_001_*.json")):
            jobname = scores_path.name.split("_scores_rank_001_", maxsplit=1)[0]
            if jobname in scores_by_jobname:
                raise ValueError(
                    f"Multiple rank-001 scores JSON files for {jobname}: "
                    f"{scores_by_jobname[jobname]} and {scores_path}"
                )
            scores_by_jobname[jobname] = scores_path

        output_pairs: dict[str, Rank001Pair] = {}
        for pdb_path in sorted(output_dir.glob("*_unrelaxed_rank_001_*.pdb")):
            jobname = pdb_path.name.split("_unrelaxed_rank_001_", maxsplit=1)[0]
            if jobname in output_pairs:
                previous = output_pairs[jobname].pdb_path
                raise ValueError(
                    f"Duplicate rank-001 PDB for {jobname}: {previous} and {pdb_path}"
                )
            output_pairs[jobname] = Rank001Pair(
                pdb_path=pdb_path,
                scores_path=scores_by_jobname.get(jobname),
            )

        model_candidates: dict[str, list[Rank001Pair]] = {}
        for model_dir in sorted(output_dir.glob("model_*")):
            if not model_dir.is_dir():
                continue
            for jobname, candidates in collect_model_candidates(model_dir).items():
                model_candidates.setdefault(jobname, []).extend(candidates)

        for jobname, candidates in model_candidates.items():
            if jobname in output_pairs:
                previous = output_pairs[jobname].pdb_path
                raise ValueError(
                    f"Multiple rank-001 sources for {jobname}: "
                    f"{previous} and {candidates[0].pdb_path}"
                )
            output_pairs[jobname] = max(
                candidates,
                key=lambda pair: (
                    *rank001_quality(pair.scores_path),
                    str(pair.pdb_path),
                ),
            )

        for jobname, pair in output_pairs.items():
            if jobname in pairs:
                raise ValueError(
                    f"Duplicate rank-001 PDB for {jobname}: "
                    f"{pairs[jobname].pdb_path} and {pair.pdb_path}"
                )
            pairs[jobname] = pair
    return pairs


def check_pdb_for_nan(pdb_path: Path) -> bool:
    """Return whether a PDB contains NaN coordinates or B-factors."""
    with pdb_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            coordinates = line[30:54]
            bfactor = line[60:66] if len(line) > 66 else ""
            if "nan" in coordinates.lower() or "nan" in bfactor.lower():
                return True
    return False


def load_and_validate_scores(
    scores_json_path: Path,
    expected_length: int,
) -> PreparedScores:
    """Validate a ColabFold scores JSON and convert it to DPAM schema."""
    with scores_json_path.open() as handle:
        scores = json.load(handle)

    pae = scores.get("pae")
    plddt = scores.get("plddt")
    max_pae = scores.get("max_pae")
    if not isinstance(pae, list) or not pae:
        raise ValueError("missing or empty pae matrix")
    if not isinstance(plddt, list) or not plddt:
        raise ValueError("missing or empty plddt vector")
    if not isinstance(max_pae, (int, float)):
        raise TypeError("missing or non-numeric max_pae")

    sequence_length = len(plddt)
    if sequence_length != expected_length:
        raise ValueError(
            f"pLDDT length {sequence_length} != FASTA length {expected_length}"
        )
    if len(pae) != sequence_length:
        raise ValueError(
            f"PAE row count {len(pae)} != pLDDT length {sequence_length}"
        )
    malformed_row = next(
        (
            index
            for index, row in enumerate(pae)
            if not isinstance(row, list) or len(row) != sequence_length
        ),
        None,
    )
    if malformed_row is not None:
        raise ValueError(
            f"PAE row {malformed_row} does not have length {sequence_length}"
        )

    data = [
        {
            "predicted_aligned_error": pae,
            "max_predicted_aligned_error": max_pae,
            "plddt": plddt,
        }
    ]
    return PreparedScores(data=data, sequence_length=sequence_length)


def write_dpam_json(scores: PreparedScores, output_path: Path) -> None:
    """Write validated scores using the JSON schema expected by DPAM."""
    with output_path.open("w") as handle:
        json.dump(scores.data, handle, separators=(",", ":"))


def write_manifest(rows: list[dict[str, Any]], path: Path) -> None:
    """Write the input manifest with a stable column order."""
    pl.DataFrame(rows).select(MANIFEST_COLUMNS).write_csv(path, separator="\t")


def validate_dataset_name(dataset: str) -> None:
    """Reject unsafe dataset names before constructing paths."""
    if not DATASET_PATTERN.fullmatch(dataset):
        raise ValueError(
            "Dataset must start with an alphanumeric character and contain only "
            "letters, numbers, dots, underscores, or hyphens"
        )


def ensure_output_is_available(
    dataset_dir: Path,
    output_files: list[Path],
    force: bool,
) -> None:
    """Refuse implicit overwrite and remove outputs only with explicit force."""
    existing = [path for path in [dataset_dir, *output_files] if path.exists()]
    if existing and not force:
        formatted = "\n  ".join(str(path) for path in existing)
        raise FileExistsError(
            f"DPAM input outputs already exist; use --force to replace:\n  {formatted}"
        )
    if not force:
        return
    if dataset_dir.exists():
        shutil.rmtree(dataset_dir)
    for path in output_files:
        if path.is_file() or path.is_symlink():
            path.unlink()


@app.command()
def main(
    reference_fasta: Annotated[
        Path,
        typer.Option(
            help="Reference FASTA defining canonical Protein IDs and expected targets."
        ),
    ] = DEFAULT_REFERENCE_FASTA,
    output_dirs: Annotated[
        list[Path],
        typer.Option(
            "--output-dir",
            help="ColabFold output directory; repeat for all main and retry outputs.",
        ),
    ] = DEFAULT_OUTPUT_DIRS,
    a3m_dirs: Annotated[
        list[Path] | None,
        typer.Option(
            "--a3m-dir",
            help="Optional A3M directories used only to select a smoke-test subset.",
        ),
    ] = None,
    dpam_input_dir: Annotated[
        Path,
        typer.Option(help="Directory in which the canonical DPAM dataset is created."),
    ] = DEFAULT_DPAM_INPUT_DIR,
    dataset: Annotated[str, typer.Option(help="DPAM dataset name.")] = "bm_ppi",
    force: Annotated[
        bool,
        typer.Option(
            "--force",
            help="Replace an existing canonical dataset and its manifests.",
        ),
    ] = False,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run",
            help="Validate every candidate and report counts without writing outputs.",
        ),
    ] = False,
) -> None:
    """Prepare a validated, canonical DPAM dataset from ColabFold outputs."""
    try:
        validate_dataset_name(dataset)
        if not reference_fasta.is_file():
            raise ValueError(f"Reference FASTA not found: {reference_fasta}")

        references = load_reference_proteins(reference_fasta)
        reference_by_jobname = {record.jobname: record for record in references}
        if a3m_dirs:
            selected_jobnames = collect_selected_jobnames(a3m_dirs)
            unknown = sorted(selected_jobnames - reference_by_jobname.keys())
            if unknown:
                raise ValueError(
                    f"{len(unknown)} A3M jobnames do not match the reference FASTA, "
                    f"e.g. {unknown[:5]}"
                )
            references = [
                record for record in references if record.jobname in selected_jobnames
            ]

        pairs = collect_rank001_pairs(output_dirs)
        logger.info(
            f"Reference targets: {len(references)}; rank-001 structures found: "
            f"{len(pairs)}"
        )

        dataset_dir = dpam_input_dir / dataset
        structure_list = dpam_input_dir / f"{dataset}_struc.list"
        failed_list = dpam_input_dir / f"{dataset}_failed.list"
        manifest_path = dpam_input_dir / f"{dataset}_input_manifest.tsv"
        staging_dir = dpam_input_dir / f".{dataset}.staging"
        output_files = [structure_list, failed_list, manifest_path]

        if not dry_run:
            dpam_input_dir.mkdir(parents=True, exist_ok=True)
            if staging_dir.exists():
                raise FileExistsError(
                    f"Staging directory already exists; inspect it before retrying: "
                    f"{staging_dir}"
                )
            ensure_output_is_available(dataset_dir, output_files, force)
            staging_dir.mkdir()

        rows: list[dict[str, Any]] = []
        prepared_ids: list[str] = []
        failed_rows: list[tuple[str, str]] = []

        for reference in tqdm(references, desc="Preparing DPAM input"):
            pair = pairs.get(reference.jobname)
            status = "prepared"
            reason = ""
            source_dir = ""
            source_pdb = ""
            source_json = ""
            final_pdb = dataset_dir / f"{reference.protein_id}.pdb"
            final_json = dataset_dir / f"{reference.protein_id}.json"

            if pair is None:
                status = "missing_rank001"
                reason = "rank-001 PDB not found"
            elif pair.scores_path is None:
                status = "missing_scores_json"
                reason = "matching rank-001 scores JSON not found"
                source_dir = str(pair.pdb_path.parent.resolve())
                source_pdb = str(pair.pdb_path.resolve())
            else:
                source_dir = str(pair.pdb_path.parent.resolve())
                source_pdb = str(pair.pdb_path.resolve())
                source_json = str(pair.scores_path.resolve())
                if check_pdb_for_nan(pair.pdb_path):
                    status = "invalid_pdb"
                    reason = "PDB contains NaN coordinates or B-factors"
                else:
                    try:
                        scores = load_and_validate_scores(
                            pair.scores_path,
                            expected_length=reference.length,
                        )
                    except (json.JSONDecodeError, OSError, TypeError, ValueError) as error:
                        status = "invalid_scores_json"
                        reason = str(error)
                    else:
                        if not dry_run:
                            shutil.copy2(
                                pair.pdb_path,
                                staging_dir / f"{reference.protein_id}.pdb",
                            )
                            write_dpam_json(
                                scores,
                                staging_dir / f"{reference.protein_id}.json",
                            )
                        prepared_ids.append(reference.protein_id)

            if status != "prepared":
                failed_rows.append((reference.protein_id, reason))

            rows.append(
                {
                    "fasta_order": reference.fasta_order,
                    "protein_id": reference.protein_id,
                    "colabfold_jobname": reference.jobname,
                    "length": reference.length,
                    "status": status,
                    "reason": reason,
                    "source_output_dir": source_dir,
                    "source_pdb": source_pdb,
                    "source_json": source_json,
                    "dpam_pdb": str(final_pdb.resolve()),
                    "dpam_json": str(final_json.resolve()),
                }
            )

        if not dry_run:
            staging_dir.rename(dataset_dir)
            structure_list.write_text("".join(f"{name}\n" for name in prepared_ids))
            write_manifest(rows, manifest_path)
            if failed_rows:
                failed_list.write_text(
                    "".join(f"{protein_id}\t{reason}\n" for protein_id, reason in failed_rows)
                )

        status_counts = Counter(row["status"] for row in rows)
        logger.info(
            "Input status: "
            + ", ".join(
                f"{status}={count}" for status, count in sorted(status_counts.items())
            )
        )
        logger.info(
            f"Summary: expected={len(rows)}, prepared={len(prepared_ids)}, "
            f"failed={len(failed_rows)}"
        )
        if dry_run:
            logger.success("DPAM input dry-run validation completed; no files written")
        else:
            logger.success(f"Canonical DPAM input prepared at: {dataset_dir}")
    except (FileExistsError, OSError, ValueError) as error:
        logger.error(str(error))
        raise typer.Exit(code=1) from error


if __name__ == "__main__":
    app()
