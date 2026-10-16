"""Validate, audit, and merge isolated DPAM chunk outputs."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
import json
import os
from pathlib import Path
import re
from typing import Annotated, Any

from loguru import logger
import polars as pl
import typer

from bm_ppi.config import INTERIM_DATA_DIR

app = typer.Typer(add_completion=False)

DEFAULT_CHUNKS_PATH = INTERIM_DATA_DIR / "DPAM_runs" / "full_v1" / "chunks.tsv"
DEFAULT_INPUT_MANIFEST_PATH = (
    INTERIM_DATA_DIR / "DPAM_input" / "bm_ppi_input_manifest.tsv"
)
DEFAULT_DOMAINS_OUTPUT = INTERIM_DATA_DIR / "DPAM_input" / "bm_ppi_domains"
DOMAIN_COLUMNS = [
    "Protein",
    "Domain",
    "Range",
    "ECOD_num",
    "ECOD_key",
    "T-group",
    "DPAM_prob",
    "HH_prob",
    "DALI_zscore",
    "Hit_cov",
    "Tgroup_cov",
    "Judge",
    "Hcount",
    "Scount",
]
STEP_COMMAND_PATTERN = re.compile(r"step(?P<step>\d+)_[^\s]+\.py\s+\S+\s+(?P<protein>\S+)")


@dataclass(slots=True)
class ValidationResult:
    """Machine-readable validation state for one DPAM chunk."""

    dataset: str
    input_dir: str
    hard_success: bool = False
    status: str = "hard_failed"
    n_targets: int = 0
    n_pdb: int = 0
    n_json: int = 0
    n_domain_rows: int = 0
    n_domain_proteins: int = 0
    missing_step_logs: list[int] = field(default_factory=list)
    non_done_steps: list[int] = field(default_factory=list)
    warning_records: list[dict[str, Any]] = field(default_factory=list)
    unexpected_domain_ids: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def read_structure_list(path: Path) -> list[str]:
    """Read non-empty Protein IDs from a DPAM structure list."""
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def read_domains(path: Path) -> pl.DataFrame:
    """Read a DPAM domains table while preserving every column as text."""
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if header != DOMAIN_COLUMNS:
        raise ValueError(
            f"Unexpected domains header in {path}: expected {DOMAIN_COLUMNS}, got {header}"
        )
    return pl.read_csv(
        path,
        separator="\t",
        schema_overrides={column: pl.Utf8 for column in DOMAIN_COLUMNS},
    )


def parse_step_warning(step: int, text: str) -> list[dict[str, Any]]:
    """Extract per-protein failures from a non-done DPAM step log."""
    records: list[dict[str, Any]] = []
    for match in STEP_COMMAND_PATTERN.finditer(text):
        records.append(
            {
                "step": step,
                "protein_id": match.group("protein"),
                "message": " ".join(text.splitlines()[:10]),
            }
        )
    if not records:
        records.append(
            {
                "step": step,
                "protein_id": "",
                "message": " ".join(text.splitlines()[:10]),
            }
        )
    return records


def validate_chunk(input_dir: Path, dataset: str) -> ValidationResult:
    """Validate one DPAM chunk independently of the container exit code."""
    result = ValidationResult(dataset=dataset, input_dir=str(input_dir.resolve()))
    dataset_dir = input_dir / dataset
    structure_list = input_dir / f"{dataset}_struc.list"
    domains_path = input_dir / f"{dataset}_domains"

    if not dataset_dir.is_dir():
        result.errors.append(f"dataset directory missing: {dataset_dir}")
    if not structure_list.is_file():
        result.errors.append(f"structure list missing: {structure_list}")
    if result.errors:
        return result

    target_ids = read_structure_list(structure_list)
    result.n_targets = len(target_ids)
    duplicates = sorted(
        protein_id
        for protein_id, count in Counter(target_ids).items()
        if count > 1
    )
    if not target_ids:
        result.errors.append("structure list is empty")
    if duplicates:
        result.errors.append(f"duplicate structure-list IDs: {duplicates[:10]}")

    target_set = set(target_ids)
    pdb_ids = {path.stem for path in dataset_dir.glob("*.pdb")}
    cif_ids = {path.stem for path in dataset_dir.glob("*.cif")}
    json_ids = {path.stem for path in dataset_dir.glob("*.json")}
    structure_ids = pdb_ids | cif_ids
    result.n_pdb = len(structure_ids)
    result.n_json = len(json_ids)
    if structure_ids != target_set:
        result.errors.append(
            "structure IDs do not match struc.list: "
            f"missing={sorted(target_set - structure_ids)[:10]}, "
            f"unexpected={sorted(structure_ids - target_set)[:10]}"
        )
    if json_ids != target_set:
        result.errors.append(
            "JSON IDs do not match struc.list: "
            f"missing={sorted(target_set - json_ids)[:10]}, "
            f"unexpected={sorted(json_ids - target_set)[:10]}"
        )

    for step in range(1, 25):
        log_path = input_dir / f"{dataset}_step{step}.log"
        if not log_path.is_file():
            result.missing_step_logs.append(step)
            continue
        text = log_path.read_text(errors="replace").strip()
        if text != "done":
            result.non_done_steps.append(step)
            result.warning_records.extend(parse_step_warning(step, text))

    if result.missing_step_logs:
        result.errors.append(f"missing step logs: {result.missing_step_logs}")
    step24 = input_dir / f"{dataset}_step24.log"
    if not step24.is_file() or step24.read_text(errors="replace").strip() != "done":
        result.errors.append("step24 is not marked done")

    if not domains_path.is_file() or domains_path.stat().st_size == 0:
        result.errors.append(f"domains output missing or empty: {domains_path}")
    else:
        try:
            domains = read_domains(domains_path)
        except (OSError, ValueError) as error:
            result.errors.append(str(error))
        else:
            result.n_domain_rows = domains.height
            domain_ids = set(domains.get_column("Protein").to_list())
            result.n_domain_proteins = len(domain_ids)
            result.unexpected_domain_ids = sorted(domain_ids - target_set)
            if result.unexpected_domain_ids:
                result.errors.append(
                    f"domain output contains IDs outside this chunk: "
                    f"{result.unexpected_domain_ids[:10]}"
                )

    if result.non_done_steps:
        result.warnings.append(
            f"non-done per-step logs: {result.non_done_steps}; "
            "per-protein failures retained for audit"
        )
    result.hard_success = not result.errors
    if result.hard_success and result.warnings:
        result.status = "complete_with_warnings"
    elif result.hard_success:
        result.status = "complete"
    return result


def write_validation(result: ValidationResult, path: Path) -> None:
    """Atomically write one chunk validation JSON."""
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(asdict(result), indent=2) + "\n")
    os.replace(temporary, path)


def read_chunks(path: Path) -> pl.DataFrame:
    """Read and validate the chunk manifest schema."""
    dataframe = pl.read_csv(path, separator="\t")
    required = {
        "array_index",
        "chunk_id",
        "dataset",
        "input_dir",
        "n_targets",
        "total_aa",
        "max_length",
        "total_weight",
    }
    missing = required - set(dataframe.columns)
    if missing:
        raise ValueError(f"Chunk manifest is missing columns: {sorted(missing)}")
    duplicate_indices = (
        dataframe.group_by("array_index")
        .len()
        .filter(pl.col("len") > 1)
        .get_column("array_index")
        .to_list()
    )
    if duplicate_indices:
        raise ValueError(f"Duplicate array indices: {duplicate_indices}")
    return dataframe.sort("array_index")


def build_chunk_summary(
    chunks: pl.DataFrame,
    validations: dict[int, ValidationResult],
) -> pl.DataFrame:
    """Build one audit row per expected chunk."""
    rows: list[dict[str, Any]] = []
    for row in chunks.iter_rows(named=True):
        index = int(row["array_index"])
        result = validations[index]
        rows.append(
            {
                **row,
                "validation_status": result.status,
                "hard_success": result.hard_success,
                "n_domain_rows": result.n_domain_rows,
                "n_domain_proteins": result.n_domain_proteins,
                "missing_step_logs": ",".join(map(str, result.missing_step_logs)),
                "non_done_steps": ",".join(map(str, result.non_done_steps)),
                "n_warning_records": len(result.warning_records),
                "errors": " | ".join(result.errors),
                "warnings": " | ".join(result.warnings),
            }
        )
    return pl.DataFrame(rows)


def merge_domain_tables(
    chunks: pl.DataFrame,
    input_manifest: pl.DataFrame,
) -> pl.DataFrame:
    """Merge validated chunk tables and sort them by reference FASTA order."""
    tables: list[pl.DataFrame] = []
    protein_to_chunk: dict[str, int] = {}

    for row in chunks.iter_rows(named=True):
        index = int(row["array_index"])
        input_dir = Path(row["input_dir"])
        dataset = str(row["dataset"])
        target_ids = set(read_structure_list(input_dir / f"{dataset}_struc.list"))
        overlap = set(protein_to_chunk) & target_ids
        if overlap:
            raise ValueError(
                f"Protein IDs assigned to multiple completed chunks: {sorted(overlap)[:10]}"
            )
        protein_to_chunk.update({protein_id: index for protein_id in target_ids})

        domains = read_domains(input_dir / f"{dataset}_domains")
        unexpected = set(domains.get_column("Protein").to_list()) - target_ids
        if unexpected:
            raise ValueError(
                f"Chunk {index} contains unexpected domain IDs: {sorted(unexpected)[:10]}"
            )
        tables.append(domains)

    merged = pl.concat(tables, how="vertical")
    duplicate_rows = merged.group_by(DOMAIN_COLUMNS).len().filter(pl.col("len") > 1)
    if duplicate_rows.height:
        raise ValueError(
            f"Merged domains contain {duplicate_rows.height} exact duplicate rows"
        )

    order = input_manifest.select(
        pl.col("protein_id").cast(pl.Utf8).alias("Protein"),
        pl.col("fasta_order").cast(pl.Int64),
    )
    merged = merged.join(order, on="Protein", how="left")
    if merged.get_column("fasta_order").null_count():
        unknown = (
            merged.filter(pl.col("fasta_order").is_null())
            .get_column("Protein")
            .unique()
            .to_list()
        )
        raise ValueError(f"Merged domains contain unknown Protein IDs: {unknown[:10]}")
    return (
        merged.with_columns(
            pl.col("Domain")
            .str.extract(r"(\d+)$", 1)
            .cast(pl.Int64, strict=False)
            .fill_null(0)
            .alias("_domain_order")
        )
        .sort(["fasta_order", "_domain_order", "Domain"])
        .select(DOMAIN_COLUMNS)
    )


def build_protein_status(
    input_manifest: pl.DataFrame,
    domains: pl.DataFrame,
    validations: dict[int, ValidationResult],
) -> pl.DataFrame:
    """Classify every reference protein after DPAM chunk processing."""
    domain_counts = Counter(domains.get_column("Protein").to_list())
    warnings_by_protein: dict[str, set[int]] = defaultdict(set)
    for result in validations.values():
        for warning in result.warning_records:
            protein_id = str(warning.get("protein_id", ""))
            if protein_id:
                warnings_by_protein[protein_id].add(int(warning["step"]))

    rows: list[dict[str, Any]] = []
    for row in input_manifest.sort("fasta_order").iter_rows(named=True):
        protein_id = str(row["protein_id"])
        input_status = str(row["status"])
        n_domains = domain_counts.get(protein_id, 0)
        warning_steps = sorted(warnings_by_protein.get(protein_id, set()))
        if input_status != "prepared":
            dpam_status = "input_unavailable"
        elif warning_steps:
            dpam_status = "dpam_warning"
        elif n_domains:
            dpam_status = "domain_output"
        else:
            dpam_status = "processed_no_domain"
        rows.append(
            {
                "fasta_order": int(row["fasta_order"]),
                "protein_id": protein_id,
                "length": int(row["length"]),
                "input_status": input_status,
                "input_reason": str(row.get("reason") or ""),
                "dpam_status": dpam_status,
                "n_domains": n_domains,
                "warning_steps": ",".join(map(str, warning_steps)),
            }
        )
    return pl.DataFrame(rows)


def atomic_write_domains(dataframe: pl.DataFrame, output_path: Path) -> None:
    """Atomically write the final merged DPAM domains table."""
    temporary = output_path.with_name(f".{output_path.name}.tmp")
    dataframe.write_csv(temporary, separator="\t")
    os.replace(temporary, output_path)


@app.command("validate-chunk")
def validate_chunk_command(
    input_dir: Annotated[
        Path,
        typer.Option(help="Isolated DPAM chunk root."),
    ],
    dataset: Annotated[
        str,
        typer.Option(help="DPAM dataset name inside the chunk root."),
    ],
    validation_path: Annotated[
        Path | None,
        typer.Option(
            help="Output JSON path; defaults to <input-dir>/validation.json."
        ),
    ] = None,
) -> None:
    """Validate one DPAM chunk and return nonzero only for hard failures."""
    if validation_path is None:
        validation_path = input_dir / "validation.json"
    result = validate_chunk(input_dir, dataset)
    write_validation(result, validation_path)
    logger.info(
        f"{dataset}: status={result.status}, targets={result.n_targets}, "
        f"domain_rows={result.n_domain_rows}, "
        f"domain_proteins={result.n_domain_proteins}"
    )
    if result.warnings:
        for warning in result.warnings:
            logger.warning(warning)
    if not result.hard_success:
        for error in result.errors:
            logger.error(error)
        raise typer.Exit(code=1)


@app.command("merge")
def merge_command(
    chunks_path: Annotated[
        Path,
        typer.Option(help="Chunk manifest produced by prepare_dpam_chunks."),
    ] = DEFAULT_CHUNKS_PATH,
    input_manifest_path: Annotated[
        Path,
        typer.Option(
            help="Canonical DPAM input manifest covering the full reference FASTA."
        ),
    ] = DEFAULT_INPUT_MANIFEST_PATH,
    output_path: Annotated[
        Path,
        typer.Option(help="Final merged DPAM domains table."),
    ] = DEFAULT_DOMAINS_OUTPUT,
    report_dir: Annotated[
        Path | None,
        typer.Option(
            help="Audit output directory; defaults to the chunk-manifest directory."
        ),
    ] = None,
    force: Annotated[
        bool,
        typer.Option(
            "--force",
            help="Replace an existing merged domains file.",
        ),
    ] = False,
) -> None:
    """Audit all expected chunks and merge only a complete hard-success set."""
    try:
        if not chunks_path.is_file():
            raise ValueError(f"Chunk manifest not found: {chunks_path}")
        if not input_manifest_path.is_file():
            raise ValueError(f"Input manifest not found: {input_manifest_path}")
        if output_path.exists() and not force:
            raise FileExistsError(
                f"Merged domains output already exists; use --force to replace: "
                f"{output_path}"
            )
        if report_dir is None:
            report_dir = chunks_path.parent
        report_dir.mkdir(parents=True, exist_ok=True)

        chunks = read_chunks(chunks_path)
        input_manifest = pl.read_csv(
            input_manifest_path,
            separator="\t",
            schema_overrides={"protein_id": pl.Utf8},
        )
        validations: dict[int, ValidationResult] = {}
        for row in chunks.iter_rows(named=True):
            index = int(row["array_index"])
            input_dir = Path(row["input_dir"])
            result = validate_chunk(input_dir, str(row["dataset"]))
            validations[index] = result
            write_validation(result, input_dir / "validation.json")

        summary = build_chunk_summary(chunks, validations)
        summary.write_csv(report_dir / "dpam_chunk_summary.tsv", separator="\t")
        failed_indices = sorted(
            index for index, result in validations.items() if not result.hard_success
        )
        (report_dir / "rerun_array_indices.txt").write_text(
            "".join(f"{index}\n" for index in failed_indices)
        )
        if failed_indices:
            raise ValueError(
                f"Refusing incomplete merge; failed/not-run array indices: "
                f"{failed_indices}"
            )

        domains = merge_domain_tables(chunks, input_manifest)
        protein_status = build_protein_status(input_manifest, domains, validations)
        protein_status.write_csv(
            report_dir / "dpam_protein_status.tsv",
            separator="\t",
        )
        protein_status.filter(
            pl.col("dpam_status").is_in(["input_unavailable", "dpam_warning"])
        ).write_csv(report_dir / "dpam_failures.tsv", separator="\t")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        atomic_write_domains(domains, output_path)
        judge_counts = domains.group_by("Judge").len().sort("Judge")
        logger.info(f"Judge distribution:\n{judge_counts}")
        logger.success(
            f"Merged {domains.height} domain rows covering "
            f"{domains.get_column('Protein').n_unique()} proteins into {output_path}"
        )
    except (FileExistsError, OSError, ValueError) as error:
        logger.error(str(error))
        raise typer.Exit(code=1) from error


if __name__ == "__main__":
    app()
