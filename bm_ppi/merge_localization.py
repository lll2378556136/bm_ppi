"""Merge pilot subcellular-localization predictions into one protein table."""

from __future__ import annotations

from collections import Counter
import csv
import json
from pathlib import Path
import re

from Bio import SeqIO
from loguru import logger
import typer

from bm_ppi.config import INTERIM_DATA_DIR

app = typer.Typer()

LOCALIZATION_DIR = INTERIM_DATA_DIR / "localization"
PILOT_DIR = LOCALIZATION_DIR / "pilot"
RAW_DIR = LOCALIZATION_DIR / "raw"

DEFAULT_REFERENCE_FASTA = PILOT_DIR / "protein_pilot100.faa"
DEFAULT_DEEPLOC_CSV = RAW_DIR / "deeploc_2.1/pilot/deeploc_2.1_pilot100.csv"
DEFAULT_DEEPTMHMM_3LINE = (
    RAW_DIR / "deeptmhmm_1.0/pilot/deeptmhmm_1.0_pilot100_topologies.3line"
)
DEFAULT_DEEPTMHMM_GFF3 = RAW_DIR / "deeptmhmm_1.0/pilot/deeptmhmm_1.0_pilot100_tmrs.gff3"
DEFAULT_SIGNALP_SUMMARY = RAW_DIR / "signalp_6.0/pilot/signalp_6.0_pilot100_results.txt"
DEFAULT_SIGNALP_JSON = RAW_DIR / "signalp_6.0/pilot/signalp_6.0_pilot100_output.json"
DEFAULT_OUTPUT = LOCALIZATION_DIR / "processed/protein_localization_master.tsv"

DEEPLOC_PROBABILITY_COLUMNS = {
    "Cytoplasm": "deeploc_cytoplasm_probability",
    "Nucleus": "deeploc_nucleus_probability",
    "Extracellular": "deeploc_extracellular_probability",
    "Cell membrane": "deeploc_cell_membrane_probability",
    "Mitochondrion": "deeploc_mitochondrion_probability",
    "Plastid": "deeploc_plastid_probability",
    "Endoplasmic reticulum": "deeploc_endoplasmic_reticulum_probability",
    "Lysosome/Vacuole": "deeploc_lysosome_vacuole_probability",
    "Golgi apparatus": "deeploc_golgi_apparatus_probability",
    "Peroxisome": "deeploc_peroxisome_probability",
    "Peripheral": "deeploc_peripheral_probability",
    "Transmembrane": "deeploc_transmembrane_probability",
    "Lipid anchor": "deeploc_lipid_anchor_probability",
    "Soluble": "deeploc_soluble_probability",
}


def normalize_id(value: str) -> str:
    """Return the first whitespace-delimited token from a sequence header."""
    stripped = value.strip().lstrip(">")
    if not stripped:
        return ""
    return stripped.split()[0]


def parse_reference_fasta(path: Path) -> list[dict[str, str | int]]:
    """Read reference protein IDs, headers, and lengths from a FASTA file."""
    records: list[dict[str, str | int]] = []
    seen: set[str] = set()
    for record in SeqIO.parse(path, "fasta"):
        protein_id = normalize_id(record.id)
        if not protein_id:
            raise ValueError(f"Empty protein ID in FASTA: {path}")
        if protein_id in seen:
            raise ValueError(f"Duplicate protein ID in FASTA: {protein_id}")
        seen.add(protein_id)
        records.append(
            {
                "protein_id": protein_id,
                "fasta_header": record.description,
                "length": len(record.seq),
            }
        )
    if not records:
        raise ValueError(f"No FASTA records found: {path}")
    return records


def format_probability(value: str, column: str, protein_id: str) -> str:
    """Validate and format a probability value from a tool result."""
    if value == "":
        return ""
    try:
        return f"{float(value):.6f}"
    except ValueError as exc:
        raise ValueError(f"Invalid probability in {column} for {protein_id}: {value}") from exc


def parse_deeploc(path: Path) -> dict[str, dict[str, str]]:
    """Parse the DeepLoc CSV keyed by normalized protein ID."""
    required_columns = {"Protein_ID", "Localizations", "Signals", "Membrane types"}
    parsed: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not required_columns.issubset(reader.fieldnames):
            missing = sorted(required_columns - set(reader.fieldnames or []))
            raise ValueError(f"DeepLoc CSV is missing columns: {missing}")
        for row in reader:
            protein_id = normalize_id(row.get("Protein_ID", ""))
            if not protein_id:
                raise ValueError("DeepLoc CSV contains an empty Protein_ID")
            if protein_id in parsed:
                raise ValueError(f"Duplicate protein ID in DeepLoc CSV: {protein_id}")
            parsed_row = {
                "deeploc_localizations": row.get("Localizations", ""),
                "deeploc_signals": row.get("Signals", ""),
                "deeploc_membrane_types": row.get("Membrane types", ""),
            }
            for source_column, output_column in DEEPLOC_PROBABILITY_COLUMNS.items():
                parsed_row[output_column] = format_probability(
                    row.get(source_column, ""), source_column, protein_id
                )
            parsed[protein_id] = parsed_row
    return parsed


def find_regions(topology: str, symbols: set[str]) -> list[str]:
    """Convert contiguous topology symbols into one-based coordinate ranges."""
    regions: list[str] = []
    start: int | None = None
    for index, symbol in enumerate(topology, start=1):
        if symbol in symbols and start is None:
            start = index
        elif symbol not in symbols and start is not None:
            regions.append(f"{start}-{index - 1}")
            start = None
    if start is not None:
        regions.append(f"{start}-{len(topology)}")
    return regions


def parse_deeptmhmm_3line(path: Path) -> dict[str, dict[str, str | int]]:
    """Parse DeepTMHMM's three-line topology output."""
    lines = [line.rstrip("\n") for line in path.open(encoding="utf-8") if line.strip()]
    if len(lines) % 3 != 0:
        raise ValueError(f"DeepTMHMM 3line file does not contain complete triplets: {path}")

    parsed: dict[str, dict[str, str | int]] = {}
    for index in range(0, len(lines), 3):
        header, sequence, topology = lines[index : index + 3]
        if not header.startswith(">") or " | " not in header:
            raise ValueError(f"Malformed DeepTMHMM header near line {index + 1}: {header}")
        protein_id = normalize_id(header.split(" | ", 1)[0])
        protein_type = header.split(" | ", 1)[1].strip()
        if not protein_id or not sequence or not topology:
            raise ValueError(f"Incomplete DeepTMHMM record for: {protein_id}")
        if len(sequence) != len(topology):
            raise ValueError(f"Sequence/topology length mismatch for: {protein_id}")
        if protein_id in parsed:
            raise ValueError(f"Duplicate protein ID in DeepTMHMM 3line: {protein_id}")
        tm_regions = find_regions(topology, {"M", "B"})
        signal_regions = find_regions(topology, {"S"})
        parsed[protein_id] = {
            "deeptmhmm_protein_type": protein_type,
            "deeptmhmm_sequence_length": len(sequence),
            "deeptmhmm_topology": topology,
            "deeptmhmm_tm_regions_3line": ";".join(tm_regions),
            "deeptmhmm_signal_regions_3line": ";".join(signal_regions),
            "deeptmhmm_n_tm_regions_3line": len(tm_regions),
            "deeptmhmm_has_signal_region": bool(signal_regions),
        }
    return parsed


def parse_deeptmhmm_gff3(path: Path) -> dict[str, dict[str, str | int]]:
    """Parse DeepTMHMM GFF3-like records and extract feature coordinates."""
    records: dict[str, dict[str, str | int]] = {}
    current: dict[str, str | int | list[str]] | None = None

    def flush() -> None:
        if current is None:
            return
        protein_id = str(current["protein_id"])
        if protein_id in records:
            raise ValueError(f"Duplicate protein ID in DeepTMHMM GFF3: {protein_id}")
        records[protein_id] = {
            "deeptmhmm_gff3_length": int(current["deeptmhmm_gff3_length"]),
            "deeptmhmm_gff3_reported_tmr_count": int(
                current["deeptmhmm_gff3_reported_tmr_count"]
            ),
            "deeptmhmm_gff3_feature_count": int(current["deeptmhmm_gff3_feature_count"]),
            "deeptmhmm_gff3_tm_regions": ";".join(
                current["deeptmhmm_gff3_tm_regions"]
            ),
            "deeptmhmm_gff3_signal_regions": ";".join(
                current["deeptmhmm_gff3_signal_regions"]
            ),
        }

    length_pattern = re.compile(r"^#\s+(\S+)\s+Length:\s+(\d+)")
    tmr_pattern = re.compile(r"^#\s+(\S+)\s+Number of predicted TMRs:\s+(\d+)")

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            length_match = length_pattern.match(line)
            if length_match:
                if current is not None:
                    flush()
                current = {
                    "protein_id": length_match.group(1),
                    "deeptmhmm_gff3_length": int(length_match.group(2)),
                    "deeptmhmm_gff3_reported_tmr_count": 0,
                    "deeptmhmm_gff3_feature_count": 0,
                    "deeptmhmm_gff3_tm_regions": [],
                    "deeptmhmm_gff3_signal_regions": [],
                }
                continue
            tmr_match = tmr_pattern.match(line)
            if tmr_match:
                if current is None or str(current["protein_id"]) != tmr_match.group(1):
                    raise ValueError(f"Malformed DeepTMHMM GFF3 record: {line}")
                current["deeptmhmm_gff3_reported_tmr_count"] = int(tmr_match.group(2))
                continue
            if line == "//":
                flush()
                current = None
                continue
            if not line or line.startswith("#") or current is None:
                continue

            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(f"Malformed DeepTMHMM GFF3 feature line: {line}")
            protein_id, feature_type = fields[0], fields[1]
            if protein_id != str(current["protein_id"]):
                raise ValueError(f"GFF3 feature ID mismatch: {line}")
            try:
                start, end = int(fields[2]), int(fields[3])
            except ValueError as exc:
                raise ValueError(f"Invalid GFF3 coordinates: {line}") from exc
            region = f"{start}-{end}"
            current["deeptmhmm_gff3_feature_count"] = (
                int(current["deeptmhmm_gff3_feature_count"]) + 1
            )
            if feature_type == "TMhelix":
                current["deeptmhmm_gff3_tm_regions"].append(region)
            elif feature_type == "signal":
                current["deeptmhmm_gff3_signal_regions"].append(region)

    if current is not None:
        flush()
    return records


def parse_signalp_summary(path: Path) -> dict[str, dict[str, str]]:
    """Parse SignalP's one-line-per-protein summary."""
    parsed: dict[str, dict[str, str]] = {}
    cleavage_pattern = re.compile(r"CS pos:\s*([0-9]+-[0-9]+)")
    probability_pattern = re.compile(r"Pr:\s*([0-9.]+)")

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                raise ValueError(f"Malformed SignalP summary line: {line}")
            protein_id = normalize_id(fields[0])
            if not protein_id:
                raise ValueError("SignalP summary contains an empty protein ID")
            if protein_id in parsed:
                raise ValueError(f"Duplicate protein ID in SignalP summary: {protein_id}")
            cleavage_match = cleavage_pattern.search(fields[4])
            probability_match = probability_pattern.search(fields[4])
            parsed[protein_id] = {
                "signalp_prediction": fields[1],
                "signalp_other_probability": fields[2],
                "signalp_sp_probability": fields[3],
                "signalp_cleavage_position": cleavage_match.group(1)
                if cleavage_match
                else "",
                "signalp_cleavage_probability": probability_match.group(1)
                if probability_match
                else "",
            }
    return parsed


def parse_signalp_json(path: Path) -> tuple[set[str], int, int]:
    """Read SignalP JSON metadata and sequence IDs for an independent check."""
    with path.open(encoding="utf-8") as handle:
        data = json.load(handle)
    info = data.get("INFO", {})
    failed_jobs = int(info.get("failedjobs", 0))
    reported_size = int(info.get("size", 0))
    sequence_ids = {
        normalize_id(sequence_name) for sequence_name in data.get("SEQUENCES", {})
    }
    if failed_jobs:
        raise ValueError(f"SignalP reports failed jobs: {failed_jobs}")
    return sequence_ids, reported_size, failed_jobs


def validate_ids(name: str, ids: list[str] | set[str], reference_ids: set[str]) -> None:
    """Require an exact, unique ID match against the reference FASTA."""
    id_list = list(ids)
    duplicates = sorted(
        protein_id for protein_id, count in Counter(id_list).items() if count > 1
    )
    observed = set(id_list)
    missing = sorted(reference_ids - observed)
    extra = sorted(observed - reference_ids)
    if duplicates or missing or extra:
        raise ValueError(
            f"{name} ID mismatch: duplicates={duplicates[:5]}, "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )


def split_labels(value: str) -> set[str]:
    """Split tool labels while treating an empty value as no label."""
    return {label for label in value.split("|") if label}


def build_rows(
    reference_records: list[dict[str, str | int]],
    deeploc: dict[str, dict[str, str]],
    deeptmhmm_3line: dict[str, dict[str, str | int]],
    deeptmhmm_gff3: dict[str, dict[str, str | int]],
    signalp: dict[str, dict[str, str]],
) -> list[dict[str, str | int | bool]]:
    """Join all parsed results in the reference FASTA order."""
    rows: list[dict[str, str | int | bool]] = []
    for reference in reference_records:
        protein_id = str(reference["protein_id"])
        dl = deeploc[protein_id]
        tm = deeptmhmm_3line[protein_id]
        gff = deeptmhmm_gff3[protein_id]
        sp = signalp[protein_id]

        deeploc_labels = split_labels(dl["deeploc_localizations"])
        deeploc_membrane_types = split_labels(dl["deeploc_membrane_types"])
        tm_regions = str(tm["deeptmhmm_tm_regions_3line"])
        gff_tm_regions = str(gff["deeptmhmm_gff3_tm_regions"])
        is_transmembrane = bool(tm_regions or gff_tm_regions) or "Transmembrane" in (
            deeploc_membrane_types
        )
        is_membrane_associated = bool(
            deeploc_membrane_types - {"Soluble"}
        ) or is_transmembrane
        is_secretory_pathway = (
            sp["signalp_prediction"] != "OTHER"
            or bool(tm["deeptmhmm_signal_regions_3line"])
            or "Signal peptide" in split_labels(dl["deeploc_signals"])
        )
        if len(deeploc_labels) > 1:
            localization_status = "multi_localized"
        elif deeploc_labels:
            localization_status = "single_localized"
        else:
            localization_status = "unknown"

        row: dict[str, str | int | bool] = {
            "protein_id": protein_id,
            "fasta_header": reference["fasta_header"],
            "length": reference["length"],
            **dl,
            **tm,
            **gff,
            **sp,
            "is_multilocalized": len(deeploc_labels) > 1,
            "is_transmembrane": is_transmembrane,
            "is_membrane_associated": is_membrane_associated,
            "is_secretory_pathway": is_secretory_pathway,
            "localization_status": localization_status,
        }
        if int(tm["deeptmhmm_sequence_length"]) != int(reference["length"]):
            raise ValueError(f"DeepTMHMM sequence length mismatch for {protein_id}")
        if int(gff["deeptmhmm_gff3_length"]) != int(reference["length"]):
            raise ValueError(f"DeepTMHMM GFF3 length mismatch for {protein_id}")
        gff_tm_region_text = str(gff["deeptmhmm_gff3_tm_regions"])
        gff_tm_region_count = (
            len(gff_tm_region_text.split(";")) if gff_tm_region_text else 0
        )
        if int(gff["deeptmhmm_gff3_reported_tmr_count"]) != gff_tm_region_count:
            logger.warning(f"DeepTMHMM GFF3 TMR count differs for {protein_id}")
        rows.append(row)
    return rows


def write_tsv(path: Path, rows: list[dict[str, str | int | bool]]) -> None:
    """Write merged rows as a tab-separated table."""
    if not rows:
        raise ValueError("No rows to write")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@app.command()
def main(
    reference_fasta: Path = DEFAULT_REFERENCE_FASTA,
    deeploc_csv: Path = DEFAULT_DEEPLOC_CSV,
    deeptmhmm_3line: Path = DEFAULT_DEEPTMHMM_3LINE,
    deeptmhmm_gff3: Path = DEFAULT_DEEPTMHMM_GFF3,
    signalp_summary: Path = DEFAULT_SIGNALP_SUMMARY,
    signalp_json: Path = DEFAULT_SIGNALP_JSON,
    output: Path = DEFAULT_OUTPUT,
    force: bool = False,
) -> None:
    """Merge DeepLoc, DeepTMHMM, and SignalP results for one FASTA batch."""
    input_paths = [
        reference_fasta,
        deeploc_csv,
        deeptmhmm_3line,
        deeptmhmm_gff3,
        signalp_summary,
        signalp_json,
    ]
    missing_paths = [str(path) for path in input_paths if not path.exists()]
    if missing_paths:
        logger.error(f"Missing input files: {missing_paths}")
        raise typer.Exit(code=1)
    if output.exists() and not force:
        logger.error(f"Output already exists; use --force to replace it: {output}")
        raise typer.Exit(code=1)

    reference_records = parse_reference_fasta(reference_fasta)
    reference_ids = {str(record["protein_id"]) for record in reference_records}
    deeploc = parse_deeploc(deeploc_csv)
    deeptmhmm_3line_data = parse_deeptmhmm_3line(deeptmhmm_3line)
    deeptmhmm_gff3_data = parse_deeptmhmm_gff3(deeptmhmm_gff3)
    signalp = parse_signalp_summary(signalp_summary)
    signalp_json_ids, signalp_json_size, signalp_failed_jobs = parse_signalp_json(signalp_json)

    validate_ids("DeepLoc", list(deeploc), reference_ids)
    validate_ids("DeepTMHMM 3line", list(deeptmhmm_3line_data), reference_ids)
    validate_ids("DeepTMHMM GFF3", list(deeptmhmm_gff3_data), reference_ids)
    validate_ids("SignalP summary", list(signalp), reference_ids)
    validate_ids("SignalP JSON", signalp_json_ids, reference_ids)
    if signalp_json_size != len(signalp) or signalp_failed_jobs != 0:
        raise ValueError("SignalP JSON metadata does not match the summary file")

    rows = build_rows(
        reference_records,
        deeploc,
        deeptmhmm_3line_data,
        deeptmhmm_gff3_data,
        signalp,
    )
    write_tsv(output, rows)

    membrane_count = sum(bool(row["is_membrane_associated"]) for row in rows)
    secretory_count = sum(bool(row["is_secretory_pathway"]) for row in rows)
    multi_count = sum(bool(row["is_multilocalized"]) for row in rows)
    logger.info(f"Merged proteins: {len(rows)}")
    logger.info(f"Membrane-associated proteins: {membrane_count}")
    logger.info(f"Secretory-pathway proteins: {secretory_count}")
    logger.info(f"DeepLoc multi-localized proteins: {multi_count}")
    logger.success(f"Wrote merged table: {output}")


if __name__ == "__main__":
    app()
