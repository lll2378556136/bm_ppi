"""Quality-control script for the Bombyx mori non-redundant protein set.

This module performs QC on:

1. ``data/raw/protein_nr95.faa`` — the CD-HIT non-redundant representative
   protein FASTA.
2. ``data/raw/protein_nr95.faa.clstr`` — the CD-HIT cluster report.
3. Optionally ``data/raw/protein.faa`` — the original full proteome, if present.

Outputs:

- ``data/interim/qc_summary.json`` — aggregated summary statistics.
- ``data/interim/qc_report.tsv`` — per-sequence/per-cluster QC table.

Run with:

    pixi run python -m bm_ppi.qc

"""

from __future__ import annotations

import csv
import json
import re
import statistics
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import typer
from loguru import logger
from tqdm import tqdm

from bm_ppi.config import INTERIM_DATA_DIR, RAW_DATA_DIR

app = typer.Typer()

# Standard/IUPAC amino-acid letters plus the stop codon symbol.
# Ambiguous characters (B, Z, X) and selenocysteine/pyrrolysine (U, O) are
# treated as acceptable for this QC.
ALLOWED_AA = set("ACEDGFIHKMLNQPSRTVWYOUBZX*")

# Regex for one CD-HIT member line, e.g.
#   0\t1122aa, >NP_001040203.1... at 95.10%
#   1\t18125aa, >XP_021202734.2... *
CLSTR_MEMBER_RE = re.compile(
    r"^\d+\s+"                       # index number and whitespace
    r"(?P<length>\d+)aa,\s+"         # length followed by 'aa, '
    r">(?P<seqid>\S+)\.\.\.\s+"     # '>' then sequence id then '... '
    r"(?P<flag>.*)$"                 # remainder (percentage or '*')
)


def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    """Parse a FASTA file into a list of (header, sequence) tuples.

    The header is the full text after the leading '>' character.  The
    sequence has all inter-line whitespace removed and is upper-cased.
    Empty sequences are preserved so that downstream QC can flag them.

    Args:
        path: Path to the FASTA file.

    Returns:
        A list of (header, sequence) tuples in file order.
    """
    records: List[Tuple[str, str]] = []
    current_header: Optional[str] = None
    current_seq_parts: List[str] = []

    logger.info(f"Parsing FASTA file: {path}")
    with path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                # Blank lines are ignored; they do not terminate a sequence.
                continue
            if line.startswith(">"):
                # Store the previous record before starting a new one.
                if current_header is not None:
                    records.append((current_header, "".join(current_seq_parts).upper()))
                current_header = line[1:].strip()
                current_seq_parts = []
            else:
                if current_header is None:
                    # Sequence data appears before any header — malformed file.
                    logger.warning("FASTA sequence data found before the first header.")
                    continue
                current_seq_parts.append(line)

    # Flush the final record.
    if current_header is not None:
        records.append((current_header, "".join(current_seq_parts).upper()))

    return records


def compute_n50(lengths: List[int]) -> int:
    """Compute the N50 length statistic for a collection of sequence lengths.

    N50 is the length of the shortest sequence for which the cumulative
    length of all sequences of that length or longer is at least half of
    the total length.

    Args:
        lengths: List of sequence lengths.

    Returns:
        The N50 value.  Returns 0 if ``lengths`` is empty.
    """
    if not lengths:
        return 0
    total = sum(lengths)
    half = total / 2.0
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= half:
            return length
    return 0


def qc_fasta(records: List[Tuple[str, str]]) -> Dict:
    """Run QC checks on a parsed FASTA file.

    Args:
        records: List of (header, sequence) tuples from ``parse_fasta``.

    Returns:
        A dictionary containing summary statistics and per-record flags.
    """
    total = len(records)
    headers = [header for header, _ in records]
    seq_ids = [header.split()[0] if header.split() else "" for header, _ in records]

    # Duplicate detection
    header_counts = Counter(headers)
    duplicate_headers = sum(1 for count in header_counts.values() if count > 1)
    seq_id_counts = Counter(seq_ids)
    duplicate_ids = sum(1 for count in seq_id_counts.values() if count > 1)

    lengths = []
    empty_or_malformed = 0
    non_standard_count = 0
    stop_codon_count = 0
    short_count = 0

    per_record_flags = []

    for header, seq in records:
        seq_id = header.split()[0] if header.split() else ""
        length = len(seq)
        lengths.append(length)

        # Empty or malformed: header existed but produced no sequence data.
        is_empty = length == 0
        if is_empty:
            empty_or_malformed += 1

        # Non-standard amino acids (anything outside the allowed alphabet).
        has_non_standard = any(char not in ALLOWED_AA for char in seq)
        if has_non_standard:
            non_standard_count += 1

        # Internal stop codon: '*' anywhere in the sequence.
        has_stop = "*" in seq
        if has_stop:
            stop_codon_count += 1

        # Sequences shorter than 50 amino acids.
        is_short = 0 < length < 50
        if is_short:
            short_count += 1

        per_record_flags.append(
            {
                "seq_id": seq_id,
                "header": header,
                "length": length,
                "has_non_standard_aa": has_non_standard,
                "has_internal_stop": has_stop,
                "is_short": is_short,
                "is_empty_or_malformed": is_empty,
            }
        )

    # Length statistics
    if lengths:
        length_stats = {
            "min": min(lengths),
            "max": max(lengths),
            "mean": round(statistics.mean(lengths), 2),
            "median": round(statistics.median(lengths), 2),
            "n50": compute_n50(lengths),
            "total_aa": sum(lengths),
        }
    else:
        length_stats = {"min": 0, "max": 0, "mean": 0.0, "median": 0.0, "n50": 0, "total_aa": 0}

    summary = {
        "total_sequences": total,
        "duplicate_headers": duplicate_headers,
        "duplicate_ids": duplicate_ids,
        "length_stats": length_stats,
        "non_standard_aa": {
            "count": non_standard_count,
            "fraction": round(non_standard_count / total, 6) if total else 0.0,
        },
        "stop_codon_count": stop_codon_count,
        "short_sequences_count": short_count,
        "empty_or_malformed_count": empty_or_malformed,
    }

    return {"summary": summary, "per_record": per_record_flags}


def parse_cdhit_clusters(path: Path) -> Dict:
    """Parse a CD-HIT ``.clstr`` file.

    Args:
        path: Path to the ``.clstr`` file.

    Returns:
        A dictionary with cluster metadata and a per-cluster list of members.
        Cluster ids are the numeric ids written by CD-HIT.
    """
    clusters: List[Dict] = []
    current_cluster: Optional[Dict] = None

    logger.info(f"Parsing CD-HIT cluster file: {path}")
    with path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                # Start a new cluster; store the previous one if it exists.
                if current_cluster is not None:
                    clusters.append(current_cluster)
                cluster_id = line.replace(">Cluster", "").strip()
                current_cluster = {
                    "cluster_id": cluster_id,
                    "members": [],
                    "representative_id": None,
                }
            else:
                if current_cluster is None:
                    logger.warning("CD-HIT member line found before first cluster header.")
                    continue
                match = CLSTR_MEMBER_RE.match(line)
                if not match:
                    logger.warning(f"Could not parse CD-HIT member line: {line}")
                    continue
                seq_id = match.group("seqid")
                reported_length = int(match.group("length"))
                flag = match.group("flag").strip()
                is_representative = flag.endswith("*")
                member = {
                    "seq_id": seq_id,
                    "reported_length": reported_length,
                    "is_representative": is_representative,
                }
                current_cluster["members"].append(member)
                if is_representative:
                    current_cluster["representative_id"] = seq_id

    if current_cluster is not None:
        clusters.append(current_cluster)

    return {"clusters": clusters}


def qc_cdhit(clusters_data: Dict, seq_lengths: Dict[str, int]) -> Dict:
    """Run QC checks on parsed CD-HIT clusters.

    Args:
        clusters_data: Output from ``parse_cdhit_clusters``.
        seq_lengths: Mapping from sequence id to actual FASTA sequence length.

    Returns:
        A dictionary containing summary statistics and mismatch reports.
    """
    clusters = clusters_data["clusters"]
    total_clusters = len(clusters)
    member_counts = [len(c["members"]) for c in clusters]
    singletons = sum(1 for count in member_counts if count == 1)
    multi_member = sum(1 for count in member_counts if count > 1)
    size_distribution = dict(Counter(member_counts))

    reps_in_clstr = set()
    members_in_clstr = set()
    cluster_lookup: Dict[str, Dict] = {}
    length_mismatches: List[Dict] = []
    missing_representative: List[str] = []

    for cluster in clusters:
        rep_id = cluster["representative_id"]
        if rep_id is None:
            logger.warning(f"Cluster {cluster['cluster_id']} has no representative marked with *.")
            continue
        reps_in_clstr.add(rep_id)
        cluster_lookup[rep_id] = cluster
        for member in cluster["members"]:
            members_in_clstr.add(member["seq_id"])
            if member["is_representative"]:
                actual_length = seq_lengths.get(rep_id)
                if actual_length is not None and actual_length != member["reported_length"]:
                    length_mismatches.append(
                        {
                            "seq_id": rep_id,
                            "cluster_id": cluster["cluster_id"],
                            "clstr_reported_length": member["reported_length"],
                            "actual_fasta_length": actual_length,
                        }
                    )

    headers_in_faa = set(seq_lengths.keys())

    missing_representatives_in_faa = sorted(reps_in_clstr - headers_in_faa)
    missing_headers_in_clstr = sorted(headers_in_faa - reps_in_clstr)

    # Report clusters without a representative separately.
    for cluster in clusters:
        if cluster["representative_id"] is None:
            missing_representative.append(cluster["cluster_id"])

    summary = {
        "total_clusters": total_clusters,
        "singletons": singletons,
        "multi_member_clusters": multi_member,
        "cluster_size_distribution": size_distribution,
        "clusters_without_representative": missing_representative,
        "representative_mismatches": {
            "reps_in_clstr_not_in_faa": missing_representatives_in_faa,
            "headers_in_faa_not_in_clstr": missing_headers_in_clstr,
            "count_reps_in_clstr_not_in_faa": len(missing_representatives_in_faa),
            "count_headers_in_faa_not_in_clstr": len(missing_headers_in_clstr),
        },
        "length_mismatches": {
            "count": len(length_mismatches),
            "examples": length_mismatches[:20],  # limit to first 20 to keep JSON readable
        },
    }

    return {"summary": summary, "cluster_lookup": cluster_lookup}


def optional_full_proteome_check(full_faa_path: Path, representative_ids: List[str]) -> Dict:
    """Optionally check whether representative IDs are present in the full proteome.

    Args:
        full_faa_path: Path to the original full proteome FASTA.
        representative_ids: List of representative sequence ids.

    Returns:
        A dictionary describing the result of the subset check.
    """
    if not full_faa_path.exists():
        return {"checked": False, "reason": "file not found"}

    logger.info(f"Optional cross-check against full proteome: {full_faa_path}")
    try:
        full_records = parse_fasta(full_faa_path)
    except Exception as exc:  # noqa: BLE001
        logger.error(f"Could not read full proteome file: {exc}")
        return {"checked": False, "reason": f"read error: {exc}"}

    full_ids = {header.split()[0] if header.split() else "" for header, _ in full_records}
    rep_set = set(representative_ids)
    in_full = sorted(rep_set & full_ids)
    not_in_full = sorted(rep_set - full_ids)

    return {
        "checked": True,
        "total_sequences": len(full_records),
        "unique_ids": len(full_ids),
        "representatives_in_full": len(in_full),
        "representatives_not_in_full": not_in_full,
        "count_representatives_not_in_full": len(not_in_full),
    }


def write_json_summary(path: Path, summary: Dict) -> None:
    """Write the QC summary to a JSON file with stable key ordering."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False, sort_keys=True)
    logger.info(f"Wrote JSON summary to {path}")


def write_tsv_report(path: Path, per_record_rows: List[Dict]) -> None:
    """Write the per-record QC table to a TSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not per_record_rows:
        path.touch()
        logger.info(f"Created empty TSV report at {path}")
        return

    fieldnames = list(per_record_rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(per_record_rows)
    logger.info(f"Wrote per-record TSV report to {path}")


@app.command()
def main(
    nr_faa_path: Path = RAW_DATA_DIR / "protein_nr95.faa",
    clstr_path: Path = RAW_DATA_DIR / "protein_nr95.faa.clstr",
    full_faa_path: Path = RAW_DATA_DIR / "protein.faa",
    json_output_path: Path = INTERIM_DATA_DIR / "qc_summary.json",
    tsv_output_path: Path = INTERIM_DATA_DIR / "qc_report.tsv",
) -> None:
    """Run QC on the Bombyx mori non-redundant protein set and CD-HIT mapping."""
    logger.info("Starting QC for Bombyx mori non-redundant protein set.")

    # Ensure the output directory exists; this is a non-destructive operation.
    INTERIM_DATA_DIR.mkdir(parents=True, exist_ok=True)

    # ---- 1. FASTA QC ----
    if not nr_faa_path.exists():
        logger.error(f"Representative FASTA not found: {nr_faa_path}")
        raise typer.Exit(code=1)

    records = parse_fasta(nr_faa_path)
    fasta_qc_result = qc_fasta(records)
    fasta_summary = fasta_qc_result["summary"]
    per_record_flags = fasta_qc_result["per_record"]

    logger.info(f"FASTA total sequences: {fasta_summary['total_sequences']}")
    logger.info(
        f"FASTA duplicate headers / ids: "
        f"{fasta_summary['duplicate_headers']} / {fasta_summary['duplicate_ids']}"
    )
    logger.info(
        f"FASTA length stats (min/max/mean/median/N50): "
        f"{fasta_summary['length_stats']['min']}/"
        f"{fasta_summary['length_stats']['max']}/"
        f"{fasta_summary['length_stats']['mean']}/"
        f"{fasta_summary['length_stats']['median']}/"
        f"{fasta_summary['length_stats']['n50']}"
    )
    logger.info(
        f"FASTA non-standard AA: {fasta_summary['non_standard_aa']['count']} "
        f"({fasta_summary['non_standard_aa']['fraction']})"
    )
    logger.info(f"FASTA internal stop codons: {fasta_summary['stop_codon_count']}")
    logger.info(f"FASTA short (<50 aa) sequences: {fasta_summary['short_sequences_count']}")
    logger.info(
        f"FASTA empty or malformed sequences: {fasta_summary['empty_or_malformed_count']}"
    )

    # ---- 2. CD-HIT cluster QC ----
    if not clstr_path.exists():
        logger.error(f"CD-HIT cluster file not found: {clstr_path}")
        raise typer.Exit(code=1)

    clusters_data = parse_cdhit_clusters(clstr_path)
    seq_lengths = {row["seq_id"]: row["length"] for row in per_record_flags}
    cdhit_qc_result = qc_cdhit(clusters_data, seq_lengths)
    cdhit_summary = cdhit_qc_result["summary"]
    cluster_lookup = cdhit_qc_result["cluster_lookup"]

    logger.info(f"CD-HIT total clusters: {cdhit_summary['total_clusters']}")
    logger.info(f"CD-HIT singletons: {cdhit_summary['singletons']}")
    logger.info(f"CD-HIT multi-member clusters: {cdhit_summary['multi_member_clusters']}")
    logger.info(
        f"CD-HIT representatives missing from FASTA: "
        f"{cdhit_summary['representative_mismatches']['count_reps_in_clstr_not_in_faa']}"
    )
    logger.info(
        f"CD-HIT FASTA headers missing as representatives: "
        f"{cdhit_summary['representative_mismatches']['count_headers_in_faa_not_in_clstr']}"
    )
    logger.info(f"CD-HIT representative length mismatches: {cdhit_summary['length_mismatches']['count']}")

    # ---- 3. Optional full-proteome cross-check ----
    full_proteome_result = optional_full_proteome_check(
        full_faa_path, list(seq_lengths.keys())
    )
    if full_proteome_result["checked"]:
        logger.info(
            f"Full proteome total sequences: {full_proteome_result['total_sequences']}"
        )
        logger.info(
            f"Representatives found in full proteome: "
            f"{full_proteome_result['representatives_in_full']} / {len(seq_lengths)}"
        )
    else:
        logger.info(
            f"Skipping full proteome cross-check: {full_proteome_result.get('reason')}"
        )

    # ---- 4. Build per-record TSV rows ----
    per_record_rows: List[Dict] = []
    for flags in tqdm(per_record_flags, desc="Building TSV rows"):
        seq_id = flags["seq_id"]
        cluster = cluster_lookup.get(seq_id)
        if cluster is not None:
            cluster_id = cluster["cluster_id"]
            is_representative = True
            member_count = len(cluster["members"])
            rep_member = next(
                (m for m in cluster["members"] if m["seq_id"] == seq_id), None
            )
            clstr_length = rep_member["reported_length"] if rep_member else None
        else:
            cluster_id = None
            is_representative = False
            member_count = None
            clstr_length = None

        actual_length = flags["length"]
        length_matches = (
            None
            if clstr_length is None
            else (actual_length == clstr_length)
        )

        per_record_rows.append(
            {
                "seq_id": seq_id,
                "header": flags["header"],
                "length": actual_length,
                "cluster_id": cluster_id,
                "is_representative": is_representative,
                "cluster_member_count": member_count,
                "has_non_standard_aa": flags["has_non_standard_aa"],
                "has_internal_stop": flags["has_internal_stop"],
                "is_short": flags["is_short"],
                "is_empty_or_malformed": flags["is_empty_or_malformed"],
                "clstr_reported_length": clstr_length,
                "length_matches_clstr": length_matches,
            }
        )

    # ---- 5. Write outputs ----
    final_summary = {
        "protein_nr95_faa": fasta_summary,
        "cdhit_clusters": cdhit_summary,
        "full_proteome_check": full_proteome_result,
        "input_files": {
            "nr_faa": str(nr_faa_path.resolve()),
            "clstr": str(clstr_path.resolve()),
            "full_faa": str(full_faa_path.resolve()) if full_faa_path.exists() else None,
        },
        "output_files": {
            "json_summary": str(json_output_path.resolve()),
            "tsv_report": str(tsv_output_path.resolve()),
        },
    }

    write_json_summary(json_output_path, final_summary)
    write_tsv_report(tsv_output_path, per_record_rows)

    logger.success("QC complete.")


if __name__ == "__main__":
    app()
