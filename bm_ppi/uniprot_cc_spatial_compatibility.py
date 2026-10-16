"""Generate Stage 1 Bombyx mori spatially compatible pairs from UniProt CC keywords."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from math import comb
from pathlib import Path

from Bio import SeqIO
from loguru import logger
from tqdm import tqdm
import typer

from bm_ppi.config import INTERIM_DATA_DIR, RAW_DATA_DIR

app = typer.Typer(add_completion=False)

DEFAULT_INPUT = INTERIM_DATA_DIR / "uniprot" / "cc_mapping" / "uniprot_cc_mapping.tsv"
DEFAULT_MAPPING_SUMMARY = (
    INTERIM_DATA_DIR / "uniprot" / "cc_mapping" / "uniprot_cc_summary.json"
)
DEFAULT_FASTA = RAW_DATA_DIR / "GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa"
DEFAULT_OUTPUT_DIR = INTERIM_DATA_DIR / "uniprot" / "cc_spatial_compatibility" / "stage1"
TOTAL_PROTEOME_PAIRS = 90_565_611
EXPECTED_TOTAL_PROTEINS = 13_459
EXPECTED_KNOWN_CC_PROTEINS = 4_481
EXPECTED_CC_KEYWORDS = 54

KNOWN_COLUMNS = [
    "gene_id",
    "ncbi_accession",
    "accession_version",
    "uniprot_accession",
    "mapping_status",
    "cc_keyword_ids",
    "cc_keyword_names",
    "is_ambiguous",
    "source_evidence",
]
EXCLUDED_COLUMNS = [
    "gene_id",
    "ncbi_accession",
    "accession_version",
    "uniprot_accession",
    "mapping_status",
    "exclude_reason",
]
PAIR_COLUMNS = [
    "protein_A",
    "protein_B",
    "ncbi_accession_A",
    "ncbi_accession_B",
    "shared_cc_ids",
    "shared_cc_names",
    "pair_class",
    "mapping_status_A",
    "mapping_status_B",
]
REQUIRED_INPUT_COLUMNS = {
    "gene_id",
    "ncbi_accession",
    "accession_version",
    "uniprot_accession",
    "mapping_status",
    "cc_keyword_id",
    "cc_keyword_name",
    "has_cc_label",
    "is_ambiguous",
    "unknown_reason",
}


def sha256sum(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_cc(row: dict[str, str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Parse and validate paired CC keyword IDs and names."""
    ids = tuple(value.strip() for value in row["cc_keyword_id"].split(";") if value.strip())
    names = tuple(value.strip() for value in row["cc_keyword_name"].split(";") if value.strip())
    if len(ids) != len(names):
        raise ValueError(f"CC ID/name count mismatch for gene {row['gene_id']}")
    if len(ids) != len(set(ids)) or len(names) != len(set(names)):
        raise ValueError(f"Duplicate CC ID or name for gene {row['gene_id']}")
    return ids, names


def exclude_reason(row: dict[str, str]) -> str:
    """Map upstream localization status to a Stage 1 exclusion reason."""
    reason = row["unknown_reason"]
    if reason == "not_mapped_to_bombyx_uniprotkb":
        return "unmapped_to_uniprot"
    if reason == "mapped_but_no_cc_keyword":
        return "mapped_without_cc"
    if reason == "inconsistent_cc_among_multiple_uniprot_entries":
        return "ambiguous_cc_mapping"
    if reason == "non_bombyx_mapping":
        return "non_bombyx_mapping"
    raise ValueError(f"Unsupported exclusion reason for gene {row['gene_id']}: {reason!r}")


def read_input(path: Path) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    """Validate the mapping table and partition known and excluded proteins."""
    known: list[dict[str, object]] = []
    excluded: list[dict[str, str]] = []
    seen: set[str] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(REQUIRED_INPUT_COLUMNS - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Input table is missing columns: {missing}")
        for row in reader:
            gene_id = row["gene_id"]
            if not gene_id or gene_id in seen:
                raise ValueError(f"Empty or duplicate gene ID: {gene_id!r}")
            seen.add(gene_id)
            has_cc = row["has_cc_label"].lower()
            if has_cc not in {"true", "false"}:
                raise ValueError(f"Invalid has_cc_label for gene {gene_id}: {has_cc!r}")
            if has_cc == "true":
                ids, names = split_cc(row)
                if not ids or row["unknown_reason"]:
                    raise ValueError(f"Invalid known-CC record for gene {gene_id}")
                known.append({**row, "cc_ids": ids, "cc_names": names})
            else:
                if row["cc_keyword_id"] or row["cc_keyword_name"]:
                    raise ValueError(f"Excluded gene {gene_id} unexpectedly has CC keywords")
                excluded.append({**row, "stage1_exclude_reason": exclude_reason(row)})
    if len(seen) != EXPECTED_TOTAL_PROTEINS:
        raise ValueError(f"Expected {EXPECTED_TOTAL_PROTEINS} records, found {len(seen)}")
    if len(known) != EXPECTED_KNOWN_CC_PROTEINS:
        raise ValueError(f"Expected {EXPECTED_KNOWN_CC_PROTEINS} known-CC proteins, found {len(known)}")
    return known, excluded


def validate_fasta_links(records: list[dict[str, object]], fasta_path: Path) -> int:
    """Return the number of mapping records absent from the source FASTA."""
    fasta_ids = {record.id for record in SeqIO.parse(fasta_path, "fasta")}
    return sum(str(record["gene_id"]) not in fasta_ids for record in records)


def build_pairs(
    known: list[dict[str, object]],
) -> tuple[dict[int, int], list[str], dict[str, str], Counter[str], int]:
    """Build unique unordered pairs and aggregate shared CC keywords as bit masks."""
    keyword_names: dict[str, str] = {}
    groups: dict[str, list[int]] = defaultdict(list)
    for protein_index, protein in enumerate(known):
        for keyword_id, keyword_name in zip(protein["cc_ids"], protein["cc_names"], strict=True):
            previous = keyword_names.setdefault(keyword_id, keyword_name)
            if previous != keyword_name:
                raise ValueError(f"Conflicting names for {keyword_id}: {previous!r}, {keyword_name!r}")
            groups[keyword_id].append(protein_index)
    keyword_ids = sorted(groups)
    if len(keyword_ids) != EXPECTED_CC_KEYWORDS:
        raise ValueError(f"Expected {EXPECTED_CC_KEYWORDS} CC keywords, found {len(keyword_ids)}")
    pair_masks: dict[int, int] = {}
    group_sizes = Counter({keyword_id: len(groups[keyword_id]) for keyword_id in keyword_ids})
    raw_pair_count = 0
    protein_count = len(known)
    for bit, keyword_id in enumerate(tqdm(keyword_ids, desc="Generating CC-group pairs")):
        members = groups[keyword_id]
        raw_pair_count += comb(len(members), 2)
        mask = 1 << bit
        for left_position, left in enumerate(members[:-1]):
            base = left * protein_count
            for right in members[left_position + 1 :]:
                key = base + right
                pair_masks[key] = pair_masks.get(key, 0) | mask
    return pair_masks, keyword_ids, keyword_names, group_sizes, raw_pair_count


def write_known(path: Path, known: list[dict[str, object]]) -> None:
    """Write the Stage 1 known-localization protein table."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=KNOWN_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in known:
            writer.writerow(
                {
                    "gene_id": row["gene_id"],
                    "ncbi_accession": row["ncbi_accession"],
                    "accession_version": row["accession_version"],
                    "uniprot_accession": row["uniprot_accession"],
                    "mapping_status": row["mapping_status"],
                    "cc_keyword_ids": ";".join(row["cc_ids"]),
                    "cc_keyword_names": ";".join(row["cc_names"]),
                    "is_ambiguous": row["is_ambiguous"],
                    "source_evidence": "UniProt keyword category=Cellular component",
                }
            )


def write_excluded(path: Path, excluded: list[dict[str, str]]) -> None:
    """Write proteins excluded from the known-localization Stage 1 branch."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=EXCLUDED_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in excluded:
            writer.writerow(
                {
                    "gene_id": row["gene_id"],
                    "ncbi_accession": row["ncbi_accession"],
                    "accession_version": row["accession_version"],
                    "uniprot_accession": row["uniprot_accession"],
                    "mapping_status": row["mapping_status"],
                    "exclude_reason": row["stage1_exclude_reason"],
                }
            )


def write_pairs(
    path: Path,
    known: list[dict[str, object]],
    pair_masks: dict[int, int],
    keyword_ids: list[str],
    keyword_names: dict[str, str],
) -> None:
    """Write unique spatially compatible protein pairs."""
    protein_count = len(known)
    with gzip.open(path, "wt", newline="", encoding="utf-8", compresslevel=6) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(PAIR_COLUMNS)
        for key in tqdm(sorted(pair_masks), desc="Writing compatible pairs"):
            left, right = divmod(key, protein_count)
            mask = pair_masks[key]
            shared_ids = [keyword_id for bit, keyword_id in enumerate(keyword_ids) if mask & (1 << bit)]
            protein_a = known[left]
            protein_b = known[right]
            writer.writerow(
                [
                    protein_a["gene_id"],
                    protein_b["gene_id"],
                    protein_a["ncbi_accession"],
                    protein_b["ncbi_accession"],
                    ";".join(shared_ids),
                    ";".join(keyword_names[keyword_id] for keyword_id in shared_ids),
                    "Same locality",
                    protein_a["mapping_status"],
                    protein_b["mapping_status"],
                ]
            )


def render_report(summary: dict[str, object]) -> str:
    """Render a concise Markdown summary."""
    counts = summary["counts"]
    rates = summary["rates"]
    return f"""# Stage 1：UniProt CC-known spatial compatibility result

## 方法

仅纳入具有经核验 UniProt Cellular component keyword 的家蚕蛋白。使用 CC keyword ID 比较定位集合；至少共享一个 keyword 的无序蛋白对判定为 `Same locality`。本阶段未使用 DeepLoc、GO CC、跨物种同源定位或人工相容规则。

## 统计

| 指标 | 数值 |
|---|---:|
| NCBI 蛋白总数 | {counts['total_proteins']:,} |
| CC 已知定位蛋白 | {counts['known_cc_proteins']:,} |
| 本阶段排除或定位未知蛋白 | {counts['excluded_or_unknown_proteins']:,} |
| CC keyword 种类 | {counts['cc_keyword_count']:,} |
| 已知定位蛋白全部无序组合 | {counts['all_known_cc_pairs']:,} |
| 去重前标签内蛋白对 | {counts['raw_within_keyword_pairs']:,} |
| 去重后 Same locality 蛋白对 | {counts['same_locality_pairs']:,} |
| 推导的 Different locality 蛋白对 | {counts['different_locality_pairs']:,} |
| 多标签重复折叠数 | {counts['duplicate_pair_occurrences_collapsed']:,} |
| Stage 1 内部蛋白对减少率 | {rates['stage1_reduction_fraction']:.2%} |
| Same locality 占全蛋白组组合比例 | {rates['same_pairs_vs_full_proteome_fraction']:.2%} |
| FASTA 关联缺失记录 | {counts['missing_fasta_links']:,} |

## 输出解释

`stage1_spatially_compatible_pairs.tsv.gz` 只包含 `Same locality` 候选。`Different locality` 数量由全部已知定位组合减去 Same locality 数量获得，不写出明细。定位未知蛋白保留在独立清单中，本结果不是完整家蚕 PPI 网络。
"""


@app.command()
def main(
    input_path: Path = DEFAULT_INPUT,
    mapping_summary_path: Path = DEFAULT_MAPPING_SUMMARY,
    fasta_path: Path = DEFAULT_FASTA,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> None:
    """Generate the Stage 1 UniProt CC-known spatial compatibility result."""
    for path in (input_path, mapping_summary_path, fasta_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    known, excluded = read_input(input_path)
    all_records = known + excluded
    missing_fasta_links = validate_fasta_links(all_records, fasta_path)
    if missing_fasta_links:
        raise ValueError(f"{missing_fasta_links} mapping records cannot be linked to the FASTA")
    pair_masks, keyword_ids, keyword_names, group_sizes, raw_pair_count = build_pairs(known)
    output_dir.mkdir(parents=True, exist_ok=True)
    known_path = output_dir / "stage1_cc_known_proteins.tsv"
    pair_path = output_dir / "stage1_spatially_compatible_pairs.tsv.gz"
    excluded_path = output_dir / "stage1_excluded_or_unknown_proteins.tsv"
    summary_path = output_dir / "stage1_cc_spatial_summary.json"
    report_path = output_dir / "stage1_cc_spatial_summary.md"
    write_known(known_path, known)
    write_excluded(excluded_path, excluded)
    write_pairs(pair_path, known, pair_masks, keyword_ids, keyword_names)
    mapping_summary = json.loads(mapping_summary_path.read_text(encoding="utf-8"))
    all_known_pairs = comb(len(known), 2)
    same_pairs = len(pair_masks)
    different_pairs = all_known_pairs - same_pairs
    exclusion_counts = Counter(row["stage1_exclude_reason"] for row in excluded)
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "Stage 1: UniProt CC-known spatial compatibility result",
        "method": {
            "evidence": "UniProt keywords with category == Cellular component",
            "comparison_key": "UniProt CC keyword ID",
            "pair_rule": "Same locality if the two CC keyword ID sets intersect",
            "pair_scope": "unordered non-self pairs among known-CC proteins",
            "excluded_evidence": ["DeepLoc", "GO Cellular Component", "cross-species transfer"],
        },
        "uniprot_metadata": mapping_summary.get("query_metadata", {}),
        "inputs": {
            "mapping_table": {"path": str(input_path.resolve()), "sha256": sha256sum(input_path)},
            "mapping_summary": {
                "path": str(mapping_summary_path.resolve()),
                "sha256": sha256sum(mapping_summary_path),
            },
            "source_fasta": {"path": str(fasta_path.resolve()), "sha256": sha256sum(fasta_path)},
        },
        "counts": {
            "total_proteins": len(all_records),
            "known_cc_proteins": len(known),
            "excluded_or_unknown_proteins": len(excluded),
            "cc_keyword_count": len(keyword_ids),
            "all_known_cc_pairs": all_known_pairs,
            "raw_within_keyword_pairs": raw_pair_count,
            "same_locality_pairs": same_pairs,
            "different_locality_pairs": different_pairs,
            "unknown_locality_pairs_in_stage1": 0,
            "duplicate_pair_occurrences_collapsed": raw_pair_count - same_pairs,
            "missing_fasta_links": missing_fasta_links,
            "exclusion_reason_counts": dict(sorted(exclusion_counts.items())),
        },
        "rates": {
            "stage1_retained_fraction": same_pairs / all_known_pairs,
            "stage1_reduction_fraction": different_pairs / all_known_pairs,
            "same_pairs_vs_full_proteome_fraction": same_pairs / TOTAL_PROTEOME_PAIRS,
        },
        "cc_keyword_group_sizes": {
            keyword_id: {"name": keyword_names[keyword_id], "protein_count": group_sizes[keyword_id]}
            for keyword_id in keyword_ids
        },
        "outputs": {
            "known_proteins": str(known_path.resolve()),
            "compatible_pairs": str(pair_path.resolve()),
            "excluded_or_unknown": str(excluded_path.resolve()),
            "summary_json": str(summary_path.resolve()),
            "summary_markdown": str(report_path.resolve()),
        },
        "validation": {
            "known_plus_excluded_equals_total": len(known) + len(excluded) == len(all_records),
            "same_plus_different_equals_all_known_pairs": same_pairs + different_pairs
            == all_known_pairs,
            "candidate_pairs_do_not_exceed_pair_upper_bound": same_pairs <= all_known_pairs,
        },
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    report_path.write_text(render_report(summary), encoding="utf-8")
    logger.success(
        "Stage 1 complete: proteins={}, Same locality pairs={}, Different locality pairs={}",
        len(known),
        same_pairs,
        different_pairs,
    )


if __name__ == "__main__":
    app()
