"""Map the working set of Bombyx mori NCBI proteins to UniProtKB and extract
Cellular Component (CC) localization keywords.

FASTA IDs in ``GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa`` are
NCBI Gene IDs, not protein accessions. NCBI protein accessions are therefore
re-attached in two steps: (1) parse the assembly GFF3 for the official
gene -> protein accession links, (2) accept the accessions of a record's own
gene whose sequences match the record exactly; a global sequence match is
used only as an explicitly flagged fallback.

All UniProt queries use the official REST API of uniprot.org (endpoints, query
date and release are recorded in the summary JSON):

1. ID Mapping: ``POST /idmapping/run`` with ``from=RefSeq_Protein`` and
   ``to=UniProtKB`` on version-stripped accessions. Only mappings whose
   target entry belongs to Bombyx mori (taxonomy ID 7091) are kept; mappings
   into other species are counted and discarded (no cross-species inheritance
   of localization labels).
2. Entry details: ``GET /uniprotkb/search`` in JSON format to collect
   keyword IDs, names and categories. Keywords with category
   ``Cellular component`` are treated as CC localization tags. Reviewed and
   unreviewed entries are both accepted. GO cellular component terms and the
   ``SUBCELLULAR LOCATION`` free-text comment are deliberately NOT used.

Classification rules (a consistent partition of all records):

- has_cc / localization known:
    * unique UniProt entry with >= 1 CC keyword, or
    * multiple UniProt entries whose CC keyword sets are identical and
      non-empty (kept, but flagged as ambiguous).
- localization unknown (everything else):
    * record not mapped to any Bombyx mori UniProtKB entry,
    * mapped but the UniProt entry(ies) carry no CC keyword,
    * multiple mapped entries with inconsistent CC keyword sets
      (flagged as mapping_ambiguous).

The analysis is deterministic (no random sampling), so no random seed is set.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.error
import urllib.parse
import urllib.request

from Bio import SeqIO
from loguru import logger
from tqdm import tqdm
import typer

from bm_ppi.config import EXTERNAL_DATA_DIR, INTERIM_DATA_DIR, RAW_DATA_DIR

app = typer.Typer(add_completion=False)

TAXON_ID = 7091  # Bombyx mori
CC_CATEGORY = "Cellular component"
UNIPROT_REST = "https://rest.uniprot.org"
IDMAPPING_FROM = "RefSeq_Protein"
IDMAPPING_TO = "UniProtKB"
IDMAPPING_FIELDS = "accession,id,reviewed,organism_id,organism_name,protein_name,keyword"
ENTRY_URL_TEMPLATE = "https://www.uniprot.org/uniprotkb/{acc}"
GFF_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/269/925/"
    "GCF_030269925.1_ASM3026992v2/GCF_030269925.1_ASM3026992v2_genomic.gff.gz"
)
MAPPING_METHOD = f"UniProt_ID_Mapping_{IDMAPPING_FROM}_to_{IDMAPPING_TO}"
RETRY_DELAYS = [2, 5, 10, 20, 40, 60]

DETAIL_COLUMNS = [
    "gene_id",
    "description",
    "sequence_length",
    "ncbi_accession",
    "accession_version",
    "candidate_accessions_versioned",
    "n_candidate_accessions",
    "accession_derivation",
    "is_gene_longest_confirmed",
    "uniprot_accession",
    "n_uniprot_entries",
    "mapping_method",
    "mapping_status",
    "uniprot_reviewed",
    "uniprot_protein_name",
    "cc_keyword_id",
    "cc_keyword_name",
    "has_cc_label",
    "is_ambiguous",
    "unknown_reason",
    "uniprot_entry_url",
]

LONG_COLUMNS = [
    "gene_id",
    "ncbi_accession",
    "accession_version",
    "ncbi_accession_versioned",
    "accession_derivation",
    "uniprot_accession",
    "uniprot_reviewed",
    "uniprot_protein_name",
    "cc_keyword_id",
    "cc_keyword_name",
    "entry_has_cc",
    "record_status",
    "record_has_cc",
]

# Populated from the X-UniProt-Release / X-UniProt-Release-Date response headers.
_RELEASE_INFO: dict[str, str] = {}


class NoRedirectHandler(urllib.request.HTTPRedirectHandler):
    """Prevent urllib from following redirects for status polling."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


@dataclass(frozen=True, slots=True)
class ProteinRecord:
    """One FASTA record of the working protein set (one per NCBI gene)."""

    fasta_order: int
    gene_id: str
    description: str
    sequence: str


# --------------------------------------------------------------------------- #
# Generic helpers
# --------------------------------------------------------------------------- #


def http_request(
    url: str,
    data: bytes | None = None,
    timeout: int = 60,
    retries: int = 6,
    follow_redirects: bool = True,
) -> tuple[bytes, dict[str, str]]:
    """GET or POST a urlencoded form; retry transient failures."""
    last_error: Exception | None = None
    opener = urllib.request.build_opener(
        *([] if follow_redirects else [NoRedirectHandler()])
    )
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                data=data,
                headers={
                    "Accept": "*/*",
                    "Connection": "close",
                    "User-Agent": "bm_ppi-uniprot-cc-mapping/1.0",
                },
            )
            with opener.open(req, timeout=timeout) as resp:
                body = resp.read()
                headers = {k: v for k, v in resp.headers.items()}
                for key in ("X-UniProt-Release", "X-UniProt-Release-Date"):
                    if key in headers and key not in _RELEASE_INFO:
                        _RELEASE_INFO[key] = headers[key]
                return body, headers
        except urllib.error.HTTPError as exc:
            error_body = exc.read()
            if not follow_redirects and exc.code in (301, 302, 303, 307, 308):
                headers = {k: v for k, v in exc.headers.items()}
                return error_body, headers
            if exc.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                logger.warning(f"HTTP {exc.code} from {url[:100]} ... retrying")
                time.sleep(RETRY_DELAYS[min(attempt, len(RETRY_DELAYS) - 1)])
                continue
            raise RuntimeError(f"HTTP {exc.code} for {url}: {error_body[:200]!r}") from exc
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last_error = exc
            if attempt == retries - 1:
                break
            logger.warning(f"Network error on {url[:100]}: {exc}; retrying")
            time.sleep(RETRY_DELAYS[min(attempt, len(RETRY_DELAYS) - 1)])
    raise RuntimeError(f"Request failed after {retries} attempts: {url}") from last_error


def sha256_of_file(path: Path) -> str:
    """Return the SHA-256 hex digest of a file (streamed for large files)."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def split_version(accession: str) -> tuple[str, str]:
    """Split a versioned accession into (base, version)."""
    base, _, version = accession.partition(".")
    return base, version or "NA"


def gff_attr(attributes: str, key: str) -> str | None:
    """Extract a single GFF3 attribute value (first occurrence)."""
    match = re.search(rf"(?:^|;){key}=([^;]+)", attributes)
    return match.group(1) if match else None


# --------------------------------------------------------------------------- #
# Step 1: inputs and GFF-derived gene -> protein accession links
# --------------------------------------------------------------------------- #


def load_records(fasta_path: Path) -> list[ProteinRecord]:
    """Load the working FASTA and validate that IDs are unique."""
    records: list[ProteinRecord] = []
    seen_ids: set[str] = set()
    for index, record in enumerate(SeqIO.parse(str(fasta_path), "fasta")):
        if record.id in seen_ids:
            raise ValueError(f"Duplicate FASTA ID in working set: {record.id}")
        seen_ids.add(record.id)
        records.append(
            ProteinRecord(
                fasta_order=index,
                gene_id=str(record.id),
                description=str(record.description),
                sequence=str(record.seq).upper(),
            )
        )
    if not records:
        raise ValueError(f"No records found in {fasta_path}")
    logger.info(f"Loaded {len(records)} protein records from {fasta_path.name}")
    return records


def load_full_proteome(fasta_path: Path) -> dict[str, str]:
    """Load the full RefSeq proteome as {versioned accession: sequence}."""
    accessions: dict[str, str] = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        if record.id in accessions:
            raise ValueError(f"Duplicate accession in full proteome: {record.id}")
        accessions[str(record.id)] = str(record.seq).upper()
    logger.info(f"Loaded {len(accessions)} proteins from {fasta_path.name}")
    return accessions


def download_gff_if_needed(gff_path: Path) -> None:
    """Download the assembly GFF3 unless it already exists locally."""
    if gff_path.exists() and gff_path.stat().st_size > 0:
        logger.info(f"GFF already present: {gff_path} ({gff_path.stat().st_size / 1e6:.1f} MB)")
        return
    gff_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Downloading {GFF_URL}")
    body, _ = http_request(GFF_URL, timeout=600)
    gff_path.write_bytes(body)
    logger.info(f"Saved GFF to {gff_path} ({gff_path.stat().st_size / 1e6:.1f} MB)")


def parse_gff_gene_proteins(gff_path: Path) -> dict[str, set[str]]:
    """Map NCBI Gene ID -> set of versioned protein accessions from GFF3.

    CDS features carry protein_id attributes; their genes are resolved either
    from the GeneID Db_xref on the CDS/mRNA line or via the mRNA -> gene
    Parent chain. Both passes are streamed to keep memory low.
    """
    rna2gene: dict[str, str] = {}
    with gzip.open(gff_path, "rt") as fh:
        for line in fh:
            if "\tmRNA\t" not in line:
                continue
            attributes = line.rstrip("\n").split("\t")[-1]
            rna_id = gff_attr(attributes, "ID")
            if not rna_id:
                continue
            match = re.search(r"GeneID:(\d+)", attributes)
            if match:
                rna2gene[rna_id] = match.group(1)
                continue
            parent = gff_attr(attributes, "Parent") or ""
            parent_match = re.search(r"gene-(\d+)", parent)
            if parent_match:
                rna2gene[rna_id] = parent_match.group(1)
    logger.info(f"GFF pass 1: {len(rna2gene)} mRNA transcripts linked to genes")

    gene2accs: dict[str, set[str]] = defaultdict(set)
    n_cds = 0
    with gzip.open(gff_path, "rt") as fh:
        for line in fh:
            if "\tCDS\t" not in line:
                continue
            attributes = line.rstrip("\n").split("\t")[-1]
            protein_id = gff_attr(attributes, "protein_id")
            if not protein_id:
                continue
            n_cds += 1
            gene_ids: set[str] = set()
            match = re.search(r"GeneID:(\d+)", attributes)
            if match:
                gene_ids.add(match.group(1))
            for parent in (gff_attr(attributes, "Parent") or "").split(","):
                if parent in rna2gene:
                    gene_ids.add(rna2gene[parent])
            if not gene_ids:
                logger.debug(f"CDS without gene link skipped: {protein_id}")
                continue
            for gene_id in gene_ids:
                gene2accs[gene_id].add(protein_id)
    logger.info(
        f"GFF pass 2: {n_cds} CDS features -> {len(gene2accs)} genes with proteins"
    )
    return gene2accs


# --------------------------------------------------------------------------- #
# Step 2: attach NCBI protein accessions to every record
# --------------------------------------------------------------------------- #


def derive_candidates(
    records: list[ProteinRecord],
    acc_seqs: dict[str, str],
    gene2accs: dict[str, set[str]],
) -> list[dict]:
    """Attach candidate NCBI protein accessions to each record.

    Primary route: accessions of the record's own gene with an exact sequence
    match. Fallbacks (explicitly flagged): global sequence match when the gene
    is missing from the GFF or when none of the gene's accessions match.
    """
    seq2accs: dict[str, list[str]] = defaultdict(list)
    for acc, seq in acc_seqs.items():
        seq2accs[seq].append(acc)
    for accs in seq2accs.values():
        accs.sort()

    gff_accs = {acc for accs in gene2accs.values() for acc in accs}
    missing_in_proteome = sorted(acc for acc in gff_accs if acc not in acc_seqs)
    if missing_in_proteome:
        logger.warning(
            f"{len(missing_in_proteome)} GFF protein accessions missing from the full "
            f"proteome FASTA, e.g. {missing_in_proteome[:3]}"
        )

    rows: list[dict] = []
    for record in tqdm(records, desc="Deriving candidate accessions"):
        gene_accs = gene2accs.get(record.gene_id, set())
        matched = sorted(acc for acc in gene_accs if acc_seqs.get(acc) == record.sequence)
        if matched:
            derivation = "gff_gene_sequence_match"
        elif gene_accs:
            derivation = "global_sequence_fallback_gene_mismatch"
            matched = list(seq2accs.get(record.sequence, []))
        else:
            derivation = "global_sequence_fallback_no_gff_gene"
            matched = list(seq2accs.get(record.sequence, []))

        if not matched:
            derivation = "no_accession_match"

        # Integrity check: is the record the longest protein of its own gene?
        lengths = {acc: len(acc_seqs[acc]) for acc in gene_accs if acc in acc_seqs}
        if lengths:
            max_len = max(lengths.values())
            longest_seqs = {acc_seqs[acc] for acc, ln in lengths.items() if ln == max_len}
            longest_confirmed = record.sequence in longest_seqs
        else:
            longest_confirmed = None

        rows.append(
            {
                "record": record,
                "candidate_accessions": matched,
                "accession_derivation": derivation,
                "is_gene_longest_confirmed": longest_confirmed,
            }
        )

    derivation_counts = Counter(row["accession_derivation"] for row in rows)
    logger.info(f"Accession derivation: {dict(derivation_counts)}")
    confirmed = sum(1 for row in rows if row["is_gene_longest_confirmed"] is True)
    logger.info(
        f"Gene-longest-isoform confirmation: {confirmed}/{len(rows)} records "
        "match the longest protein of their gene in the GFF"
    )
    return rows


# --------------------------------------------------------------------------- #
# Step 3: UniProt ID Mapping (official API)
# --------------------------------------------------------------------------- #


def parse_idmapping_stream(tsv_text: str, taxon_id: int) -> tuple[dict[str, set[str]], int, int]:
    """Parse an ID Mapping results stream TSV.

    Returns ({stripped accession: {UniProt accession, ...}},
             n_bombyx_rows, n_non_bombyx_rows).
    """
    mapping: dict[str, set[str]] = defaultdict(set)
    n_kept = 0
    n_dropped = 0
    lines = tsv_text.splitlines()
    if not lines:
        return {}, 0, 0
    header = lines[0].split("\t")
    idx_from = header.index("From")
    idx_entry = header.index("Entry")
    idx_org = header.index("Organism (ID)")
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) <= max(idx_from, idx_entry, idx_org):
            continue
        from_acc = cols[idx_from].strip()
        entry_acc = cols[idx_entry].strip()
        if not from_acc or not entry_acc:
            continue
        if cols[idx_org].strip() == str(taxon_id):
            mapping[from_acc].add(entry_acc)
            n_kept += 1
        else:
            n_dropped += 1
    return dict(mapping), n_kept, n_dropped


def run_id_mapping(
    stripped_accs: list[str],
    chunk_size: int,
    cache_dir: Path,
    taxon_id: int,
    refresh: bool,
) -> tuple[dict[str, set[str]], dict]:
    """Map accessions to UniProtKB via the official ID Mapping API.

    The async ID Mapping backend can be slow for large jobs, so chunks are
    submitted in parallel and then polled together; results are streamed and
    cached per chunk as soon as each job finishes.
    """
    stripped_accs = sorted(stripped_accs)
    chunks = [stripped_accs[i : i + chunk_size] for i in range(0, len(stripped_accs), chunk_size)]
    logger.info(
        f"Submitting {len(stripped_accs)} unique stripped accessions in {len(chunks)} "
        f"job(s) of <= {chunk_size} IDs"
    )
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Phase 1: submit all pending chunks and collect their job ids.
    chunk_jobs: dict[int, str] = {}
    for chunk_idx, chunk in enumerate(chunks):
        cache_path = cache_dir / f"idmapping_chunk_{chunk_idx:02d}.tsv"
        if cache_path.exists() and not refresh:
            logger.info(f"Chunk {chunk_idx}: results already cached")
            continue
        form = urllib.parse.urlencode(
            {"from": IDMAPPING_FROM, "to": IDMAPPING_TO, "ids": ",".join(chunk)}
        ).encode()
        body, _ = http_request(f"{UNIPROT_REST}/idmapping/run", data=form)
        job_id = json.loads(body)["jobId"]
        chunk_jobs[chunk_idx] = job_id
        logger.info(f"ID Mapping job {job_id} submitted for chunk {chunk_idx} ({len(chunk)} ids)")

    # Phase 2: poll all jobs round-robin until every one of them finishes.
    unfinished = dict(chunk_jobs)
    for _round in range(720):  # poll up to ~60 min overall
        if not unfinished:
            break
        for chunk_idx in sorted(unfinished):
            job_id = unfinished[chunk_idx]
            status_body, _ = http_request(
                f"{UNIPROT_REST}/idmapping/status/{job_id}", follow_redirects=False
            )
            status = json.loads(status_body)
            if status.get("jobStatus") == "FINISHED":
                logger.info(f"ID Mapping job {job_id} (chunk {chunk_idx}) finished")
                del unfinished[chunk_idx]
            elif status.get("jobStatus") in ("ERROR",):
                raise RuntimeError(f"ID Mapping job {job_id} failed: {status}")
        if unfinished:
            time.sleep(5)
    if unfinished:
        raise RuntimeError(f"ID Mapping jobs did not finish in ~60 min: {unfinished}")

    mapping: dict[str, set[str]] = defaultdict(set)
    jobs_info: list[dict] = []
    n_dropped_total = 0
    for chunk_idx, chunk in enumerate(tqdm(chunks, desc="Retrieving ID Mapping results")):
        cache_path = cache_dir / f"idmapping_chunk_{chunk_idx:02d}.tsv"
        if cache_path.exists() and not refresh:
            tsv_text = cache_path.read_text()
            job_id = "cached"
        else:
            job_id = chunk_jobs.get(chunk_idx, f"unknown_job_{chunk_idx}")
            stream_url = (
                f"{UNIPROT_REST}/idmapping/uniprotkb/results/stream/{job_id}"
                f"?fields={IDMAPPING_FIELDS}&format=tsv"
            )
            tsv_text = http_request(stream_url, timeout=300)[0].decode()
            cache_path.write_text(tsv_text)

        chunk_map, n_kept, n_dropped = parse_idmapping_stream(tsv_text, taxon_id)
        for from_acc, targets in chunk_map.items():
            mapping[from_acc].update(targets)
        n_dropped_total += n_dropped
        jobs_info.append(
            {
                "chunk": chunk_idx,
                "job_id": job_id,
                "n_ids_submitted": len(chunk),
                "n_bombyx_mapping_rows": n_kept,
                "n_non_bombyx_mapping_rows_dropped": n_dropped,
            }
        )
        logger.info(
            f"Job {chunk_idx} ({job_id}): {len(chunk)} submitted -> {n_kept} Bombyx rows, "
            f"{n_dropped} non-Bombyx rows dropped"
        )

    n_hit = sum(1 for acc in stripped_accs if mapping.get(acc))
    logger.info(
        f"ID Mapping done: {n_hit}/{len(stripped_accs)} accessions map to Bombyx mori "
        f"UniProtKB ({n_dropped_total} non-Bombyx mapping rows discarded)"
    )
    return dict(mapping), {"jobs": jobs_info, "n_non_bombyx_rows_dropped": n_dropped_total}


# --------------------------------------------------------------------------- #
# Step 4: UniProt entry details (keywords with category)
# --------------------------------------------------------------------------- #


def parse_entry(entry: dict) -> dict:
    """Extract reviewed status, protein name and CC keywords from JSON."""
    entry_type = entry.get("entryType", "")
    description = entry.get("proteinDescription", {})
    name = ""
    recommended = description.get("recommendedName", {})
    if recommended.get("fullName", {}).get("value"):
        name = recommended["fullName"]["value"]
    else:
        submissions = description.get("submissionNames", [])
        if submissions and submissions[0].get("fullName", {}).get("value"):
            name = submissions[0]["fullName"]["value"]
    cc_keywords = [
        (kw["id"], kw["name"])
        for kw in entry.get("keywords", [])
        if kw.get("category") == CC_CATEGORY
    ]
    entry_type_lower = entry_type.lower()
    if "unreviewed" in entry_type_lower:
        reviewed_status = "unreviewed"
    elif "reviewed" in entry_type_lower:
        reviewed_status = "reviewed"
    else:
        reviewed_status = "NA"
    return {
        "reviewed": reviewed_status,
        "protein_name": name or "NA",
        "cc": cc_keywords,
        "taxon_id": entry.get("organism", {}).get("taxonId"),
    }


def record_uniprot_release(taxon_id: int) -> None:
    """Record UniProt release headers even when entry responses are cached."""
    if _RELEASE_INFO:
        return
    params = urllib.parse.urlencode(
        {"query": f"organism_id:{taxon_id}", "format": "json", "size": 1}
    )
    http_request(f"{UNIPROT_REST}/uniprotkb/search?{params}", timeout=120)


def fetch_uniprot_entries(
    uniprot_accs: list[str],
    batch_size: int,
    cache_dir: Path,
    taxon_id: int,
    refresh: bool,
) -> tuple[dict[str, dict], int]:
    """Fetch UniProtKB entries in JSON format, in accession batches."""
    accs = sorted(set(uniprot_accs))
    batches = [accs[i : i + batch_size] for i in range(0, len(accs), batch_size)]
    logger.info(f"Fetching details for {len(accs)} UniProt entries in {len(batches)} batches")
    cache_dir.mkdir(parents=True, exist_ok=True)

    entries: dict[str, dict] = {}
    n_taxon_mismatch = 0
    for batch_idx, batch in enumerate(tqdm(batches, desc="Fetching UniProt entries")):
        cache_path = cache_dir / f"uniprot_entries_batch_{batch_idx:04d}.json"
        if cache_path.exists() and not refresh:
            data = json.loads(cache_path.read_text())
        else:
            query = "({}) AND organism_id:{}".format(
                " OR ".join(f"accession:{acc}" for acc in batch), taxon_id
            )
            params = urllib.parse.urlencode({"query": query, "format": "json", "size": 500})
            body, _ = http_request(f"{UNIPROT_REST}/uniprotkb/search?{params}", timeout=120)
            data = json.loads(body)
            cache_path.write_text(json.dumps(data))

        for entry in data.get("results", []):
            parsed = parse_entry(entry)
            if parsed["taxon_id"] not in (None, taxon_id):
                n_taxon_mismatch += 1
                continue
            entries[entry["primaryAccession"]] = parsed
        time.sleep(0.3)

    missing = [acc for acc in accs if acc not in entries]
    if missing:
        logger.warning(f"{len(missing)} entries not returned by search, e.g. {missing[:3]}")
    logger.info(f"Retrieved {len(entries)} entries ({n_taxon_mismatch} taxon mismatches)")
    return entries, len(missing)


# --------------------------------------------------------------------------- #
# Step 5: classification
# --------------------------------------------------------------------------- #


def classify_record(
    candidate_accs: list[str],
    id_mapping: dict[str, set[str]],
    entries: dict[str, dict],
) -> dict:
    """Classify one record according to the CC-keyword rules."""
    uniprot_accs = sorted(
        {
            acc
            for cand in candidate_accs
            for acc in id_mapping.get(split_version(cand)[0], set())
        }
    )
    result = {
        "uniprot_accs": uniprot_accs,
        "n_uniprot": len(uniprot_accs),
        "status": "unmapped",
        "has_cc": False,
        "cc_keywords": [],
        "is_ambiguous": False,
        "unknown_reason": "not_mapped_to_bombyx_uniprotkb",
        "entry_details_missing": 0,
    }
    if not uniprot_accs:
        return result

    cc_sets: dict[str, frozenset[str]] = {}
    for acc in uniprot_accs:
        if acc in entries:
            cc_sets[acc] = frozenset(kw_id for kw_id, _ in entries[acc]["cc"])
        else:
            result["entry_details_missing"] += 1
            cc_sets[acc] = frozenset()
    result["is_ambiguous"] = len(uniprot_accs) > 1
    sets = list(cc_sets.values())
    consistent = all(s == sets[0] for s in sets)

    if len(uniprot_accs) == 1:
        result["status"] = "unique_mapped"
    elif consistent:
        result["status"] = "multi_mapped_consistent"
    else:
        result["status"] = "multi_mapped_inconsistent"

    if result["status"] == "multi_mapped_inconsistent":
        result["unknown_reason"] = "inconsistent_cc_among_multiple_uniprot_entries"
    elif sets[0]:
        result["has_cc"] = True
        result["unknown_reason"] = ""
    else:
        result["unknown_reason"] = "mapped_but_no_cc_keyword"

    if result["has_cc"]:
        cc_by_id = {kw_id: kw_name for acc in uniprot_accs if acc in entries
                    for kw_id, kw_name in entries[acc]["cc"]}
        result["cc_keywords"] = sorted((i, cc_by_id[i]) for i in sets[0])
    return result


# --------------------------------------------------------------------------- #
# Step 6: outputs
# --------------------------------------------------------------------------- #


def join(values: list) -> str:
    """Join list values with ';' for TSV columns."""
    return ";".join(str(v) for v in values)


def build_stats(rows: list[dict]) -> dict:
    """Aggregate the classification counts."""
    total = len(rows)
    mapped = sum(1 for row in rows if row["n_uniprot"] >= 1)
    has_cc = sum(1 for row in rows if row["has_cc"])
    inconsistent = sum(1 for row in rows if row["status"] == "multi_mapped_inconsistent")
    mapped_no_cc = sum(
        1
        for row in rows
        if row["n_uniprot"] >= 1 and not row["has_cc"] and row["status"] != "multi_mapped_inconsistent"
    )
    unmapped = total - mapped
    multi = sum(1 for row in rows if row["n_uniprot"] > 1)
    multi_consistent_cc = sum(1 for row in rows if row["status"] == "multi_mapped_consistent" and row["has_cc"])
    multi_consistent_no_cc = sum(
        1 for row in rows if row["status"] == "multi_mapped_consistent" and not row["has_cc"]
    )
    return {
        "total_records": total,
        "mapped_to_bombyx_uniprotkb": mapped,
        "has_cc_label": has_cc,
        "mapped_but_no_cc_label": mapped_no_cc,
        "unmapped": unmapped,
        "mapping_ambiguous_inconsistent": inconsistent,
        "localization_unknown_total": unmapped + mapped_no_cc + inconsistent,
        "multi_mapped_total": multi,
        "multi_mapped_consistent_with_cc": multi_consistent_cc,
        "multi_mapped_consistent_no_cc": multi_consistent_no_cc,
        "multi_mapped_inconsistent": inconsistent,
        "records_with_multiple_candidate_accessions": sum(
            1 for row in rows if len(row["candidate_accessions"]) > 1
        ),
        "records_with_multiple_uniprot_entries": multi,
    }


def build_duplicate_stats(
    records: list[ProteinRecord],
    rows: list[dict],
    id_mapping: dict[str, set[str]],
) -> dict:
    """Summarize sequence/accession/mapping multiplicities."""
    seq_counts = Counter(record.sequence for record in records)
    gene_counts = Counter(record.gene_id for record in records)
    acc_to_records: Counter[str] = Counter()
    for row in rows:
        for acc in row["candidate_accessions"]:
            acc_to_records[acc] += 1
    uni_to_records: Counter[str] = Counter()
    for row in rows:
        for acc in row["uniprot_accs"]:
            uni_to_records[acc] += 1
    return {
        "duplicate_sequence_groups": sum(1 for count in seq_counts.values() if count > 1),
        "records_in_duplicate_sequence_groups": sum(
            count for count in seq_counts.values() if count > 1
        ),
        "duplicate_gene_ids": sum(1 for count in gene_counts.values() if count > 1),
        "accessions_shared_by_multiple_records": sum(
            1 for count in acc_to_records.values() if count > 1
        ),
        "uniprot_entries_shared_by_multiple_records": sum(
            1 for count in uni_to_records.values() if count > 1
        ),
    }


def write_detail_tables(
    output_dir: Path,
    rows: list[dict],
    id_mapping: dict[str, set[str]],
    entries: dict[str, dict],
) -> None:
    """Write the per-record detail table and the long-format mapping table."""
    detail_path = output_dir / "uniprot_cc_mapping.tsv"
    long_path = output_dir / "uniprot_cc_mapping_long.tsv"

    with open(detail_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(DETAIL_COLUMNS)
        for row in rows:
            record = row["record"]
            primary = row["candidate_accessions"][0] if row["candidate_accessions"] else None
            base, version = split_version(primary) if primary else ("NA", "NA")
            longest_flag = row["is_gene_longest_confirmed"]
            writer.writerow(
                [
                    record.gene_id,
                    record.description,
                    len(record.sequence),
                    base,
                    version,
                    join(row["candidate_accessions"]),
                    len(row["candidate_accessions"]),
                    row["accession_derivation"],
                    "" if longest_flag is None else str(longest_flag).lower(),
                    join(row["uniprot_accs"]),
                    row["n_uniprot"],
                    MAPPING_METHOD if row["candidate_accessions"] else "not_applicable",
                    row["status"],
                    join(entries[acc]["reviewed"] for acc in row["uniprot_accs"] if acc in entries),
                    join(entries[acc]["protein_name"] for acc in row["uniprot_accs"] if acc in entries),
                    join(kw_id for kw_id, _ in row["cc_keywords"]),
                    join(kw_name for _, kw_name in row["cc_keywords"]),
                    str(row["has_cc"]).lower(),
                    str(row["is_ambiguous"]).lower(),
                    row["unknown_reason"],
                    join(ENTRY_URL_TEMPLATE.format(acc=acc) for acc in row["uniprot_accs"]),
                ]
            )

    n_long_rows = 0
    with open(long_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(LONG_COLUMNS)
        for row in rows:
            record = row["record"]
            for cand in row["candidate_accessions"]:
                mapped = sorted(id_mapping.get(split_version(cand)[0], set()))
                if not mapped:
                    writer.writerow(
                        [record.gene_id, *split_version(cand), cand, row["accession_derivation"],
                         "", "", "", "", "", "", row["status"], str(row["has_cc"]).lower()]
                    )
                    n_long_rows += 1
                    continue
                for uni_acc in mapped:
                    entry = entries.get(uni_acc, {})
                    writer.writerow(
                        [
                            record.gene_id,
                            *split_version(cand),
                            cand,
                            row["accession_derivation"],
                            uni_acc,
                            entry.get("reviewed", "NA"),
                            entry.get("protein_name", "NA"),
                            join(kw_id for kw_id, _ in entry.get("cc", [])),
                            join(kw_name for _, kw_name in entry.get("cc", [])),
                            str(bool(entry.get("cc"))).lower(),
                            row["status"],
                            str(row["has_cc"]).lower(),
                        ]
                    )
                    n_long_rows += 1

    logger.info(f"Wrote {detail_path} ({len(rows)} rows)")
    logger.info(f"Wrote {long_path} ({n_long_rows} rows)")


def write_report(output_dir: Path, summary: dict, stats: dict) -> None:
    """Write a human-readable Markdown report (Chinese) with the key table."""
    total = stats["total_records"]

    def pct(n: int) -> str:
        return f"{100 * n / total:.2f}%"

    release = summary["query_metadata"]["uniprot_release"] or "unknown"
    release_date = summary["query_metadata"]["uniprot_release_date"] or "unknown"
    meta = summary["query_metadata"]
    counts = summary["counts"]
    dups = summary["duplicates"]
    lines = [
        "# 家蚕 NCBI 蛋白 → UniProtKB 映射与 CC 定位标签统计",
        "",
        f"- 数据查询日期：{meta['query_date']}",
        f"- UniProt 版本：{release}（{release_date}），来源 {UNIPROT_REST} 官方 REST API",
        f"- ID Mapping：POST `/idmapping/run`，`from={IDMAPPING_FROM}` → `to={IDMAPPING_TO}`（去除版本号后提交）",
        "- 条目详情：GET `/uniprotkb/search`（JSON），关键词类别取自条目 keywords 字段",
        (
            f"- 物种过滤：仅保留 organism_id = {TAXON_ID}（Bombyx mori）条目，其他物种的映射行已剔除"
            f"（共 {meta['id_mapping']['n_non_bombyx_rows_dropped']} 行）"
        ),
        (
            f"- CC 标签定义：UniProt keywords 中 category = '{CC_CATEGORY}' 的关键词"
            "（不使用 GO-CC 或 SUBCELLULAR LOCATION 评论字段，不使用 DeepLoc 预测）"
        ),
        "- reviewed 与 unreviewed 条目均接受",
        "",
        "## 统计结果",
        "",
        "| 统计项目 | 数量 | 占比 |",
        "|---|---:|---:|",
        f"| NCBI蛋白总数 | {total} | 100.00% |",
        (
            f"| 成功映射到家蚕UniProtKB | {stats['mapped_to_bombyx_uniprotkb']} | "
            f"{pct(stats['mapped_to_bombyx_uniprotkb'])} |"
        ),
        f"| 具有至少一个CC标签 | {stats['has_cc_label']} | {pct(stats['has_cc_label'])} |",
        (
            f"| 已映射但无CC标签 | {stats['mapped_but_no_cc_label']} | "
            f"{pct(stats['mapped_but_no_cc_label'])} |"
        ),
        f"| 未成功映射 | {stats['unmapped']} | {pct(stats['unmapped'])} |",
        (
            f"| 映射歧义 | {stats['mapping_ambiguous_inconsistent']} | "
            f"{pct(stats['mapping_ambiguous_inconsistent'])} |"
        ),
        (
            f"| 定位未知总数 | {stats['localization_unknown_total']} | "
            f"{pct(stats['localization_unknown_total'])} |"
        ),
        "",
        "说明：",
        "",
        (
            "- 各行含义按既定规则：`映射歧义` = 多个 UniProt 条目且 CC 标签集不一致"
            "（计入定位未知）；多映射但 CC 标签集一致且非空者计入`具有至少一个CC标签`"
            "并在明细表中保留 `is_ambiguous=true` 标记。"
        ),
        (
            "- 加和关系：`具有至少一个CC标签 + 定位未知总数 = NCBI蛋白总数`；"
            "`成功映射 + 未成功映射 = NCBI蛋白总数`；"
            "`成功映射 = 具有至少一个CC标签 + 已映射但无CC标签 + 映射歧义`。"
        ),
        "",
        "## 补充统计",
        "",
        (
            f"- 多映射记录（对应 >1 个 UniProt 条目）：{stats['multi_mapped_total']}"
            f"（其中标签一致且有 CC：{stats['multi_mapped_consistent_with_cc']}；"
            f"标签一致但无 CC：{stats['multi_mapped_consistent_no_cc']}；"
            f"标签不一致（映射歧义）：{stats['multi_mapped_inconsistent']}）"
        ),
        f"- 候选 accession 数 >1 的记录：{stats['records_with_multiple_candidate_accessions']}",
        f"- accession 推导来源：{counts['accession_derivation_counts']}",
        (
            f"- `longest protein by gene` 校验：{counts['gene_longest_confirmed']}"
            f"/{counts['gene_longest_checkable']} 条记录与其基因在 GFF 中的最长蛋白一致"
            f"（不一致 {counts['gene_longest_conflict']}，"
            f"基因不在 GFF {counts['records_without_gff_gene']}）"
        ),
        (
            f"- 重复序列：{dups['duplicate_sequence_groups']} 组涉及 "
            f"{dups['records_in_duplicate_sequence_groups']} 条记录"
        ),
        f"- 被多条记录共享的 NCBI accession：{dups['accessions_shared_by_multiple_records']}",
        f"- 被多条记录共享的 UniProt 条目：{dups['uniprot_entries_shared_by_multiple_records']}",
        "",
        "## 输出文件",
        "",
        f"- 逐条明细表（每条 NCBI 蛋白一行）：`{output_dir / 'uniprot_cc_mapping.tsv'}`",
        f"- 长表（记录 × accession × UniProt 条目）：`{output_dir / 'uniprot_cc_mapping_long.tsv'}`",
        f"- 汇总 JSON：`{output_dir / 'uniprot_cc_summary.json'}`",
        f"- 原始 API 响应缓存：`{output_dir / 'cache'}`",
    ]
    report_path = output_dir / "report.md"
    report_path.write_text("\n".join(lines) + "\n")
    logger.info(f"Wrote {report_path}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #


@app.command()
def main(
    input_fasta: Path = RAW_DATA_DIR / "GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa",
    full_faa: Path = RAW_DATA_DIR / "protein.faa",
    gff_path: Path = EXTERNAL_DATA_DIR / "ncbi" / "GCF_030269925.1_ASM3026992v2_genomic.gff.gz",
    output_dir: Path = INTERIM_DATA_DIR / "uniprot" / "cc_mapping",
    taxon_id: int = TAXON_ID,
    id_chunk_size: int = 1000,
    entry_batch_size: int = 100,
    refresh_idmapping: bool = typer.Option(
        False, help="Re-run ID Mapping jobs even if cached results exist."
    ),
    refresh_entries: bool = typer.Option(
        False, help="Re-fetch UniProt entries even if cached results exist."
    ),
) -> None:
    """Map NCBI proteins to Bombyx mori UniProtKB and annotate CC keywords."""
    logger.info("Step 1/6: loading inputs")
    records = load_records(input_fasta)
    acc_seqs = load_full_proteome(full_faa)

    logger.info("Step 2/6: gene -> protein accession links from GFF3")
    download_gff_if_needed(gff_path)
    gene2accs = parse_gff_gene_proteins(gff_path)

    logger.info("Step 3/6: attaching candidate accessions to records")
    rows = derive_candidates(records, acc_seqs, gene2accs)

    stripped_to_versioned: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        for acc in row["candidate_accessions"]:
            stripped_to_versioned[split_version(acc)[0]].append(acc)
    all_stripped = sorted(stripped_to_versioned)
    logger.info(f"Step 4/6: UniProt ID Mapping for {len(all_stripped)} unique accessions")
    id_mapping, id_job_info = run_id_mapping(
        all_stripped, id_chunk_size, output_dir / "cache", taxon_id, refresh_idmapping
    )

    mapped_accs = sorted({acc for targets in id_mapping.values() for acc in targets})
    logger.info("Step 5/6: fetching UniProt entries and CC keywords")
    entries, n_entries_missing = fetch_uniprot_entries(
        mapped_accs, entry_batch_size, output_dir / "cache", taxon_id, refresh_entries
    )
    record_uniprot_release(taxon_id)

    logger.info("Step 6/6: classifying records and writing outputs")
    output_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        row.update(classify_record(row["candidate_accessions"], id_mapping, entries))

    stats = build_stats(rows)
    dup_stats = build_duplicate_stats(records, rows, id_mapping)
    derivation_counts = Counter(row["accession_derivation"] for row in rows)
    longest_confirmed = sum(1 for row in rows if row["is_gene_longest_confirmed"] is True)
    longest_conflict = sum(1 for row in rows if row["is_gene_longest_confirmed"] is False)
    gene_checkable = longest_confirmed + longest_conflict

    summary = {
        "query_metadata": {
            "query_date": datetime.now(timezone.utc).date().isoformat(),
            "query_date_note": "UTC date of the UniProt REST API queries",
            "uniprot_release": _RELEASE_INFO.get("X-UniProt-Release"),
            "uniprot_release_date": _RELEASE_INFO.get("X-UniProt-Release-Date"),
            "id_mapping": {
                "api": "POST https://rest.uniprot.org/idmapping/run",
                "from": IDMAPPING_FROM,
                "to": IDMAPPING_TO,
                "poll": "GET /idmapping/status/{jobId}",
                "results": "GET /idmapping/uniprotkb/results/stream/{jobId}",
                "fields": IDMAPPING_FIELDS,
                "chunk_size": id_chunk_size,
                **id_job_info,
            },
            "entry_details": {
                "api": "GET https://rest.uniprot.org/uniprotkb/search",
                "format": "json",
                "batch_size": entry_batch_size,
                "n_entries_not_returned": n_entries_missing,
            },
            "taxon_filter": f"organism_id == {taxon_id}",
            "cc_definition": f"UniProt keywords with category == '{CC_CATEGORY}'",
            "reviewed_policy": "reviewed and unreviewed entries both accepted",
        },
        "inputs": {
            "working_fasta": {
                "path": str(input_fasta),
                "records": len(records),
                "sha256": sha256_of_file(input_fasta),
            },
            "full_proteome": {"path": str(full_faa), "sha256": sha256_of_file(full_faa)},
            "gff": {"path": str(gff_path), "sha256": sha256_of_file(gff_path)},
        },
        "counts": {
            **stats,
            "accession_derivation_counts": dict(derivation_counts),
            "gene_longest_confirmed": longest_confirmed,
            "gene_longest_conflict": longest_conflict,
            "gene_longest_checkable": gene_checkable,
            "records_without_gff_gene": derivation_counts.get(
                "global_sequence_fallback_no_gff_gene", 0
            ),
            "records_with_missing_entry_details": sum(
                1 for row in rows if row["entry_details_missing"] > 0
            ),
        },
        "duplicates": dup_stats,
        "verification": {},
    }

    # Verify the counting identities before writing outputs.
    identities = {
        "has_cc_plus_unknown_equals_total": (
            stats["has_cc_label"] + stats["localization_unknown_total"] == len(records)
        ),
        "mapped_plus_unmapped_equals_total": (
            stats["mapped_to_bombyx_uniprotkb"] + stats["unmapped"] == len(records)
        ),
        "mapped_partition": (
            stats["mapped_to_bombyx_uniprotkb"]
            == stats["has_cc_label"]
            + stats["mapped_but_no_cc_label"]
            + stats["mapping_ambiguous_inconsistent"]
        ),
        "rows_equal_total": len(rows) == len(records),
    }
    summary["verification"] = identities
    if not all(identities.values()):
        raise RuntimeError(f"Counting identity check failed: {identities}")
    logger.info(f"All counting identities verified: {identities}")

    write_detail_tables(output_dir, rows, id_mapping, entries)
    summary_path = output_dir / "uniprot_cc_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    logger.info(f"Wrote {summary_path}")
    write_report(output_dir, summary, stats)

    logger.success(
        f"Done. mapped={stats['mapped_to_bombyx_uniprotkb']}, "
        f"has_cc={stats['has_cc_label']}, "
        f"unknown={stats['localization_unknown_total']} / {len(records)}"
    )


if __name__ == "__main__":
    app()
