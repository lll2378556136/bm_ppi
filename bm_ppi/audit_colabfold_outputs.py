"""Audit ColabFold structure-prediction outputs for the B. mori proteome.

Scans every ``data/processed/colabfold_output_gcf_*`` directory (main
batches plus ``*_5m`` per-model retry dirs) and checks, per protein,
whether the expected five ranked models and their companion files were
produced:

- ``<name>_unrelaxed_rank_00N_*.pdb``  ranks 1-5 (the five models)
- ``<name>_scores_rank_00N_*.json``    per-model confidence scores
- ``<name>_predicted_aligned_error_v1.json``  PAE matrix
- ``<name>.done.txt``                  ColabFold completion marker

Outputs ``data/interim/colabfold_audit.tsv`` (one row per protein) and
logs per-directory and proteome-level summaries.

Run with:

    pixi run python -m bm_ppi.audit_colabfold_outputs
"""

from __future__ import annotations

from collections import defaultdict
import csv
from pathlib import Path
import re

from Bio import SeqIO
from loguru import logger
import typer

from bm_ppi.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR, RAW_DATA_DIR

app = typer.Typer()

PDB_RE = re.compile(r"^(?P<name>.+)_unrelaxed_rank_(?P<rank>\d{3})_.+\.pdb$")
SCORES_RE = re.compile(r"^(?P<name>.+)_scores_rank_(?P<rank>\d{3})_.+\.json$")
PAE_SUFFIX = "_predicted_aligned_error_v1.json"
DONE_SUFFIX = ".done.txt"
A3M_SUFFIX = ".a3m"


def collect(dirpath: Path, proteins: dict[str, dict]) -> None:
    """Record per-protein artifacts found in one output directory."""
    for f in dirpath.iterdir():
        if not f.is_file():
            continue
        fname = f.name
        m = PDB_RE.match(fname)
        if m:
            rec = proteins[m["name"]]
            rec["dirs"].add(dirpath.name)
            rec["pdb_ranks"].add(m["rank"])
            if f.stat().st_size == 0:
                rec["empty_files"].append(fname)
            continue
        m = SCORES_RE.match(fname)
        if m:
            rec = proteins[m["name"]]
            rec["dirs"].add(dirpath.name)
            rec["scores_ranks"].add(m["rank"])
            if f.stat().st_size == 0:
                rec["empty_files"].append(fname)
            continue
        if fname.endswith(PAE_SUFFIX):
            rec = proteins[fname[: -len(PAE_SUFFIX)]]
            rec["dirs"].add(dirpath.name)
            rec["has_pae"] = True
        elif fname.endswith(DONE_SUFFIX):
            rec = proteins[fname[: -len(DONE_SUFFIX)]]
            rec["dirs"].add(dirpath.name)
            rec["has_done"] = True
        elif fname.endswith(A3M_SUFFIX):
            rec = proteins[fname[: -len(A3M_SUFFIX)]]
            rec["dirs"].add(dirpath.name)
            rec["has_a3m"] = True


@app.command()
def main(
    input_fasta: Path = RAW_DATA_DIR / "GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa",
    output_glob: str = "colabfold_output_gcf_*",
    report_path: Path = INTERIM_DATA_DIR / "colabfold_audit.tsv",
) -> None:
    # Map protein id -> sanitized query name. ColabFold replaces runs of
    # characters outside [A-Za-z0-9._-] with a single underscore; verify the
    # reconstruction against real .a3m names below.
    id_to_name: dict[str, str] = {}
    for rec in SeqIO.parse(str(input_fasta), "fasta"):
        name = re.sub(r"[^A-Za-z0-9._-]", "_", rec.description)
        id_to_name[rec.id] = name
    logger.info(f"{len(id_to_name)} expected proteins from {input_fasta.name}")

    dirs = sorted(d for d in PROCESSED_DATA_DIR.glob(output_glob) if d.is_dir())
    logger.info(f"Scanning {len(dirs)} output directories")

    proteins: dict[str, dict] = defaultdict(
        lambda: {
            "dirs": set(),
            "pdb_ranks": set(),
            "scores_ranks": set(),
            "has_pae": False,
            "has_done": False,
            "has_a3m": False,
            "empty_files": [],
        }
    )
    per_dir = {}
    for d in dirs:
        before = sum(1 for v in proteins.values() if v["pdb_ranks"])
        collect(d, proteins)
        after = sum(1 for v in proteins.values() if v["pdb_ranks"])
        per_dir[d.name] = after - before

    name_to_id = {v: k for k, v in id_to_name.items()}
    unexpected = sorted(n for n in proteins if n not in name_to_id)
    if unexpected:
        logger.warning(
            f"{len(unexpected)} output names do not match any FASTA header, e.g. {unexpected[:5]}"
        )

    # ---- Per-protein status ----
    rows = []
    status_counts: dict[str, int] = defaultdict(int)
    for seq_id, name in id_to_name.items():
        rec = proteins.get(name)
        if rec is None or not rec["pdb_ranks"]:
            status = "absent" if rec is None else "no_models"
            missing = "1-5"
            dirs_seen = ";".join(sorted(rec["dirs"])) if rec else ""
        else:
            missing_ranks = sorted(set("12345") - {r.lstrip("0") for r in rec["pdb_ranks"]})
            missing = ",".join(missing_ranks)
            if missing_ranks:
                status = "partial_models"
            elif rec["empty_files"]:
                status = "empty_file"
            elif not rec["has_pae"] or not rec["has_done"]:
                status = "missing_companion"
            else:
                status = "complete"
            dirs_seen = ";".join(sorted(rec["dirs"]))
        status_counts[status] += 1
        rows.append(
            [
                seq_id,
                name,
                status,
                len(rec["pdb_ranks"]) if rec else 0,
                len(rec["scores_ranks"]) if rec else 0,
                int(rec["has_pae"]) if rec else 0,
                int(rec["has_done"]) if rec else 0,
                missing,
                dirs_seen,
            ]
        )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "protein_id",
                "query_name",
                "status",
                "n_pdb_models",
                "n_scores_json",
                "has_pae_json",
                "has_done_txt",
                "missing_ranks",
                "output_dirs",
            ]
        )
        writer.writerows(rows)

    logger.info("Per-directory protein counts (with >=1 model):")
    for dname, n in per_dir.items():
        logger.info(f"  {dname}: {n}")
    logger.info("Proteome-level status:")
    for status, n in sorted(status_counts.items()):
        logger.info(f"  {status}: {n}")

    problem = [r for r in rows if r[2] != "complete"]
    for r in problem[:20]:
        logger.info(f"  {r[0]}\t{r[2]}\tmissing={r[7]}\tdirs={r[8]}")
    logger.success(f"Report written to {report_path}")


if __name__ == "__main__":
    app()
