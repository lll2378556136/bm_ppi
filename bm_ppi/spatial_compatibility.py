"""Pre-filter protein pairs by spatial compatibility from DeepLoc 2.1 predictions.

Implements references/DeepLoc2.1_空间相容性预筛选方案.md (v1.0).  Every protein
pair is assigned one of five space classes, in the fixed priority order

    DIRECT > PROBABILITY_RESCUE > INTERFACE_RESCUE > UNCERTAIN > INCOMPATIBLE

and only INCOMPATIBLE pairs are discarded.  Classification runs at the
whole-protein level; fragments inherit their parent's class downstream when the
pair table is joined with the DPAM domain table.  L_core and L_border are
derived from the raw DeepLoc probabilities with the per-label thresholds of the
Accurate (ProtT5) model, not from the reported "Localizations" column (which
contains a fallback label when nothing passes the threshold).

Usage:
    pixi run python -m bm_ppi.spatial_compatibility --help
"""

from __future__ import annotations

import contextlib
import csv
import gzip
import json
from pathlib import Path

from loguru import logger
import numpy as np
from tqdm import tqdm
import typer

from bm_ppi.config import INTERIM_DATA_DIR
from bm_ppi.merge_localization import DEEPLOC_PROBABILITY_COLUMNS, normalize_id

app = typer.Typer()

LOCALIZATION_DIR = INTERIM_DATA_DIR / "localization"
DEFAULT_DEEPLOC_CSV = LOCALIZATION_DIR / "deeploc_2.1_accurate_prott5.csv"
DEFAULT_OUTPUT_DIR = LOCALIZATION_DIR / "spatial_compatibility"

# DeepLoc 2.1 localization labels, in the CSV column order.
LOCALIZATION_LABELS = [
    "Cytoplasm",
    "Nucleus",
    "Extracellular",
    "Cell membrane",
    "Mitochondrion",
    "Plastid",
    "Endoplasmic reticulum",
    "Lysosome/Vacuole",
    "Golgi apparatus",
    "Peroxisome",
]
CYT, NUC, EXT, CM = range(4)  # label indices used by the interface rules

# Per-label decision thresholds of the Accurate (ProtT5) model: label_threshold[1:]
# from scratch/deeploc21_src/deeploc2_package/DeepLoc2/deeploc2.py, line 109.
# DeepLoc itself compares `p > tau` strictly (utils.convert_label2string); this
# module follows the scheme's `p >= tau`, which differs only on exact float
# equality and therefore never triggers in practice.
ACCURATE_LABEL_THRESHOLDS = np.array(
    [
        0.47612305,  # Cytoplasm
        0.50136719,  # Nucleus
        0.61728516,  # Extracellular
        0.56464844,  # Cell membrane
        0.62197266,  # Mitochondrion
        0.63945312,  # Plastid
        0.60898438,  # Endoplasmic reticulum
        0.58476562,  # Lysosome/Vacuole
        0.64941406,  # Golgi apparatus
        0.73642578,  # Peroxisome
    ]
)

# Membrane-type labels and their Accurate-model thresholds (deeploc2.py, line 111).
MEMBRANE_TYPES = ["Peripheral", "Transmembrane", "Lipid anchor", "Soluble"]
ACCURATE_MEMTYPE_THRESHOLDS = np.array([0.60, 0.51, 0.82, 0.50])

(
    CLASS_DIRECT,
    CLASS_PROBABILITY_RESCUE,
    CLASS_INTERFACE_RESCUE,
    CLASS_UNCERTAIN,
    CLASS_INCOMPATIBLE,
) = range(5)
CLASS_NAMES = [
    "DIRECT",
    "PROBABILITY_RESCUE",
    "INTERFACE_RESCUE",
    "UNCERTAIN",
    "INCOMPATIBLE",
]

PAIR_COLUMNS = [
    "protein_A",
    "protein_B",
    "shared_core_labels",
    "probability_rescue_label",
    "interface_type",
    "space_class",
    "space_score",
    "decision",
    "reason",
]

# Precomputed "a|b" string for every 10-bit core-label mask.
MASK_STRINGS = [
    "|".join(label for k, label in enumerate(LOCALIZATION_LABELS) if (mask >> k) & 1)
    for mask in range(1 << len(LOCALIZATION_LABELS))
]


def read_deeploc_csv(path: Path) -> dict:
    """Parse a DeepLoc 2.1 result CSV into label lists and probability arrays."""
    required = ["Protein_ID", "Localizations", "Signals", "Membrane types"]
    required += LOCALIZATION_LABELS + MEMBRANE_TYPES
    ids: list[str] = []
    loc_probs: list[list[float]] = []
    mem_probs: list[list[float]] = []
    reported: list[list[str]] = []
    seen: set[str] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [c for c in required if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"DeepLoc CSV {path} is missing columns: {missing}")
        for row in reader:
            protein_id = normalize_id(row.get("Protein_ID", ""))
            if not protein_id:
                raise ValueError(f"Empty Protein_ID in {path}")
            if protein_id in seen:
                raise ValueError(f"Duplicate protein ID in DeepLoc CSV: {protein_id}")
            seen.add(protein_id)
            ids.append(protein_id)
            loc_probs.append([float(row[c]) for c in LOCALIZATION_LABELS])
            mem_probs.append([float(row[c]) for c in MEMBRANE_TYPES])
            reported.append(
                [
                    row.get("Localizations", ""),
                    row.get("Signals", ""),
                    row.get("Membrane types", ""),
                ]
            )
    if not ids:
        raise ValueError(f"No protein rows found in {path}")
    return {
        "ids": ids,
        "probs": np.asarray(loc_probs, dtype=np.float32),
        "memtype_probs": np.asarray(mem_probs, dtype=np.float32),
        "reported_localizations": [r[0] for r in reported],
        "reported_signals": [r[1] for r in reported],
        "reported_membrane_types": [r[2] for r in reported],
    }


def subset_table(table: dict, indices: list[int]) -> dict:
    """Subset every array/list of a DeepLoc table by row indices."""
    out = {}
    for key, value in table.items():
        out[key] = value[indices] if isinstance(value, np.ndarray) else [value[i] for i in indices]
    return out


def read_master_evidence(path: Path) -> dict[str, tuple[bool, bool]]:
    """Read optional (secretory, transmembrane) evidence flags keyed by protein ID.

    Expected columns come from merge_localization.py output: `protein_id`,
    `is_secretory_pathway`, `is_transmembrane`.
    """
    evidence: dict[str, tuple[bool, bool]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            protein_id = normalize_id(row.get("protein_id", ""))
            if not protein_id:
                continue
            truthy = {"true", "1", "yes"}
            evidence[protein_id] = (
                row.get("is_secretory_pathway", "").strip().lower() in truthy,
                row.get("is_transmembrane", "").strip().lower() in truthy,
            )
    return evidence


def read_fragment_counts(path: Path) -> dict[str, int]:
    """Count DPAM fragments per protein; every Range segment counts as one fragment.

    Assumes the DPAM domains schema (Protein/Range/...); Judge-based filtering is
    left to downstream pairing decisions.
    """
    counts: dict[str, int] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            protein_id = normalize_id(row.get("Protein", ""))
            ranges = row.get("Range", "")
            if not protein_id or not ranges.strip():
                continue
            segments = [seg for seg in ranges.split(",") if seg.strip()]
            counts[protein_id] = counts.get(protein_id, 0) + len(segments)
    return counts


def read_subset_ids(path: Path) -> set[str]:
    """Read one protein ID per line from a subset file."""
    with path.open(encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def derive_matrices(
    table: dict,
    tau: np.ndarray,
    delta: np.ndarray,
    ambiguity_margin: float,
    weak_margin: float,
    max_core_labels: int,
    master_evidence: dict[str, tuple[bool, bool]],
) -> dict:
    """Compute L_core/L_border masks, rescue scores, and per-protein flags.

    A protein is UNCERTAIN (scheme section 4.3) when it has no core label, a
    single weak core label, near-tied top probabilities, or too many core
    labels; such proteins never end up in INCOMPATIBLE pairs.
    """
    probs = table["probs"]
    lower = (tau - delta)[None, :]
    core = probs >= tau[None, :]
    border = (probs >= lower) & ~core
    q = np.clip((probs - lower) / delta[None, :], 0.0, 1.0).astype(np.float32)

    mem_probs = table["memtype_probs"]
    mem_interface = (
        (mem_probs[:, 0] >= ACCURATE_MEMTYPE_THRESHOLDS[0])
        | (mem_probs[:, 1] >= ACCURATE_MEMTYPE_THRESHOLDS[1])
        | (mem_probs[:, 2] >= ACCURATE_MEMTYPE_THRESHOLDS[2])
    )
    signal_peptide = np.array(
        ["Signal peptide" in s.split("|") for s in table["reported_signals"]], dtype=bool
    )
    for i, protein_id in enumerate(table["ids"]):
        secretory, transmembrane = master_evidence.get(protein_id, (False, False))
        signal_peptide[i] |= secretory
        mem_interface[i] |= transmembrane

    core_count = core.sum(axis=1)
    top2 = np.sort(probs, axis=1)[:, -2:]
    near_tied = (top2[:, 1] - top2[:, 0]) < ambiguity_margin
    only_core_prob = np.max(np.where(core, probs, -1.0), axis=1)
    only_label = np.argmax(core, axis=1)
    single_weak = (core_count == 1) & (only_core_prob < tau[only_label] + weak_margin)
    if max_core_labels > 0:
        multi_core = core_count >= max_core_labels
    else:
        multi_core = np.zeros(len(table["ids"]), dtype=bool)
    uncertain = (core_count == 0) | near_tied | single_weak | multi_core

    reasons: list[str] = []
    for i in range(len(table["ids"])):
        parts = []
        if core_count[i] == 0:
            parts.append("no_core_label")
        if near_tied[i]:
            parts.append("near_tied_top_labels")
        if single_weak[i]:
            parts.append("single_weak_core_label")
        if multi_core[i]:
            parts.append("multi_core_labels")
        reasons.append("|".join(parts))

    shifts = 1 << np.arange(len(LOCALIZATION_LABELS))
    return {
        "core": core,
        "border": border,
        "q": q,
        "core_mask_int": (core * shifts).sum(axis=1).astype(np.int16),
        "cm_core": core[:, CM],
        "cyt_core": core[:, CYT],
        "ext_core": core[:, EXT],
        "mem_interface": mem_interface,
        "signal_peptide": signal_peptide,
        "uncertain": uncertain,
        "uncertain_reasons": reasons,
    }


def open_tsv(path: Path, gzipped: bool):
    """Open a TSV file for writing, transparently gzipping when requested."""
    if gzipped:
        return gzip.open(path, "wt", newline="", encoding="utf-8")
    return path.open("w", newline="", encoding="utf-8")


def _uncertain_pair_reason(reason_a: str, reason_b: str) -> str:
    """Format a pair-level uncertainty reason from per-protein reasons."""
    parts = []
    if reason_a:
        parts.append(f"A:{reason_a}")
    if reason_b:
        parts.append(f"B:{reason_b}")
    return "uncertain(" + ";".join(parts) + ")"


def classify_pairs_and_write(
    table: dict,
    mats: dict,
    theta_probability: float,
    block_size: int,
    outdir: Path,
    gzip_output: bool,
    emit_incompatible: bool,
    frag_vector: np.ndarray | None,
) -> dict:
    """Classify every protein pair by blocks and stream kept pairs to TSV."""
    ids = table["ids"]
    n = len(ids)
    ids_np = np.asarray(ids, dtype=object)
    probs = table["probs"]
    q = mats["q"]
    core_mask_int = mats["core_mask_int"]
    cm_mem = mats["cm_core"] & mats["mem_interface"]
    ext_evidence = mats["ext_core"] & mats["signal_peptide"]
    cyt_core = mats["cyt_core"]
    uncertain = mats["uncertain"]
    prot_reasons = np.asarray(mats["uncertain_reasons"], dtype=object)
    label_arr = np.array(LOCALIZATION_LABELS, dtype=object)
    class_arr = np.array(CLASS_NAMES, dtype=object)

    n_pairs_total = n * (n - 1) // 2
    class_pair_counts = np.zeros(5, dtype=np.int64)
    frag_pair_counts = np.zeros(5, dtype=np.float64)

    pair_path = outdir / (
        "pair_spatial_compatibility.tsv.gz" if gzip_output else "pair_spatial_compatibility.tsv"
    )
    discarded_path = outdir / (
        "pair_spatial_compatibility_discarded.tsv.gz"
        if gzip_output
        else "pair_spatial_compatibility_discarded.tsv"
    )

    with contextlib.ExitStack() as stack:
        pair_handle = stack.enter_context(open_tsv(pair_path, gzip_output))
        writer = csv.writer(pair_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(PAIR_COLUMNS)
        inc_writer = None
        if emit_incompatible:
            inc_handle = stack.enter_context(open_tsv(discarded_path, gzip_output))
            inc_writer = csv.writer(inc_handle, delimiter="\t", lineterminator="\n")
            inc_writer.writerow(PAIR_COLUMNS)

        for i0 in tqdm(range(0, n - 1, block_size), desc="Classifying protein pairs"):
            rows = np.arange(i0, min(i0 + block_size, n))

            # Step 1 DIRECT: shared high-confidence labels via bitwise AND of core masks.
            shared = core_mask_int[rows][:, None] & core_mask_int[None, :]
            direct = shared != 0

            # Step 2 PROBABILITY_RESCUE: S = max_k min(q_A(k), q_B(k)) over the border zone.
            s_best = np.zeros((rows.size, n), dtype=np.float32)
            k_best = np.full((rows.size, n), -1, dtype=np.int8)
            for k in range(len(LOCALIZATION_LABELS)):
                m_k = np.minimum(q[rows, k][:, None], q[None, :, k])
                upd = m_k > s_best
                s_best[upd] = m_k[upd]
                k_best[upd] = k
            prob_rescue = s_best >= theta_probability

            # Step 3 INTERFACE_RESCUE: membrane-cytoplasm / membrane-extracellular contact.
            iface_mc = (cm_mem[rows][:, None] & cyt_core[None, :]) | (
                cyt_core[rows][:, None] & cm_mem[None, :]
            )
            iface_me = (cm_mem[rows][:, None] & ext_evidence[None, :]) | (
                ext_evidence[rows][:, None] & cm_mem[None, :]
            )

            # Steps 4/5: assign in increasing priority so higher classes overwrite.
            cls = np.full((rows.size, n), CLASS_INCOMPATIBLE, dtype=np.int8)
            cls[uncertain[rows][:, None] | uncertain[None, :]] = CLASS_UNCERTAIN
            cls[iface_mc | iface_me] = CLASS_INTERFACE_RESCUE
            cls[prob_rescue] = CLASS_PROBABILITY_RESCUE
            cls[direct] = CLASS_DIRECT

            valid = np.arange(n)[None, :] > rows[:, None]
            cls[~valid] = -1

            ra, ca = np.nonzero(valid)
            vcls = cls[ra, ca]
            pos = vcls >= 0
            vcls = vcls[pos]
            class_pair_counts += np.bincount(vcls, minlength=5)
            if frag_vector is not None:
                frag_pair_counts += np.bincount(
                    vcls,
                    weights=frag_vector[i0 + ra[pos]] * frag_vector[ca[pos]],
                    minlength=5,
                )

            ka, kb = np.nonzero((cls >= 0) & (cls <= CLASS_UNCERTAIN))
            m = ka.size
            if m:
                ga, gb = i0 + ka, kb
                cls_sel = cls[ka, kb]
                shared_vals = shared[ka, kb]
                prot_a = ids_np[ga].tolist()
                prot_b = ids_np[gb].tolist()
                shared_col = [MASK_STRINGS[v] for v in shared_vals.tolist()]

                kb_sel = k_best[ka, kb]
                rescue_col = np.where(
                    cls_sel == CLASS_PROBABILITY_RESCUE,
                    np.where(kb_sel >= 0, label_arr[np.maximum(kb_sel, 0)], ""),
                    "",
                )
                mc_sel = iface_mc[ka, kb]
                me_sel = iface_me[ka, kb]
                iface_col = np.where(
                    cls_sel == CLASS_INTERFACE_RESCUE,
                    np.where(
                        mc_sel & me_sel,
                        "membrane_cytoplasm|membrane_extracellular",
                        np.where(
                            mc_sel,
                            "membrane_cytoplasm",
                            np.where(me_sel, "membrane_extracellular", ""),
                        ),
                    ),
                    "",
                )

                scores = np.full(m, np.nan, dtype=np.float32)
                dsel = np.flatnonzero(cls_sel == CLASS_DIRECT)
                if dsel.size:
                    # Space score: max over shared core labels of min(p_A, p_B).
                    best = np.full(dsel.size, -1.0, dtype=np.float32)
                    sv = shared_vals[dsel]
                    pa, pb = probs[ga[dsel]], probs[gb[dsel]]
                    for k in range(len(LOCALIZATION_LABELS)):
                        on = ((sv >> k) & 1) > 0
                        if on.any():
                            vals = np.minimum(pa[on, k], pb[on, k])
                            best[on] = np.maximum(best[on], vals)
                    scores[dsel] = best
                psel = np.flatnonzero(cls_sel == CLASS_PROBABILITY_RESCUE)
                if psel.size:
                    scores[psel] = s_best[ka[psel], kb[psel]]
                isel = np.flatnonzero(cls_sel == CLASS_INTERFACE_RESCUE)
                if isel.size:
                    a_idx, b_idx = ga[isel], gb[isel]
                    cand = np.full((isel.size, 4), -1.0, dtype=np.float32)
                    sel = cm_mem[a_idx] & cyt_core[b_idx]
                    cand[sel, 0] = np.minimum(probs[a_idx[sel], CM], probs[b_idx[sel], CYT])
                    sel = cm_mem[b_idx] & cyt_core[a_idx]
                    cand[sel, 1] = np.minimum(probs[b_idx[sel], CM], probs[a_idx[sel], CYT])
                    sel = cm_mem[a_idx] & ext_evidence[b_idx]
                    cand[sel, 2] = np.minimum(probs[a_idx[sel], CM], probs[b_idx[sel], EXT])
                    sel = cm_mem[b_idx] & ext_evidence[a_idx]
                    cand[sel, 3] = np.minimum(probs[b_idx[sel], CM], probs[a_idx[sel], EXT])
                    scores[isel] = cand.max(axis=1)

                reasons = np.full(m, "", dtype=object)
                reasons[dsel] = "shared high-confidence localization"
                reasons[psel] = "near-threshold common label"
                reasons[isel] = "membrane-interface contact possible"
                usel = np.flatnonzero(cls_sel == CLASS_UNCERTAIN)
                if usel.size:
                    reasons[usel] = [
                        _uncertain_pair_reason(a, b)
                        for a, b in zip(prot_reasons[ga[usel]], prot_reasons[gb[usel]])
                    ]

                score_col = ["" if not np.isfinite(s) else f"{s:.4f}" for s in scores.tolist()]
                writer.writerows(
                    (a, b, sh, r, i_, c, sc, "KEEP", rs)
                    for a, b, sh, r, i_, c, sc, rs in zip(
                        prot_a,
                        prot_b,
                        shared_col,
                        rescue_col.tolist(),
                        iface_col.tolist(),
                        class_arr[cls_sel].tolist(),
                        score_col,
                        reasons.tolist(),
                    )
                )

            if inc_writer is not None:
                xa, xb = np.nonzero(cls == CLASS_INCOMPATIBLE)
                if xa.size:
                    s_vals = s_best[xa, xb]
                    inc_writer.writerows(
                        (a, b, "", "", "", "INCOMPATIBLE", "", "DISCARD", reason)
                        for a, b, sv, reason in zip(
                            ids_np[i0 + xa].tolist(),
                            ids_np[xb].tolist(),
                            s_vals.tolist(),
                            (
                                f"no core overlap; S_probability={sv:.3f} < {theta_probability}; "
                                "no membrane interface; both sides confident"
                                for sv in s_vals.tolist()
                            ),
                        )
                    )

    counted = int(class_pair_counts.sum())
    if counted != n_pairs_total:
        raise AssertionError(f"Classified {counted} pairs, expected {n_pairs_total}")
    stats = {
        "n_pairs_total": n_pairs_total,
        "pair_counts": {name: int(class_pair_counts[i]) for i, name in enumerate(CLASS_NAMES)},
        "n_pairs_kept": int(class_pair_counts[:CLASS_UNCERTAIN].sum()),
        "n_pairs_discarded": int(class_pair_counts[CLASS_INCOMPATIBLE]),
        "pair_reduction_rate": float(class_pair_counts[CLASS_INCOMPATIBLE] / n_pairs_total),
    }
    if frag_vector is not None:
        total_frag = float(frag_pair_counts.sum())
        kept_frag = float(frag_pair_counts[:CLASS_UNCERTAIN].sum())
        stats["fragment"] = {
            "n_fragment_pairs_total": total_frag,
            "n_fragment_pairs_kept": kept_frag,
            "fragment_pair_reduction_rate": float(1.0 - kept_frag / total_frag)
            if total_frag
            else None,
        }
    stats["outputs"] = {
        "pair_table": str(pair_path),
        "discarded_table": str(discarded_path) if emit_incompatible else None,
    }
    return stats


def write_protein_table(path: Path, table: dict, mats: dict) -> None:
    """Write per-protein localization evidence and derived label sets."""
    probs = table["probs"]
    mem_probs = table["memtype_probs"]
    border = mats["border"]
    columns = [
        "protein_id",
        "deeploc_localizations",
        "deeploc_signals",
        "deeploc_membrane_types",
        "loc_core",
        "loc_border",
        "has_signal_peptide",
        "has_interface_membrane_type",
        "is_uncertain",
        "uncertain_reason",
        *[DEEPLOC_PROBABILITY_COLUMNS[label] for label in LOCALIZATION_LABELS],
        *[DEEPLOC_PROBABILITY_COLUMNS[memtype] for memtype in MEMBRANE_TYPES],
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for i, protein_id in enumerate(table["ids"]):
            writer.writerow(
                [
                    protein_id,
                    table["reported_localizations"][i],
                    table["reported_signals"][i],
                    table["reported_membrane_types"][i],
                    MASK_STRINGS[mats["core_mask_int"][i]],
                    "|".join(label for k, label in enumerate(LOCALIZATION_LABELS) if border[i, k]),
                    bool(mats["signal_peptide"][i]),
                    bool(mats["mem_interface"][i]),
                    bool(mats["uncertain"][i]),
                    mats["uncertain_reasons"][i],
                    *[f"{probs[i, k]:.6f}" for k in range(len(LOCALIZATION_LABELS))],
                    *[f"{mem_probs[i, k]:.6f}" for k in range(len(MEMBRANE_TYPES))],
                ]
            )


def write_summary(path: Path, summary: dict) -> None:
    """Write the run summary as JSON."""
    with path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


@app.command()
def main(
    deeploc_csv: Path = DEFAULT_DEEPLOC_CSV,
    outdir: Path = DEFAULT_OUTPUT_DIR,
    master_tsv: Path | None = None,
    domains_table: Path | None = None,
    subset: Path | None = None,
    border_delta: float = 0.10,
    theta_probability: float = 0.50,
    ambiguity_margin: float = 0.05,
    weak_margin: float = 0.10,
    max_core_labels: int = 3,
    block_size: int = 256,
    emit_incompatible: bool = False,
    gzip_output: bool = False,
):
    """Classify protein-pair spatial compatibility from DeepLoc 2.1 results."""
    if not 0.0 < theta_probability < 1.0:
        raise typer.BadParameter("theta_probability must lie in (0, 1)")
    if border_delta <= 0:
        raise typer.BadParameter("border_delta must be positive")
    if ambiguity_margin < 0 or weak_margin < 0:
        raise typer.BadParameter("margins must be non-negative")
    if block_size <= 0:
        raise typer.BadParameter("block_size must be positive")

    table = read_deeploc_csv(deeploc_csv)
    logger.info(f"Loaded {len(table['ids'])} proteins from {deeploc_csv}")
    if subset is not None:
        wanted = read_subset_ids(subset)
        indices = [i for i, protein_id in enumerate(table["ids"]) if protein_id in wanted]
        logger.info(f"Subset filter: kept {len(indices)}/{len(table['ids'])} proteins")
        table = subset_table(table, indices)
    if len(table["ids"]) < 2:
        raise ValueError("Spatial compatibility classification needs at least two proteins")

    master_evidence = read_master_evidence(master_tsv) if master_tsv else {}
    if master_tsv:
        matched = sum(protein_id in master_evidence for protein_id in table["ids"])
        logger.info(f"Master localization evidence matched {matched} proteins")

    frag_vector = None
    if domains_table is not None:
        frag_counts = read_fragment_counts(domains_table)
        frag_vector = np.array([frag_counts.get(p, 0) for p in table["ids"]], dtype=np.float64)
        covered = int((frag_vector > 0).sum())
        logger.info(f"Fragment table covers {covered}/{len(table['ids'])} proteins")
        if covered == 0:
            frag_vector = None

    tau = ACCURATE_LABEL_THRESHOLDS.astype(np.float64)
    delta = np.full(len(LOCALIZATION_LABELS), border_delta, dtype=np.float64)
    mats = derive_matrices(
        table, tau, delta, ambiguity_margin, weak_margin, max_core_labels, master_evidence
    )
    n_uncertain = int(mats["uncertain"].sum())
    logger.info(f"Uncertain localization: {n_uncertain}/{len(table['ids'])} proteins")

    outdir.mkdir(parents=True, exist_ok=True)
    stats = classify_pairs_and_write(
        table,
        mats,
        theta_probability,
        block_size,
        outdir,
        gzip_output,
        emit_incompatible,
        frag_vector,
    )

    protein_path = outdir / "protein_localization.tsv"
    write_protein_table(protein_path, table, mats)
    summary_path = outdir / "spatial_compatibility_summary.json"
    write_summary(
        summary_path,
        {
            "n_proteins": len(table["ids"]),
            "n_uncertain_proteins": n_uncertain,
            "parameters": {
                "model": "DeepLoc 2.1 Accurate (ProtT5)",
                "label_thresholds": dict(
                    zip(LOCALIZATION_LABELS, ACCURATE_LABEL_THRESHOLDS.tolist())
                ),
                "membrane_type_thresholds": dict(
                    zip(MEMBRANE_TYPES, ACCURATE_MEMTYPE_THRESHOLDS.tolist())
                ),
                "border_delta": border_delta,
                "theta_probability": theta_probability,
                "ambiguity_margin": ambiguity_margin,
                "weak_margin": weak_margin,
                "max_core_labels": max_core_labels,
            },
            "inputs": {
                "deeploc_csv": str(deeploc_csv),
                "master_tsv": str(master_tsv) if master_tsv else None,
                "domains_table": str(domains_table) if domains_table else None,
                "subset": str(subset) if subset else None,
            },
            # Known-interaction retention and false-exclusion rate need a curated
            # Bombyx mori PPI reference set and are calibrated later (scheme section 9).
            "known_interaction_retention": None,
            "false_exclusion_rate": None,
            **stats,
            "outputs": {
                **stats["outputs"],
                "protein_table": str(protein_path),
                "summary": str(summary_path),
            },
        },
    )

    for name, count in stats["pair_counts"].items():
        logger.info(f"  {name}: {count} pairs")
    logger.info(
        f"Discarded protein pairs: {stats['n_pairs_discarded']} "
        f"({stats['pair_reduction_rate']:.2%} of {stats['n_pairs_total']})"
    )
    if "fragment" in stats and stats["fragment"]["fragment_pair_reduction_rate"] is not None:
        logger.info(
            f"Fragment-pair reduction rate: "
            f"{stats['fragment']['fragment_pair_reduction_rate']:.2%}"
        )
    logger.success(f"Spatial compatibility pre-filtering complete: {outdir}")


if __name__ == "__main__":
    app()
