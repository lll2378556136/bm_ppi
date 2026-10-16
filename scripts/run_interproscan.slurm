#!/usr/bin/env bash
#SBATCH --job-name=iprscan_bm_ppi
#SBATCH --output=/share/home/lll/projects/bm_ppi/logs/interproscan_%A_%a.out
#SBATCH --error=/share/home/lll/projects/bm_ppi/logs/interproscan_%A_%a.err
#SBATCH --partition=gpu,compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00

# Run InterProScan 5.78-109.0 (functional domain annotation) on the Bombyx
# mori longest-protein-per-gene set through the official container image,
# converted to a singularity SIF. One SLURM array task processes one chunk
# FASTA produced by bm_ppi/prepare_interproscan_input.py.
#
# Expected input layout (prepared beforehand):
#   data/interim/interproscan/input/chunk_0001.faa ... chunk_NNNN.faa
#   data/interim/interproscan/input/manifest.tsv
#
# Outputs per chunk (written to OUTPUT_DIR):
#   chunk_NNNN.tsv   primary result table
#   chunk_NNNN.json  full structured results (GO terms, pathways)
#
# Submit with:
#   sbatch --array=1-27 scripts/run_interproscan.slurm
# (chunk count comes from the manifest / prepare step output).
#
# Environment overrides:
#   SIF                singularity image (default sifs/interproscan-5.78-109.0.sif)
#   IPR_DATA           extracted interproscan data dir mounted to /opt/interproscan/data
#   DISABLE_PRECALC=1  pass -dp (disable EBI match lookup); set 0 if the
#                      compute nodes can reach the internet (much faster)
#   FORCE=1            rerun a chunk even if its .tsv already exists

set -euo pipefail

# Slurm copies this script into its spool directory before execution, so the
# project root cannot be derived from BASH_SOURCE (same convention as
# run_dpam.slurm / deeploc_accurate.slurm).
PROJECT_ROOT="/share/home/lll/projects/bm_ppi"

# ---- Configuration (each value can be overridden via environment) ----
IPR_VERSION="${IPR_VERSION:-5.78-109.0}"
SIF="${SIF:-${PROJECT_ROOT}/sifs/interproscan-${IPR_VERSION}.sif}"
IPR_DATA="${IPR_DATA:-${PROJECT_ROOT}/data/external/interproscan/interproscan-${IPR_VERSION}/data}"
INPUT_DIR="${INPUT_DIR:-${PROJECT_ROOT}/data/interim/interproscan/input}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_ROOT}/data/interim/interproscan/output}"
TMP_ROOT="${TMP_ROOT:-${PROJECT_ROOT}/scratch/interproscan_tmp}"
IPR_BIN="${IPR_BIN:-/opt/interproscan/interproscan.sh}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
DISABLE_PRECALC="${DISABLE_PRECALC:-1}"

CHUNK_NAME="chunk_$(printf '%04d' "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set; submit with --array}")"
CHUNK_FASTA="${INPUT_DIR}/${CHUNK_NAME}.faa"
OUT_BASE="${OUTPUT_DIR}/${CHUNK_NAME}"
TMP_DIR="${TMP_ROOT}/${CHUNK_NAME}"

# ---- Validation ----
command -v singularity >/dev/null 2>&1 || { echo "Error: singularity not found." >&2; exit 1; }

if [[ ! -f "${SIF}" ]]; then
    echo "Error: InterProScan singularity image not found: ${SIF}" >&2
    echo "Hint: cd ${PROJECT_ROOT}/sifs && singularity pull interproscan-${IPR_VERSION}.sif docker://interpro/interproscan:${IPR_VERSION}" >&2
    exit 1
fi

if [[ ! -d "${IPR_DATA}" ]]; then
    echo "Error: InterProScan data directory not found: ${IPR_DATA}" >&2
    echo "Hint: download and extract interproscan-data-${IPR_VERSION}.tar.gz under data/external/interproscan/" >&2
    exit 1
fi

if [[ ! -s "${CHUNK_FASTA}" ]]; then
    echo "Error: chunk FASTA missing or empty: ${CHUNK_FASTA}" >&2
    exit 1
fi

# Skip a finished chunk unless FORCE=1 is set.
if [[ -s "${OUT_BASE}.tsv" && "${FORCE:-0}" != "1" ]]; then
    echo "Output ${OUT_BASE}.tsv already exists; skipping (set FORCE=1 to rerun)."
    exit 0
fi

mkdir -p "${OUTPUT_DIR}" "${TMP_DIR}"

N_SEQS=$(grep -c '^>' "${CHUNK_FASTA}")
PRECALC_FLAG=()
if [[ "${DISABLE_PRECALC}" == "1" ]]; then
    PRECALC_FLAG=(-dp)
fi

echo "Starting InterProScan ${IPR_VERSION} for ${CHUNK_NAME} at $(date --iso-8601=seconds)"
echo "Node: ${SLURMD_NODENAME:-unknown} (job ${SLURM_ARRAY_JOB_ID:-?}_${SLURM_ARRAY_TASK_ID})"
echo "Input: ${CHUNK_FASTA} (${N_SEQS} sequences)"
echo "Output base: ${OUT_BASE}"
echo "Image: ${SIF}"
echo "Data: ${IPR_DATA}"
echo "Threads: ${THREADS}; disable-precalc: ${DISABLE_PRECALC}"

# ---- Run InterProScan inside the container ----
# Data is mounted read-only at /opt/interproscan/data (the image expects it
# there); input/output/tmp dirs are mounted read-write at their real paths.
singularity exec --fakeroot \
    --bind "${IPR_DATA}:/opt/interproscan/data:ro" \
    --bind "${INPUT_DIR}:${INPUT_DIR}:rw" \
    --bind "${OUTPUT_DIR}:${OUTPUT_DIR}:rw" \
    --bind "${TMP_DIR}:${TMP_DIR}:rw" \
    "${SIF}" \
    "${IPR_BIN}" \
        --input "${CHUNK_FASTA}" \
        --seqtype p \
        --output-file-base "${OUT_BASE}" \
        --formats TSV,JSON \
        --goterms \
        --iprlookup \
        --pathways \
        --cpu "${THREADS}" \
        --tempdir "${TMP_DIR}" \
        "${PRECALC_FLAG[@]}"

# ---- Check outputs ----
if [[ ! -s "${OUT_BASE}.tsv" && ! -s "${OUT_BASE}.json" ]]; then
    echo "Error: InterProScan produced no output for ${CHUNK_NAME}." >&2
    exit 1
fi

# A chunk can legitimately have zero annotation hits; an empty TSV with a
# non-empty run log is acceptable, so only require the files to exist.
N_HITS=0
if [[ -s "${OUT_BASE}.tsv" ]]; then
    N_HITS=$(wc -l < "${OUT_BASE}.tsv")
fi
N_PROTEINS=0
if [[ -s "${OUT_BASE}.tsv" ]]; then
    N_PROTEINS=$(cut -f1 "${OUT_BASE}.tsv" | sort -u | wc -l)
fi

echo "InterProScan finished for ${CHUNK_NAME} at $(date --iso-8601=seconds)"
echo "Summary: ${N_HITS} annotation lines covering ${N_PROTEINS}/${N_SEQS} proteins"
