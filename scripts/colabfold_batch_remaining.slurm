#!/usr/bin/env bash
#SBATCH --job-name=colabfold_remaining
#SBATCH --output=logs/colabfold_remaining_%j.out
#SBATCH --error=logs/colabfold_remaining_%j.err
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=110G
#SBATCH --gres=gpu:1
#SBATCH --time=11-00:00:00

# Resume worker for the "remaining" (unfinished, ultra-long) queries of each
# batch. Inputs come from the *_remaining directories, while results are still
# written into the original per-batch output directories so that ColabFold's
# per-query .done.txt markers keep the full history consistent.

set -euo pipefail

# Workers must use the project root exported by the launcher because
# BASH_SOURCE points to Slurm's spool copy inside an allocated job.
if [[ "${COLABFOLD_WORKER:-}" == "1" ]]; then
    if [[ -z "${PROJECT_ROOT:-}" ]]; then
        echo "Error: PROJECT_ROOT was not exported by the launcher." >&2
        exit 1
    fi

    PROJECT_ROOT="$(realpath -- "${PROJECT_ROOT}")"
    SCRIPT_PATH="${PROJECT_ROOT}/scripts/colabfold_batch_remaining.slurm"
else
    # Resolve paths from the original script on the login node.
    SCRIPT_PATH="$(realpath -- "${BASH_SOURCE[0]}")"
    SCRIPT_DIR="$(dirname -- "${SCRIPT_PATH}")"
    PROJECT_ROOT="$(realpath -- "${SCRIPT_DIR}/..")"
fi

# Launcher mode submits four independent workers with fixed node assignments.
if [[ "${1:-}" == "submit" ]]; then
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "Error: Submit mode must run outside a Slurm job." >&2
        exit 1
    fi

    mkdir -p "${PROJECT_ROOT}/logs"

    for batch_id in 001 002 003 004; do
        case "${batch_id}" in
            001|002) node="node3" ;;
            003|004) node="node4" ;;
        esac

        sbatch \
            --export="ALL,PROJECT_ROOT=${PROJECT_ROOT},COLABFOLD_WORKER=1,BATCH_ID=${batch_id}" \
            --job-name="cf_rem_${batch_id}" \
            --nodelist="${node}" \
            --chdir="${PROJECT_ROOT}" \
            --output="${PROJECT_ROOT}/logs/colabfold_remaining_${batch_id}_%j.out" \
            --error="${PROJECT_ROOT}/logs/colabfold_remaining_${batch_id}_%j.err" \
            "${SCRIPT_PATH}"
    done

    echo "Submitted four ColabFold remaining-sequence resume jobs."
    exit 0
fi

# Prevent accidental submission of the worker without using launcher mode.
if [[ "${COLABFOLD_WORKER:-}" != "1" ]]; then
    echo "Usage: bash ${SCRIPT_PATH} submit" >&2
    exit 1
fi

# Validate the batch ID strictly.
BATCH_ID="${BATCH_ID:-}"
if [[ ! "${BATCH_ID}" =~ ^(001|002|003|004)$ ]]; then
    echo "Error: BATCH_ID must be one of 001|002|003|004, got '${BATCH_ID}'" >&2
    exit 1
fi

SIF="${PROJECT_ROOT}/sifs/colabfold_1.6.2-cuda12.sif"
CACHE_DIR="/data/ssd/cache"
# Only the unfinished queries of the batch; the original input directories
# are never touched by this script.
INPUT_DIR="${PROJECT_ROOT}/data/interim/colabfold_batch${BATCH_ID}_remaining"
# Results still land in the original per-batch output directory. ColabFold
# skips queries whose <jobname>.done.txt marker already exists there, so
# completed predictions are never recomputed or overwritten.
OUTPUT_DIR="${PROJECT_ROOT}/data/processed/colabfold_output_gcf_${BATCH_ID}"

mkdir -p "${PROJECT_ROOT}/logs" "${OUTPUT_DIR}"

if [[ ! -f "${SIF}" ]]; then
    echo "Error: Singularity image not found: ${SIF}" >&2
    exit 1
fi

if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "Error: Input directory not found: ${INPUT_DIR}" >&2
    exit 1
fi

# Require at least one .a3m file (regular file or symlink) in the input directory.
A3M_COUNT=$(find "${INPUT_DIR}" -maxdepth 1 \( -type f -o -type l \) -name '*.a3m' | wc -l)
if [[ "${A3M_COUNT}" -eq 0 ]]; then
    echo "Error: No .a3m files found in ${INPUT_DIR}" >&2
    exit 1
fi

# Confirm the AlphaFold/ColabFold model parameters are present.
if [[ ! -d "${CACHE_DIR}/colabfold/params" ]]; then
    echo "Error: ColabFold params directory not found: ${CACHE_DIR}/colabfold/params" >&2
    exit 1
fi

echo "Starting ColabFold remaining-sequence resume for batch ${BATCH_ID} at $(date --iso-8601=seconds)"
echo "Node: ${SLURM_JOB_NODELIST:-unknown}"
echo "Input MSAs: ${INPUT_DIR} (${A3M_COUNT} files)"
echo "Output: ${OUTPUT_DIR}"
echo "Cache: ${CACHE_DIR}"

# Strip the project root prefix so paths are relative to the /work mount.
INPUT_IN_WORK="/work/${INPUT_DIR#"${PROJECT_ROOT}/"}"
OUTPUT_IN_WORK="/work/${OUTPUT_DIR#"${PROJECT_ROOT}/"}"

# Bind cache (AlphaFold2 weights) and project (MSAs/output) into the container.
# Memory/compute-saving flags for ultra-long queries (all supported by
# colabfold_batch 1.6.2, verified via `colabfold_batch --help`):
#   --num-recycle 1     single recycle pass (default: 3)
#   --num-models 1      only model_1 instead of all 5 models
#   --stop-at-score 85  stop early once pLDDT/pTM exceeds 85
singularity run --nv \
    -B "${CACHE_DIR}:/cache" \
    -B "${PROJECT_ROOT}:/work" \
    "${SIF}" \
    colabfold_batch \
      --model-type alphafold2_ptm \
      --num-recycle 1 \
      --num-models 1 \
      --stop-at-score 85 \
      "${INPUT_IN_WORK}" \
      "${OUTPUT_IN_WORK}"

echo "ColabFold remaining-sequence resume for batch ${BATCH_ID} completed at $(date --iso-8601=seconds)"
