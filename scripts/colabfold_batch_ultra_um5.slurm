#!/usr/bin/env bash
#SBATCH --job-name=colabfold_ultra_um5
#SBATCH --output=logs/colabfold_ultra_um5_%j.out
#SBATCH --error=logs/colabfold_ultra_um5_%j.err
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=110G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00

# Full 5-model rerun for ultra-long sequences with unified memory enabled.
#
# This variant re-predicts the 8 sequences just below 5000 aa that previously
# finished the single-model (--num-models 1) unified-memory resume run. It
# writes to a fresh, independent output directory (*_5m) so the original
# per-batch results are never touched, and keeps JAX/ColabFold's default
# unified-memory behaviour so GPU allocations can spill into host RAM.

set -euo pipefail

# Workers must use the project root exported by the launcher because
# BASH_SOURCE points to Slurm's spool copy inside an allocated job.
if [[ "${COLABFOLD_WORKER:-}" == "1" ]]; then
    if [[ -z "${PROJECT_ROOT:-}" ]]; then
        echo "Error: PROJECT_ROOT was not exported by the launcher." >&2
        exit 1
    fi

    PROJECT_ROOT="$(realpath -- "${PROJECT_ROOT}")"
    SCRIPT_PATH="${PROJECT_ROOT}/scripts/colabfold_batch_ultra_um5.slurm"
else
    # Resolve paths from the original script on the login node.
    SCRIPT_PATH="$(realpath -- "${BASH_SOURCE[0]}")"
    SCRIPT_DIR="$(dirname -- "${SCRIPT_PATH}")"
    PROJECT_ROOT="$(realpath -- "${SCRIPT_DIR}/..")"
fi

# Prevent accidental submission of the worker without using launcher mode.
if [[ "${COLABFOLD_WORKER:-}" != "1" ]]; then
    echo "Usage: export COLABFOLD_WORKER=1 BATCH_ID=001 GROUP_ID=1 PROJECT_ROOT=... && sbatch ${SCRIPT_PATH}" >&2
    exit 1
fi

# Validate the batch and group IDs strictly.
BATCH_ID="${BATCH_ID:-}"
GROUP_ID="${GROUP_ID:-}"
if [[ ! "${BATCH_ID}" =~ ^(001|002|003|004)$ ]]; then
    echo "Error: BATCH_ID must be one of 001|002|003|004, got '${BATCH_ID}'" >&2
    exit 1
fi
if [[ ! "${GROUP_ID}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: GROUP_ID must be a positive integer, got '${GROUP_ID}'" >&2
    exit 1
fi

SIF="${PROJECT_ROOT}/sifs/colabfold_1.6.2-cuda12.sif"
CACHE_DIR="/data/ssd/cache"
# Same input directories as the unified-memory resume run; symlinks point back
# to the original *_remaining directories, so no large A3M files are duplicated
# and the inputs are only read, never modified.
INPUT_DIR="${PROJECT_ROOT}/data/interim/colabfold_batch${BATCH_ID}_remaining_um${GROUP_ID}"
# Fresh, independent output directory for the 5-model rerun. ColabFold skips
# queries whose <jobname>.done.txt marker already exists in the output dir, so
# using a new directory guarantees the predictions are recomputed from scratch
# while the original per-batch results stay untouched.
OUTPUT_DIR="${PROJECT_ROOT}/data/processed/colabfold_output_gcf_${BATCH_ID}_5m"

mkdir -p "${PROJECT_ROOT}/logs" "${OUTPUT_DIR}"

if [[ ! -f "${SIF}" ]]; then
    echo "Error: Singularity image not found: ${SIF}" >&2
    exit 1
fi

if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "Error: Input directory not found: ${INPUT_DIR}" >&2
    exit 1
fi

# Require exactly one .a3m file (regular file or symlink) in the input directory.
A3M_COUNT=$(find "${INPUT_DIR}" -maxdepth 1 \( -type f -o -type l \) -name '*.a3m' | wc -l)
if [[ "${A3M_COUNT}" -ne 1 ]]; then
    echo "Error: Expected exactly 1 .a3m file in ${INPUT_DIR}, found ${A3M_COUNT}" >&2
    exit 1
fi

# Confirm the AlphaFold/ColabFold model parameters are present.
if [[ ! -d "${CACHE_DIR}/colabfold/params" ]]; then
    echo "Error: ColabFold params directory not found: ${CACHE_DIR}/colabfold/params" >&2
    exit 1
fi

echo "Starting ColabFold ultra-long 5-model rerun (unified memory) for batch ${BATCH_ID} group ${GROUP_ID} at $(date --iso-8601=seconds)"
echo "Node: ${SLURM_JOB_NODELIST:-unknown}"
echo "Input MSAs: ${INPUT_DIR} (${A3M_COUNT} files)"
echo "Output: ${OUTPUT_DIR}"
echo "Cache: ${CACHE_DIR}"

# Strip the project root prefix so paths are relative to the /work mount.
INPUT_IN_WORK="/work/${INPUT_DIR#"${PROJECT_ROOT}/"}"
OUTPUT_IN_WORK="/work/${OUTPUT_DIR#"${PROJECT_ROOT}/"}"

# Bind cache (AlphaFold2 weights) and project (MSAs/output) into the container.
# Memory/compute-saving flags for ultra-long queries:
#   --max-msa 256:512     shrink the MSA stack (default 512:1024)
#   --num-recycle 1       single recycle pass (default: 3)
#   --num-models 5        full 5-model prediction (model_1..model_5)
# No --stop-at-score so all 5 models always run to completion.
# Unified memory is intentionally left enabled so the GPU can spill into host RAM.
singularity run --nv \
    -B "${CACHE_DIR}:/cache" \
    -B "${PROJECT_ROOT}:/work" \
    "${SIF}" \
    colabfold_batch \
      --model-type alphafold2_ptm \
      --max-msa 256:512 \
      --num-recycle 1 \
      --num-models 5 \
      "${INPUT_IN_WORK}" \
      "${OUTPUT_IN_WORK}"

echo "ColabFold ultra-long 5-model rerun (unified memory) for batch ${BATCH_ID} group ${GROUP_ID} completed at $(date --iso-8601=seconds)"
