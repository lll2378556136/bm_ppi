#!/usr/bin/env bash
#SBATCH --job-name=colabfold_um5_by_model_lowmem
#SBATCH --output=logs/colabfold_um5_by_model_lowmem_%j.out
#SBATCH --error=logs/colabfold_um5_by_model_lowmem_%j.err
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=110G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00

# Low-memory rerun of model_2 only, for an ultra-long query sequence.
# The previous model_2 attempt with --max-msa 256:512 peaked around 119 GiB
# of host memory and was SIGKILLed by the cgroup limit during the final
# stage. This variant halves the MSA depth (--max-msa 128:256) to shrink the
# dominant memory consumer; everything else mirrors
# colabfold_batch_ultra_um5_by_model.slurm.
#
# Required environment variables:
#   PROJECT_ROOT       absolute path to the bm_ppi repository
#   COLABFOLD_WORKER   must be set to 1 when submitting the worker
#   BATCH_ID           one of 001|002|003|004
#   GROUP_ID           positive integer matching the *_remaining_um* directory

set -euo pipefail

if [[ "${COLABFOLD_WORKER:-}" == "1" ]]; then
    if [[ -z "${PROJECT_ROOT:-}" ]]; then
        echo "Error: PROJECT_ROOT was not exported by the launcher." >&2
        exit 1
    fi

    PROJECT_ROOT="$(realpath -- "${PROJECT_ROOT}")"
    SCRIPT_PATH="${PROJECT_ROOT}/scripts/colabfold_batch_ultra_um5_by_model_lowmem.slurm"
else
    SCRIPT_PATH="$(realpath -- "${BASH_SOURCE[0]}")"
    SCRIPT_DIR="$(dirname -- "${SCRIPT_PATH}")"
    PROJECT_ROOT="$(realpath -- "${SCRIPT_DIR}/..")"
fi

if [[ "${COLABFOLD_WORKER:-}" != "1" ]]; then
    echo "Usage: export COLABFOLD_WORKER=1 BATCH_ID=001 GROUP_ID=1 PROJECT_ROOT=... && sbatch ${SCRIPT_PATH}" >&2
    exit 1
fi

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

# This lowmem variant exists specifically to rerun model_2, so the model is
# fixed here instead of being read from the environment.
MODEL_ID="2"

SIF="${PROJECT_ROOT}/sifs/colabfold_1.6.2-cuda12.sif"
CACHE_DIR="/data/ssd/cache"
INPUT_DIR="${PROJECT_ROOT}/data/interim/colabfold_batch${BATCH_ID}_remaining_um${GROUP_ID}"

# Each model gets its own output directory. ColabFold writes a .done.txt
# marker per query, so sharing one directory between separate single-model
# jobs would make later models skip as "already completed". The model_2
# directory currently has no .done.txt, so this rerun will process the query
# again and fill in the missing model_2 result files.
OUTPUT_DIR="${PROJECT_ROOT}/data/processed/colabfold_output_gcf_${BATCH_ID}_5m/model_${MODEL_ID}"

mkdir -p "${PROJECT_ROOT}/logs" "${OUTPUT_DIR}"

if [[ ! -f "${SIF}" ]]; then
    echo "Error: Singularity image not found: ${SIF}" >&2
    exit 1
fi
if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "Error: Input directory not found: ${INPUT_DIR}" >&2
    exit 1
fi

A3M_COUNT=$(find "${INPUT_DIR}" -maxdepth 1 \( -type f -o -type l \) -name '*.a3m' | wc -l)
if [[ "${A3M_COUNT}" -ne 1 ]]; then
    echo "Error: Expected exactly 1 .a3m file in ${INPUT_DIR}, found ${A3M_COUNT}" >&2
    exit 1
fi
if [[ ! -d "${CACHE_DIR}/colabfold/params" ]]; then
    echo "Error: ColabFold params directory not found: ${CACHE_DIR}/colabfold/params" >&2
    exit 1
fi

echo "Starting ColabFold low-memory model_2 rerun: batch ${BATCH_ID} group ${GROUP_ID} model_${MODEL_ID} at $(date --iso-8601=seconds)"
echo "Node: ${SLURM_JOB_NODELIST:-unknown}"
echo "Input MSAs: ${INPUT_DIR} (${A3M_COUNT} files)"
echo "Output: ${OUTPUT_DIR}"
echo "Cache: ${CACHE_DIR}"

INPUT_IN_WORK="/work/${INPUT_DIR#"${PROJECT_ROOT}/"}"
OUTPUT_IN_WORK="/work/${OUTPUT_DIR#"${PROJECT_ROOT}/"}"

# Memory/compute-saving flags:
#   --max-msa 128:256   halve the MSA depth vs the killed 256:512 run; the MSA
#                      stack dominates host RAM for a 4928 aa query, so this is
#                      the main knob pushing peak memory below the 110G cap
#   --num-recycle 1     single recycle pass
#   --num-models 1     run exactly one model in this process
#   --model-order 2    rerun model_2 only
# Unified memory is intentionally left enabled (no --disable-unified-memory).
singularity run --nv \
    -B "${CACHE_DIR}:/cache" \
    -B "${PROJECT_ROOT}:/work" \
    "${SIF}" \
    colabfold_batch \
      --model-type alphafold2_ptm \
      --max-msa 128:256 \
      --num-recycle 1 \
      --num-models 1 \
      --model-order 2 \
      "${INPUT_IN_WORK}" \
      "${OUTPUT_IN_WORK}"

echo "ColabFold low-memory model_2 rerun completed at $(date --iso-8601=seconds)"
