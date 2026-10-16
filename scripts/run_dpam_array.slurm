#!/usr/bin/env bash
#SBATCH --job-name=dpam_array
#SBATCH --output=/share/home/lll/projects/bm_ppi/logs/dpam_array_%A_%a.out
#SBATCH --error=/share/home/lll/projects/bm_ppi/logs/dpam_array_%A_%a.err
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00

# Submit with an explicit bounded array, for example:
# sbatch --array=0-63%2 scripts/run_dpam_array.slurm

set -euo pipefail

PROJECT_ROOT="/share/home/lll/projects/bm_ppi"
CHUNKS_TSV="${CHUNKS_TSV:-${PROJECT_ROOT}/data/interim/DPAM_runs/full_v1/chunks.tsv}"
ARRAY_INDEX="${SLURM_ARRAY_TASK_ID:-}"

[[ -n "${ARRAY_INDEX}" ]] || {
    echo "Error: SLURM_ARRAY_TASK_ID is not set; submit this script as an array." >&2
    exit 1
}
[[ "${ARRAY_INDEX}" =~ ^[0-9]+$ ]] || {
    echo "Error: invalid array index: ${ARRAY_INDEX}" >&2
    exit 1
}
[[ -s "${CHUNKS_TSV}" ]] || {
    echo "Error: chunk manifest missing or empty: ${CHUNKS_TSV}" >&2
    exit 1
}

MATCH_COUNT=$(awk -F'\t' -v array_index="${ARRAY_INDEX}" 'NR > 1 && $1 == array_index {count++} END {print count + 0}' "${CHUNKS_TSV}")
if [[ "${MATCH_COUNT}" -ne 1 ]]; then
    echo "Error: expected one manifest row for array index ${ARRAY_INDEX}, found ${MATCH_COUNT}." >&2
    exit 1
fi

ROW=$(awk -F'\t' -v array_index="${ARRAY_INDEX}" 'NR > 1 && $1 == array_index {print; exit}' "${CHUNKS_TSV}")
IFS=$'\t' read -r manifest_index chunk_id dataset input_dir structure_list n_targets total_aa max_length total_weight <<<"${ROW}"

[[ "${manifest_index}" == "${ARRAY_INDEX}" ]] || {
    echo "Error: parsed manifest index does not match requested array index." >&2
    exit 1
}
[[ -d "${input_dir}/${dataset}" && -s "${structure_list}" ]] || {
    echo "Error: chunk input is incomplete for ${dataset} under ${input_dir}." >&2
    exit 1
}

VALIDATION_JSON="${input_dir}/validation.json"
if [[ -s "${VALIDATION_JSON}" ]] && grep -q '"hard_success": true' "${VALIDATION_JSON}"; then
    echo "Chunk ${dataset} already has a successful validation; skipping."
    exit 0
fi

echo "Array index: ${ARRAY_INDEX}; chunk: ${chunk_id}; dataset: ${dataset}"
echo "Targets: ${n_targets}; total aa: ${total_aa}; max length: ${max_length}"

export DATASET="${dataset}"
export INPUT_DIR="${input_dir}"
export RESUME="${RESUME:-0}"
exec bash "${PROJECT_ROOT}/scripts/run_dpam.slurm"
