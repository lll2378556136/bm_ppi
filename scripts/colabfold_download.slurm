#!/usr/bin/env bash
#SBATCH --job-name=colabfold_download
#SBATCH --output=/share/home/gyk/projects/bm_ppi/logs/colabfold_download_%j.out
#SBATCH --error=/share/home/gyk/projects/bm_ppi/logs/colabfold_download_%j.err
#SBATCH --partition=gpu
#SBATCH --nodelist=node4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --gres=gpu:1
#SBATCH --time=5:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"
SIF="${PROJECT_ROOT}/sifs/colabfold_1.6.2-cuda12.sif"
CACHE_DIR="${PROJECT_ROOT}/cache"

mkdir -p "${PROJECT_ROOT}/logs" "${CACHE_DIR}"

if [[ ! -f "${SIF}" ]]; then
    echo "Error: Singularity image not found: ${SIF}" >&2
    exit 1
fi

echo "Starting ColabFold parameter download at $(date --iso-8601=seconds)"
echo "SIF: ${SIF}"
echo "Cache: ${CACHE_DIR}"

# Download AlphaFold/ColabFold model parameters into the bound /cache directory.
singularity run --nv \
    -B "${CACHE_DIR}:/cache" \
    "${SIF}" \
    python -m colabfold.download

echo "ColabFold parameter download completed at $(date --iso-8601=seconds)"
