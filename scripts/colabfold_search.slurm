#!/usr/bin/env bash
#SBATCH --job-name=colabfold_search
#SBATCH --output=/share/home/lll/projects/bm_ppi/logs/colabfold_search_%j.out
#SBATCH --error=/share/home/lll/projects/bm_ppi/logs/colabfold_search_%j.err
#SBATCH --partition=gpu
#SBATCH --nodelist=node4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"
WORK_DIR="${PROJECT_ROOT}"
SIF="${PROJECT_ROOT}/sifs/colabfold_1.6.2-cuda12.sif"
CACHE_DIR="${PROJECT_ROOT}/cache"
INPUT_FASTA="${WORK_DIR}/data/raw/GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa"
OUTPUT_DIR="${WORK_DIR}/data/processed/colabfold_20260831"
DB_DIR="/data/ssd/colabfold_db"

mkdir -p "${PROJECT_ROOT}/logs" "${CACHE_DIR}" "${OUTPUT_DIR}"

if [[ ! -f "${SIF}" ]]; then
    echo "Error: Singularity image not found: ${SIF}" >&2
    exit 1
fi

if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "Error: Input FASTA not found: ${INPUT_FASTA}" >&2
    exit 1
fi

if [[ ! -d "${DB_DIR}" ]]; then
    echo "Error: ColabFold database not found: ${DB_DIR}" >&2
    exit 1
fi

echo "Starting ColabFold search at $(date --iso-8601=seconds)"
echo "Input: ${INPUT_FASTA}"
echo "Database: ${DB_DIR}"
echo "Output: ${OUTPUT_DIR}"

# Use mmseqs inside the container; bind node-local SSD DB and project paths.
singularity run --nv \
    -B "${CACHE_DIR}:/cache" \
    -B "${WORK_DIR}:/work" \
    -B /data/ssd:/data/ssd \
    "${SIF}" \
    colabfold_search --gpu 1 \
        "/work/data/raw/$(basename "${INPUT_FASTA}")" \
        /data/ssd/colabfold_db \
        /work/data/processed/colabfold_20260831

echo "ColabFold search completed at $(date --iso-8601=seconds)"
