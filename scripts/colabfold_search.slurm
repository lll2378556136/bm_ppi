#!/usr/bin/env bash
#SBATCH --job-name=colabfold_search
#SBATCH --output=colabfold_search_%j.out
#SBATCH --error=colabfold_search_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=128G
#SBATCH --gres=gpu2:1
#SBATCH --time=24:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"
WORK_DIR="${PROJECT_ROOT}"
SIF="${PROJECT_ROOT}/sifs/colabfold_1.6.2-cuda12.sif"
INPUT_FASTA="${WORK_DIR}/data/raw/GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa"
OUTPUT_DIR="${WORK_DIR}/data/processed/colabfold_20260831"

mkdir -p "${OUTPUT_DIR}"

if [[ ! -f "${SIF}" ]]; then
    echo "Error: Singularity image not found: ${SIF}" >&2
    exit 1
fi

if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "Error: Input FASTA not found: ${INPUT_FASTA}" >&2
    exit 1
fi

echo "Starting ColabFold at $(date --iso-8601=seconds)"
echo "Input: ${INPUT_FASTA}"
echo "Output: ${OUTPUT_DIR}"

singularity run --nv \
    -B "${WORK_DIR}:/work" \
    "${SIF}" \
   colabfold_search --mmseqs /share/home/lll/.pixi/bin/mmseqs /work/data/raw/GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa. /data/ssd/colabfold_db msas --gpu 1 
singularity run -B /local/path/to/cache:/cache \
  colabfold_1.6.0-cuda12.sif \
  python -m colabfold.download
echo "ColabFold completed at $(date --iso-8601=seconds)"
