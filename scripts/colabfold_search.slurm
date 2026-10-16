#!/usr/bin/env bash
#SBATCH --job-name=colabfold_search
#SBATCH --output=logs/colabfold_search_%j.out
#SBATCH --error=logs/colabfold_search_%j.err
#SBATCH --partition=gpu
#SBATCH --nodelist=node4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --mem=110G
#SBATCH --gres=gpu:2
#SBATCH --time=72:00:00

set -euo pipefail

# Use SLURM_SUBMIT_DIR so the script works even if SLURM copies it to a spool dir.
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-/share/home/lll/projects/bm_ppi}"
SIF="${PROJECT_ROOT}/sifs/colabfold_1.6.2-cuda12.sif"

# Use a short run ID so that the MSA directory and the later batch output match.
RUN_ID="gcf"

# Input FASTA. Use the pilot for testing; change this to the full GCF file for production.
INPUT_FASTA="${PROJECT_ROOT}/data/raw/GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa" 
# INPUT_FASTA="${PROJECT_ROOT}/data/raw/GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa"

# Local MMseqs2 databases on node4's SSD.
DB_DIR="/data/ssd/colabfold_db"
MSA_DIR="${PROJECT_ROOT}/data/interim/colabfold_msas_${RUN_ID}"

mkdir -p "${PROJECT_ROOT}/logs" "${MSA_DIR}"

if [[ ! -f "${SIF}" ]]; then
    echo "Error: Singularity image not found: ${SIF}" >&2
    exit 1
fi

if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "Error: Input FASTA not found: ${INPUT_FASTA}" >&2
    exit 1
fi

if [[ ! -d "${DB_DIR}" ]]; then
    echo "Error: ColabFold database directory not found: ${DB_DIR}" >&2
    exit 1
fi

echo "Starting ColabFold MSA search at $(date --iso-8601=seconds)"
echo "Input: ${INPUT_FASTA}"
echo "DB: ${DB_DIR}"
echo "Output MSAs: ${MSA_DIR}"

# Strip the project root prefix so that paths are relative to the /work mount.
INPUT_IN_WORK="/work/${INPUT_FASTA#"${PROJECT_ROOT}/"}"
MSA_IN_WORK="/work/${MSA_DIR#"${PROJECT_ROOT}/"}"

# Mount project (work) and database dir (read-only) into the container.
# The ghcr.io/sokrypton/colabfold image has mmseqs2 installed, so we use the container's binary.
singularity run --nv \
    -B "${PROJECT_ROOT}:/work" \
    -B "${DB_DIR}:/db:ro" \
    "${SIF}" \
    colabfold_search --mmseqs mmseqs \
      "${INPUT_IN_WORK}" \
      /db \
      "${MSA_IN_WORK}" --gpu 1\
      --threads "${SLURM_CPUS_PER_TASK}"

echo "ColabFold MSA search completed at $(date --iso-8601=seconds)"
