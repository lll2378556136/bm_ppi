#!/usr/bin/env bash
#SBATCH --job-name=copy_cache_to_ssd
#SBATCH --output=logs/copy_cache_to_ssd_%j.out
#SBATCH --error=logs/copy_cache_to_ssd_%j.err
#SBATCH --partition=gpu
#SBATCH --nodelist=node4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00

set -euo pipefail

# Use SLURM_SUBMIT_DIR so the script works even if SLURM copies it to a spool dir.
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-/share/home/lll/projects/bm_ppi}"
SRC="${PROJECT_ROOT}/cache"
DST="/data/ssd/cache"

echo "Starting cache copy at $(date --iso-8601=seconds)"
echo "Source: ${SRC}"
echo "Destination: ${DST}"
echo "Source size: $(du -sh "${SRC}" | cut -f1)"

# Make sure the SSD parent directory exists.
mkdir -p /data/ssd

# Remove any existing destination to avoid mixed/partial cache.
if [[ -d "${DST}" ]]; then
    echo "Destination already exists; removing it before copying." >&2
    rm -rf "${DST}"
fi

# Copy the cache directory to the SSD.
cp -rL "${SRC}" "${DST}"

echo "Destination size: $(du -sh "${DST}" | cut -f1)"
echo "Cache copy completed at $(date --iso-8601=seconds)"
echo "You can now delete ${SRC} after confirming the SSD copy is usable."
