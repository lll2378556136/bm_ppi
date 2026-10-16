#!/usr/bin/env bash
#SBATCH --job-name=prep_dpam_full
#SBATCH --output=/share/home/lll/projects/bm_ppi/logs/prepare_dpam_full_%j.out
#SBATCH --error=/share/home/lll/projects/bm_ppi/logs/prepare_dpam_full_%j.err
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00

# Prepare canonical full-proteome DPAM inputs and deterministic chunks only.
# This job does not run DPAM or submit any downstream jobs.

set -euo pipefail

PROJECT_ROOT="/share/home/lll/projects/bm_ppi"
cd "${PROJECT_ROOT}"

echo "Starting canonical DPAM input preparation at $(date --iso-8601=seconds)"
pixi run python -m bm_ppi.prepare_dpam_input

echo "Starting DPAM chunk preparation at $(date --iso-8601=seconds)"
pixi run python -m bm_ppi.prepare_dpam_chunks --n-chunks 64

echo "Preparation gate finished at $(date --iso-8601=seconds)"
wc -l data/interim/DPAM_input/bm_ppi_struc.list
du -sh data/interim/DPAM_input/bm_ppi data/interim/DPAM_runs/full_v1
