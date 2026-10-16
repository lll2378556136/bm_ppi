#!/usr/bin/env bash
# Launcher for four parallel ColabFold batch jobs across node3 and node4.

set -euo pipefail

# Resolve the absolute project root from this script's location, independent of the caller's cwd.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
LOG_DIR="${PROJECT_ROOT}/logs"
WORKER="${PROJECT_ROOT}/scripts/colabfold_batch.slurm"

# Ensure log directory exists before submitting.
mkdir -p "${LOG_DIR}"

if [[ ! -f "${WORKER}" ]]; then
    echo "Error: worker script not found: ${WORKER}" >&2
    exit 1
fi

# Export the absolute project root so the worker is not tied to the submission cwd.
export PROJECT_ROOT

# Mapping: 001/002 -> node3, 003/004 -> node4 (1 GPU each).
for BID in 001 002 003 004; do
    case "${BID}" in
        001|002)
            NODE="node3"
            ;;
        003|004)
            NODE="node4"
            ;;
        *)
            echo "Error: unexpected BATCH_ID ${BID}" >&2
            exit 1
            ;;
    esac

    OUT_LOG="${LOG_DIR}/colabfold_batch_${BID}_%j.out"
    ERR_LOG="${LOG_DIR}/colabfold_batch_${BID}_%j.err"

    # Submit one single-GPU worker per batch with a unique per-job log.
    sbatch \
        --job-name="colabfold_batch_${BID}" \
        --nodelist="${NODE}" \
        --output="${OUT_LOG}" \
        --error="${ERR_LOG}" \
        --export=ALL,BATCH_ID="${BID}" \
        "${WORKER}"
done
