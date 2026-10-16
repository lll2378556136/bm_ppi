#!/usr/bin/env bash
#SBATCH --job-name=dpam_bm_ppi
#SBATCH --output=/share/home/lll/projects/bm_ppi/logs/dpam_bm_ppi_%j.out
#SBATCH --error=/share/home/lll/projects/bm_ppi/logs/dpam_bm_ppi_%j.err
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00

# Run one validated DPAM v2.0 dataset inside dpam.sif. The DaliLite overlay
# fixes unreadable permissions in the stock image and is required for step8.
# The step8 patch script is bind-mounted over the in-image copy to fix a
# deadlock where a dali.pl failure killed pool workers without resolving
# their ApplyResults (parent blocked forever on get()).

set -euo pipefail

PROJECT_ROOT="/share/home/lll/projects/bm_ppi"
DATASET="${DATASET:-bm_ppi}"
INPUT_DIR="${INPUT_DIR:-${PROJECT_ROOT}/data/interim/DPAM_input}"
DATABASES_DIR="${DATABASES_DIR:-${PROJECT_ROOT}/data/external/DPAM/databases}"
SIF="${SIF:-${PROJECT_ROOT}/sifs/dpam.sif}"
DALILITE_OVERLAY="${DALILITE_OVERLAY:-${PROJECT_ROOT}/data/external/DPAM/DaliLite.v5}"
DPAM_PATCH_DIR="${DPAM_PATCH_DIR:-${PROJECT_ROOT}/scripts/dpam_patches}"
THREADS="${SLURM_CPUS_PER_TASK:-32}"
# Cap the threads handed to DPAM. hhblits/hhsearch memory grows ~linearly
# with thread count; on gpu-cpu the scheduler raises CPUs to match --mem,
# which previously drove hhblits to 64/80 threads and OOM. 32 threads was
# validated by a probe job (peak RSS 62.5G for the worst protein).
THREADS_CAP="${THREADS_CAP:-32}"
if (( THREADS > THREADS_CAP )); then
    echo "Capping DPAM threads ${THREADS} -> ${THREADS_CAP} (THREADS_CAP)."
    THREADS="${THREADS_CAP}"
fi
MAX_TARGETS_PER_JOB="${MAX_TARGETS_PER_JOB:-500}"
ALLOW_LARGE_DATASET="${ALLOW_LARGE_DATASET:-0}"
RESUME="${RESUME:-0}"

DATASET_DIR="${INPUT_DIR}/${DATASET}"
STRUC_LIST="${INPUT_DIR}/${DATASET}_struc.list"
DOMAINS_OUT="${INPUT_DIR}/${DATASET}_domains"
JOB_TAG="${SLURM_JOB_ID:-manual}_$(date +%Y%m%dT%H%M%S)"
CONTAINER_LOG="${INPUT_DIR}/${DATASET}_container_${JOB_TAG}.log"
VALIDATION_JSON="${INPUT_DIR}/validation.json"

SING_BIN="$(command -v singularity || command -v apptainer || true)"
[[ -n "${SING_BIN}" ]] || {
    echo "Error: neither singularity nor apptainer is available." >&2
    exit 1
}
command -v pixi >/dev/null 2>&1 || {
    echo "Error: pixi is not available for post-run validation." >&2
    exit 1
}
[[ -f "${SIF}" ]] || {
    echo "Error: DPAM image not found: ${SIF}" >&2
    exit 1
}
[[ -f "${DATABASES_DIR}/all_files" ]] || {
    echo "Error: DPAM database manifest not found: ${DATABASES_DIR}/all_files" >&2
    exit 1
}
[[ -x "${DALILITE_OVERLAY}/bin/dali.pl" ]] || {
    echo "Error: executable DaliLite overlay not found: ${DALILITE_OVERLAY}" >&2
    exit 1
}
[[ -f "${DPAM_PATCH_DIR}/step8_iterative_dali.py" ]] || {
    echo "Error: patched step8 script not found: ${DPAM_PATCH_DIR}/step8_iterative_dali.py" >&2
    exit 1
}
[[ -d "${DATASET_DIR}" ]] || {
    echo "Error: dataset directory not found: ${DATASET_DIR}" >&2
    exit 1
}
[[ -s "${STRUC_LIST}" ]] || {
    echo "Error: structure list missing or empty: ${STRUC_LIST}" >&2
    exit 1
}

N_TARGETS=$(grep -cve '^[[:space:]]*$' "${STRUC_LIST}")
N_PDB=$(find "${DATASET_DIR}" -maxdepth 1 -type f \( -name '*.pdb' -o -name '*.cif' \) | wc -l)
N_JSON=$(find "${DATASET_DIR}" -maxdepth 1 -type f -name '*.json' | wc -l)

if (( N_TARGETS > MAX_TARGETS_PER_JOB )) && [[ "${ALLOW_LARGE_DATASET}" != "1" ]]; then
    echo "Error: ${N_TARGETS} targets exceed MAX_TARGETS_PER_JOB=${MAX_TARGETS_PER_JOB}." >&2
    echo "Use the chunked array workflow; set ALLOW_LARGE_DATASET=1 only intentionally." >&2
    exit 1
fi
if (( N_PDB != N_TARGETS || N_JSON != N_TARGETS )); then
    echo "Error: input counts disagree: targets=${N_TARGETS}, PDB=${N_PDB}, JSON=${N_JSON}." >&2
    exit 1
fi
if [[ -s "${DOMAINS_OUT}" && "${RESUME}" != "1" ]]; then
    echo "Error: result already exists: ${DOMAINS_OUT}. Set RESUME=1 to validate/resume." >&2
    exit 1
fi

mkdir -p "${PROJECT_ROOT}/logs"
echo "Starting DPAM at $(date --iso-8601=seconds)"
echo "Job: ${SLURM_JOB_ID:-manual}; node: ${SLURM_JOB_NODELIST:-unknown}"
echo "Dataset: ${DATASET}; targets: ${N_TARGETS}"
echo "Input root: ${INPUT_DIR}"
echo "Threads: ${THREADS}; image: ${SIF}"
echo "Step8 patch: ${DPAM_PATCH_DIR}/step8_iterative_dali.py"

BIND_ARGS=(
    --bind "${DATABASES_DIR}:/mnt/databases:ro"
    --bind "${INPUT_DIR}:/home/${DATASET}:rw"
    --bind "${DALILITE_OVERLAY}:/opt/DaliLite.v5:ro"
    --bind "${DPAM_PATCH_DIR}/step8_iterative_dali.py:/opt/DPAM/scripts/step8_iterative_dali.py:ro"
)

set +e
"${SING_BIN}" exec --fakeroot "${BIND_ARGS[@]}" "${SIF}" /bin/bash -c "
    set -e
    cd /home/${DATASET}
    run_dpam.py ${DATASET} ${THREADS}
" >"${CONTAINER_LOG}" 2>&1
CONTAINER_RC=$?
set -e

gzip -f "${CONTAINER_LOG}"
CONTAINER_LOG_GZ="${CONTAINER_LOG}.gz"

if (( CONTAINER_RC != 0 )); then
    echo "Error: DPAM container exited with code ${CONTAINER_RC}." >&2
    echo "Last 50 container-log lines:" >&2
    gzip -cd "${CONTAINER_LOG_GZ}" | tail -n 50 >&2
    exit "${CONTAINER_RC}"
fi

cd "${PROJECT_ROOT}"
set +e
pixi run python -m bm_ppi.merge_dpam_chunks validate-chunk \
    --input-dir "${INPUT_DIR}" \
    --dataset "${DATASET}" \
    --validation-path "${VALIDATION_JSON}"
VALIDATION_RC=$?
set -e
if (( VALIDATION_RC != 0 )); then
    echo "Error: DPAM post-run validation failed." >&2
    echo "Container log: ${CONTAINER_LOG_GZ}" >&2
    exit "${VALIDATION_RC}"
fi

N_DOMAIN_ROWS=$(awk 'NR > 1 && NF {count++} END {print count + 0}' "${DOMAINS_OUT}")
N_DOMAIN_PROTEINS=$(awk -F'\t' 'NR > 1 && NF {seen[$1]=1} END {print length(seen)}' "${DOMAINS_OUT}")
echo "DPAM finished at $(date --iso-8601=seconds)"
echo "Domains: ${N_DOMAIN_ROWS} rows covering ${N_DOMAIN_PROTEINS} proteins"
echo "Validation: ${VALIDATION_JSON}"
echo "Compressed container log: ${CONTAINER_LOG_GZ}"
