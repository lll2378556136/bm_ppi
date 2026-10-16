#!/usr/bin/env bash
#SBATCH --job-name=deeploc_accurate
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --gres=gpu:1
#SBATCH --time=168:00:00

set -euo pipefail

TOTAL_SEQUENCES=13459
BATCH_001_SEQUENCES=6730
BATCH_002_SEQUENCES=6729

# Workers use the project root exported by submit mode because Slurm may spool the script.
if [[ -n "${PROJECT_ROOT:-}" ]]; then
    PROJECT_ROOT=$(realpath -- "${PROJECT_ROOT}")
    SCRIPT_PATH="${PROJECT_ROOT}/scripts/deeploc_accurate.slurm"
else
    SCRIPT_PATH=$(realpath -- "${BASH_SOURCE[0]}")
    PROJECT_ROOT=$(realpath -- "$(dirname -- "${SCRIPT_PATH}")/..")
fi

INPUT_FASTA="${PROJECT_ROOT}/data/raw/GCF_030269925.1_ASM3026992v2_longest_protein_by_gene.faa"
LOCALIZATION_DIR="${PROJECT_ROOT}/data/interim/localization"
BATCH_DIR="${LOCALIZATION_DIR}/batches"
HF_HOME="${PROJECT_ROOT}/scratch/hf_cache"
FINAL_CSV="${LOCALIZATION_DIR}/deeploc_2.1_accurate_prott5.csv"
LOG_DIR="${PROJECT_ROOT}/logs"
MODE="${DEEPLOC_MODE:-${1:-}}"

submit_jobs() {
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "Error: Submit mode must run outside a Slurm job." >&2
        exit 1
    fi

    mkdir -p "${LOG_DIR}" "${LOCALIZATION_DIR}"

    if [[ ! -s "${INPUT_FASTA}" ]]; then
        echo "Error: Input FASTA is missing or empty: ${INPUT_FASTA}" >&2
        exit 1
    fi

    local sequence_count
    sequence_count=$(awk '/^>/{count++} END{print count+0}' "${INPUT_FASTA}")
    if [[ "${sequence_count}" -ne "${TOTAL_SEQUENCES}" ]]; then
        echo "Error: Expected ${TOTAL_SEQUENCES} sequences, found ${sequence_count}." >&2
        exit 1
    fi

    local run_id batch_001_job batch_002_job merge_job
    run_id="$(date +%Y%m%d-%H%M%S)_$$"

    batch_001_job=$(sbatch --parsable \
        --job-name="deeploc_001" \
        --nodelist=node3 \
        --chdir="${PROJECT_ROOT}" \
        --output="${LOG_DIR}/deeploc_accurate_001_%j.out" \
        --error="${LOG_DIR}/deeploc_accurate_001_%j.err" \
        --export="ALL,PROJECT_ROOT=${PROJECT_ROOT},DEEPLOC_MODE=worker,BATCH_ID=001,RUN_ID=${run_id}" \
        "${SCRIPT_PATH}")
    batch_001_job="${batch_001_job%%;*}"

    batch_002_job=$(sbatch --parsable \
        --job-name="deeploc_002" \
        --nodelist=node4 \
        --chdir="${PROJECT_ROOT}" \
        --output="${LOG_DIR}/deeploc_accurate_002_%j.out" \
        --error="${LOG_DIR}/deeploc_accurate_002_%j.err" \
        --export="ALL,PROJECT_ROOT=${PROJECT_ROOT},DEEPLOC_MODE=worker,BATCH_ID=002,RUN_ID=${run_id}" \
        "${SCRIPT_PATH}")
    batch_002_job="${batch_002_job%%;*}"

    merge_job=$(sbatch --parsable \
        --job-name="deeploc_merge" \
        --dependency="afterok:${batch_001_job}:${batch_002_job}" \
        --gres=none \
        --cpus-per-task=1 \
        --mem=4G \
        --time=00:30:00 \
        --chdir="${PROJECT_ROOT}" \
        --output="${LOG_DIR}/deeploc_merge_%j.out" \
        --error="${LOG_DIR}/deeploc_merge_%j.err" \
        --export="ALL,PROJECT_ROOT=${PROJECT_ROOT},DEEPLOC_MODE=merge,RUN_ID=${run_id}" \
        "${SCRIPT_PATH}")
    merge_job="${merge_job%%;*}"

    echo "Submitted DeepLoc Accurate jobs for run ${run_id}:"
    echo "  batch 001: job ${batch_001_job}, node3, ${BATCH_001_SEQUENCES} sequences"
    echo "  batch 002: job ${batch_002_job}, node4, ${BATCH_002_SEQUENCES} sequences"
    echo "  merge:     job ${merge_job}, afterok:${batch_001_job}:${batch_002_job}"
    echo "  final CSV: ${FINAL_CSV}"
}

run_worker() {
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
        echo "Error: Worker mode must run inside a Slurm job." >&2
        exit 1
    fi
    if [[ -z "${RUN_ID:-}" ]]; then
        echo "Error: RUN_ID is required in worker mode." >&2
        exit 1
    fi

    local expected_nodes start_record end_record expected_sequences
    case "${BATCH_ID:-}" in
        001)
            expected_nodes="node3 node4"
            start_record=1
            end_record="${BATCH_001_SEQUENCES}"
            expected_sequences="${BATCH_001_SEQUENCES}"
            ;;
        002)
            expected_nodes="node3 node4"
            start_record=$((BATCH_001_SEQUENCES + 1))
            end_record="${TOTAL_SEQUENCES}"
            expected_sequences="${BATCH_002_SEQUENCES}"
            ;;
        *)
            echo "Error: BATCH_ID must be 001 or 002, got '${BATCH_ID:-}'." >&2
            exit 1
            ;;
    esac

    if [[ ! " ${expected_nodes} " =~ " ${SLURMD_NODENAME:-unknown} " ]]; then
        echo "Error: Batch ${BATCH_ID} must run on one of [${expected_nodes}], got '${SLURMD_NODENAME:-unknown}'." >&2
        exit 1
    fi
    if [[ -z "${CUDA_VISIBLE_DEVICES:-}" ]]; then
        echo "Error: No GPU was assigned to batch ${BATCH_ID}." >&2
        exit 1
    fi
    if [[ ! -s "${INPUT_FASTA}" ]]; then
        echo "Error: Input FASTA is missing or empty: ${INPUT_FASTA}" >&2
        exit 1
    fi

    local pixi_bin prott5_cache batch_fasta batch_output split_tmp result_marker marker_tmp
    pixi_bin="${PIXI_BIN:-$(command -v pixi || true)}"
    prott5_cache="${HF_HOME}/hub/models--Rostlab--prot_t5_xl_uniref50"
    batch_fasta="${BATCH_DIR}/deeploc_accurate_batch${BATCH_ID}.faa"
    batch_output="${LOCALIZATION_DIR}/deeploc_accurate_batch${BATCH_ID}"
    result_marker="${batch_output}/result_${RUN_ID}.path"
    split_tmp="${batch_fasta}.tmp.${SLURM_JOB_ID}"
    marker_tmp="${result_marker}.tmp.${SLURM_JOB_ID}"

    if [[ -z "${pixi_bin}" || ! -x "${pixi_bin}" ]]; then
        echo "Error: pixi executable was not found." >&2
        exit 1
    fi
    if [[ ! -x "${PROJECT_ROOT}/.pixi/envs/deeploc/bin/deeploc2" ]]; then
        echo "Error: DeepLoc executable is missing from the pixi deeploc environment." >&2
        exit 1
    fi
    if [[ ! -d "${prott5_cache}/snapshots" ]]; then
        echo "Error: Cached ProtT5 model was not found: ${prott5_cache}" >&2
        exit 1
    fi

    mkdir -p "${BATCH_DIR}" "${batch_output}"
    trap 'rm -f -- "${split_tmp:-}" "${marker_tmp:-}"' EXIT

    awk -v start="${start_record}" -v end="${end_record}" '
        /^>/ {record++}
        record >= start && record <= end {print}
    ' "${INPUT_FASTA}" > "${split_tmp}"

    local sequence_count long_sequence_count
    sequence_count=$(awk '/^>/{count++} END{print count+0}' "${split_tmp}")
    if [[ "${sequence_count}" -ne "${expected_sequences}" ]]; then
        echo "Error: Batch ${BATCH_ID} contains ${sequence_count} sequences; expected ${expected_sequences}." >&2
        exit 1
    fi
    mv -- "${split_tmp}" "${batch_fasta}"

    long_sequence_count=$(awk '
        /^>/ {
            if (seen && seq_len > 4000) count++
            seen = 1
            seq_len = 0
            next
        }
        {
            gsub(/[[:space:]]/, "")
            seq_len += length($0)
        }
        END {
            if (seen && seq_len > 4000) count++
            print count+0
        }
    ' "${batch_fasta}")

    export HF_HOME
    export HF_HUB_OFFLINE=1
    export TRANSFORMERS_OFFLINE=1
    export TOKENIZERS_PARALLELISM=false
    export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
    export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

    local start_epoch latest_result latest_mtime result result_mtime result_rows
    start_epoch=$(date +%s)

    echo "Starting DeepLoc 2.1 Accurate batch ${BATCH_ID} at $(date --iso-8601=seconds)"
    echo "Run ID: ${RUN_ID}"
    echo "Node: ${SLURMD_NODENAME}"
    echo "GPU: ${CUDA_VISIBLE_DEVICES}"
    echo "Input: ${batch_fasta}"
    echo "Sequences: ${sequence_count}"
    echo "Output: ${batch_output}"
    echo "Notice: ${long_sequence_count} sequences exceed 4000 aa and will retain their first and last 2000 aa."

    "${pixi_bin}" run --manifest-path "${PROJECT_ROOT}/pixi.toml" -e deeploc \
        deeploc2 \
        --fasta "${batch_fasta}" \
        --output "${batch_output}" \
        --model Accurate \
        --device cuda

    latest_result=""
    latest_mtime=0
    shopt -s nullglob
    for result in "${batch_output}"/results_*.csv; do
        result_mtime=$(stat -c %Y "${result}")
        if (( result_mtime >= start_epoch && result_mtime > latest_mtime )); then
            latest_result="${result}"
            latest_mtime="${result_mtime}"
        fi
    done
    shopt -u nullglob

    if [[ -z "${latest_result}" ]]; then
        echo "Error: DeepLoc did not create a new results CSV for batch ${BATCH_ID}." >&2
        exit 1
    fi

    result_rows=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "${latest_result}")
    if [[ "${result_rows}" -ne "${expected_sequences}" ]]; then
        echo "Error: Batch ${BATCH_ID} result has ${result_rows} rows; expected ${expected_sequences}." >&2
        exit 1
    fi

    printf '%s\n' "${latest_result}" > "${marker_tmp}"
    mv -- "${marker_tmp}" "${result_marker}"
    trap - EXIT

    echo "Result: ${latest_result}"
    echo "DeepLoc batch ${BATCH_ID} completed at $(date --iso-8601=seconds)"
}

run_merge() {
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
        echo "Error: Merge mode must run inside a Slurm job." >&2
        exit 1
    fi
    if [[ -z "${RUN_ID:-}" ]]; then
        echo "Error: RUN_ID is required in merge mode." >&2
        exit 1
    fi

    local output_001 output_002 marker_001 marker_002 csv_001 csv_002
    output_001="${LOCALIZATION_DIR}/deeploc_accurate_batch001"
    output_002="${LOCALIZATION_DIR}/deeploc_accurate_batch002"
    marker_001="${output_001}/result_${RUN_ID}.path"
    marker_002="${output_002}/result_${RUN_ID}.path"

    if [[ ! -s "${marker_001}" || ! -s "${marker_002}" ]]; then
        echo "Error: Missing result marker for run ${RUN_ID}." >&2
        exit 1
    fi
    csv_001=$(<"${marker_001}")
    csv_002=$(<"${marker_002}")
    if [[ ! -s "${csv_001}" || ! -s "${csv_002}" ]]; then
        echo "Error: A batch result CSV is missing or empty." >&2
        exit 1
    fi

    local rows_001 rows_002 header_001 header_002 merged_tmp expected_ids actual_ids
    rows_001=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "${csv_001}")
    rows_002=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "${csv_002}")
    if [[ "${rows_001}" -ne "${BATCH_001_SEQUENCES}" || "${rows_002}" -ne "${BATCH_002_SEQUENCES}" ]]; then
        echo "Error: Batch result sizes are ${rows_001} and ${rows_002}; expected ${BATCH_001_SEQUENCES} and ${BATCH_002_SEQUENCES}." >&2
        exit 1
    fi

    header_001=$(head -n 1 "${csv_001}")
    header_002=$(head -n 1 "${csv_002}")
    if [[ "${header_001}" != "${header_002}" ]]; then
        echo "Error: Batch CSV headers do not match." >&2
        exit 1
    fi

    merged_tmp="${FINAL_CSV}.tmp.${SLURM_JOB_ID}"
    expected_ids="${FINAL_CSV}.expected_ids.tmp.${SLURM_JOB_ID}"
    actual_ids="${FINAL_CSV}.actual_ids.tmp.${SLURM_JOB_ID}"
    trap 'rm -f -- "${merged_tmp:-}" "${expected_ids:-}" "${actual_ids:-}"' EXIT

    {
        printf '%s\n' "${header_001}"
        tail -n +2 "${csv_001}"
        tail -n +2 "${csv_002}"
    } > "${merged_tmp}"

    local merged_rows unique_ids
    merged_rows=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "${merged_tmp}")
    unique_ids=$(awk -F, 'NR > 1 {seen[$1]++} END{print length(seen)}' "${merged_tmp}")
    if [[ "${merged_rows}" -ne "${TOTAL_SEQUENCES}" || "${unique_ids}" -ne "${TOTAL_SEQUENCES}" ]]; then
        echo "Error: Merged CSV has ${merged_rows} rows and ${unique_ids} unique IDs; expected ${TOTAL_SEQUENCES}." >&2
        exit 1
    fi

    awk '/^>/{id=$1; sub(/^>/, "", id); print id}' "${INPUT_FASTA}" > "${expected_ids}"
    awk -F, 'NR > 1 {print $1}' "${merged_tmp}" > "${actual_ids}"
    if ! cmp -s "${expected_ids}" "${actual_ids}"; then
        echo "Error: Merged Protein_ID values or order do not match the source FASTA." >&2
        exit 1
    fi

    mv -- "${merged_tmp}" "${FINAL_CSV}"
    rm -f -- "${expected_ids}" "${actual_ids}"
    trap - EXIT

    echo "Merged DeepLoc results at $(date --iso-8601=seconds)"
    echo "Run ID: ${RUN_ID}"
    echo "Batch 001: ${csv_001} (${rows_001} rows)"
    echo "Batch 002: ${csv_002} (${rows_002} rows)"
    echo "Final CSV: ${FINAL_CSV} (${merged_rows} rows)"
}

case "${MODE}" in
    submit)
        submit_jobs
        ;;
    worker)
        run_worker
        ;;
    merge)
        run_merge
        ;;
    *)
        echo "Usage: bash ${SCRIPT_PATH} submit" >&2
        exit 1
        ;;
esac
