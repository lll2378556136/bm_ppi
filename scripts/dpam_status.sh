#!/usr/bin/env bash
# Live status panel for DPAM chunk runs (bm_ppi project).
#
# Usage:
#   bash scripts/dpam_status.sh                     # chunks running or partially done
#   bash scripts/dpam_status.sh bm_ppi_c0001 ...    # named datasets
#   watch -n 60 'bash scripts/dpam_status.sh'
set -u

PROJECT_ROOT="/share/home/lll/projects/bm_ppi"
RUNS_DIR="${PROJECT_ROOT}/data/interim/DPAM_runs/full_v1"
CHUNKS_TSV="${RUNS_DIR}/chunks.tsv"
STALE_S=300   # container logs silent longer than this are flagged

now_s=$(date +%s)
echo "dpam-status  $(date '+%F %T')"

# Collect queued dpam array tasks: array index -> "jobid state elapsed node".
declare -A JOBS=()
while read -r jobid state elapsed node; do
    [[ "$jobid" == *_* ]] || continue
    JOBS["${jobid##*_}"]="$jobid $state $elapsed $node"
done < <(squeue -u "${USER:-$(whoami)}" -h -o "%i %T %.10M %N" 2>/dev/null)

# Manifest summary (total / hard-success / queued).
total=0; complete=0
while IFS=$'\t' read -r idx cid ds dir slist nt aa ml wt; do
    ((total += 1))
    if [[ -f "$dir/validation.json" ]] && grep -q '"hard_success": *true' "$dir/validation.json"; then
        ((complete += 1))
    fi
done < <(awk -F'\t' 'NR > 1' "$CHUNKS_TSV")
echo "manifest : $total chunks | $complete hard_success | ${#JOBS[@]} job(s) in queue"

# Rows to display: named datasets, else "interesting" chunks (queued job or
# partial artifacts without hard_success).
rows=()
if (( $# > 0 )); then
    for ds in "$@"; do
        row=$(awk -F'\t' -v ds="$ds" 'NR > 1 && $3 == ds {print; exit}' "$CHUNKS_TSV")
        if [[ -n "$row" ]]; then rows+=("$row"); else echo "  unknown dataset: $ds"; fi
    done
else
    while IFS=$'\t' read -r idx cid ds dir slist nt aa ml wt; do
        if [[ -n "${JOBS[$idx]:-}" ]]; then
            rows+=("$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' "$idx" "$cid" "$ds" "$dir" "$slist" "$nt" "$aa" "$ml" "$wt")")
        elif [[ -e "$dir/${ds}_step1.log" || -e "$dir/step1" ]]; then
            [[ -f "$dir/validation.json" ]] && grep -q '"hard_success": *true' "$dir/validation.json" && continue
            rows+=("$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' "$idx" "$cid" "$ds" "$dir" "$slist" "$nt" "$aa" "$ml" "$wt")")
        fi
    done < <(awk -F'\t' 'NR > 1' "$CHUNKS_TSV")
fi

for row in "${rows[@]}"; do
    IFS=$'\t' read -r idx cid ds dir slist nt aa ml wt <<< "$row"
    echo
    echo "== $ds =="
    echo "  job    : ${JOBS[$idx]:-(none in queue)}"

    # Step logs: first missing or non-'done' step is "current".
    done_steps=(); current=""; flagged=""
    for s in $(seq 1 24); do
        logf="$dir/${ds}_step${s}.log"
        if [[ -f "$logf" ]] && [[ "$(cat "$logf")" == "done" ]]; then
            done_steps+=("$s")
        elif [[ -z "$current" ]]; then
            current="$s"
            [[ -f "$logf" ]] && flagged="(log exists, not done -> will rerun)"
        fi
    done
    [[ ${#done_steps[@]} -gt 0 ]] && echo "  steps  : done -> ${done_steps[*]}"
    echo "  current: step${current:-done?} $flagged"

    # Per-protein progress inside the current step.
    if [[ -n "$current" && -d "$dir/step${current}/$ds" ]]; then
        if [[ "$current" == "8" ]]; then
            h=$(find "$dir/step${current}/$ds" -name '*_hits' | wc -l)
            dn=$(find "$dir/step${current}/$ds" -name '*.done' | wc -l)
            pc=$(( nt > 0 ? 100 * (h + dn) / nt : 0 ))
            echo "  step8  : $h hits + $dn done = $((h + dn))/$nt proteins ($pc%)"
        else
            np=$(ls "$dir/step${current}/$ds" 2>/dev/null | sed 's/\..*//' | sort -u | wc -l)
            echo "  step${current} : ~$np/$nt proteins"
        fi
    fi

    # Newest container log: size + seconds since last write (stall detector).
    log=$(ls -t "$dir/${ds}"_container_*.log 2>/dev/null | head -1)
    if [[ -n "$log" ]]; then
        sz=$(stat -c %s "$log"); mt=$(stat -c %Y "$log"); age=$((now_s - mt))
        flag=""; (( age > STALE_S )) && flag="  <-- STALE?"
        printf "  log    : %s, last write %d s ago%s\n" "$(numfmt --to=iec "$sz" 2>/dev/null || echo "${sz}B")" "$age" "$flag"
    fi

    [[ -f "$dir/${ds}_domains" ]] && echo "  domains: $(wc -l < "$dir/${ds}_domains") rows"
    [[ -f "$dir/validation.json" ]] && echo "  valid  : $(grep -o '"hard_success": *[a-z]*' "$dir/validation.json" | head -1)"
done

exit 0
