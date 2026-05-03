#!/bin/bash
# ==============================================================
#  monitor_pipeline.sh — Monitor the progress of a running
#  pipeline launched by run_pipeline.sh.
#
#  Usage:
#    bash monitor_pipeline.sh /path/to/output_dir
#    bash monitor_pipeline.sh              # uses default path
#
#    # Auto-refresh every 60 seconds:
#    watch -n 60 --color bash monitor_pipeline.sh /path/to/output_dir
#
#    # Show per-chunk detail for all chunks:
#    bash monitor_pipeline.sh /path/to/output_dir --verbose
#
# ==============================================================

OUTPUT_DIR="${1:-/data/scratch/pipeline_output}"
VERBOSE=false
for arg in "$@"; do
    case "$arg" in
        --verbose|-v) VERBOSE=true ;;
    esac
done

if [ ! -d "${OUTPUT_DIR}" ]; then
    echo "ERROR: Output directory not found: ${OUTPUT_DIR}"
    exit 1
fi

# ---------- Color support ----------
# Colors work in terminal and with 'watch --color'.
# When piped without color support, tput fails gracefully → no codes.
if [ -t 1 ] || [ "${TERM:-}" != "" ]; then
    BOLD=$(tput bold 2>/dev/null || true)
    GREEN=$(tput setaf 2 2>/dev/null || true)
    YELLOW=$(tput setaf 3 2>/dev/null || true)
    RED=$(tput setaf 1 2>/dev/null || true)
    CYAN=$(tput setaf 6 2>/dev/null || true)
    RESET=$(tput sgr0 2>/dev/null || true)
else
    BOLD="" GREEN="" YELLOW="" RED="" CYAN="" RESET=""
fi

echo "=============================================="
echo "${BOLD} Pipeline Progress Monitor${RESET}"
echo " Output dir: ${OUTPUT_DIR}"
echo " Time: $(date)"
echo "=============================================="
echo ""

# Count total SRCIDs from chunks
CHUNK_DIR="${OUTPUT_DIR}/chunks"
if [ -d "${CHUNK_DIR}" ]; then
    TOTAL=$(cat "${CHUNK_DIR}"/chunk_*.txt 2>/dev/null | wc -l | tr -d ' ')
else
    TOTAL="?"
fi

# Log directory
LOG_DIR="${OUTPUT_DIR}/chunk_logs"

# ---------- Per-chunk + aggregate counts ----------
# Uses a single awk pass per log file to count all patterns at once,
# much faster than multiple grep -c calls on large log files.
N_PROCESSED=0
N_SKIPPED=0
N_ALL_FITS=0      # total successful fits (clean + flagged)
N_FLAGGED=0       # fits with quality warnings (flag >= 5)
N_ERRORS=0        # errors during processing

CHUNK_LINES=""
N_CHUNKS_DONE=0
N_CHUNKS_RUNNING=0
N_CHUNKS_STALLED=0
N_CHUNKS_PENDING=0
NOW=$(date +%s)

# Collect recent error messages across all logs
RECENT_ERRORS_FILE=$(mktemp 2>/dev/null || echo "/tmp/monitor_errors_$$")
: > "${RECENT_ERRORS_FILE}"

if [ -d "${LOG_DIR}" ]; then
    for log_file in "${LOG_DIR}"/chunk_*.log; do
        [ -f "${log_file}" ] || continue
        chunk_name=$(basename "${log_file}" .log)

        # How many SRCIDs assigned to this chunk?
        chunk_srcid_file="${CHUNK_DIR}/${chunk_name}.txt"
        if [ -f "${chunk_srcid_file}" ]; then
            n_in_chunk=$(wc -l < "${chunk_srcid_file}" | tr -d ' ')
        else
            n_in_chunk="?"
        fi

        # Single awk pass: count all key patterns at once
        read -r n_processed n_skipped n_fits n_flagged n_errors <<< $(
            awk '
                /\] Processing SRCID/           { proc++ }
                /Skipping SRCID/                { skip++ }
                /Fit completed successfully/    { fits++ }
                /quality flag=/                 { flagged++ }
                /ERROR.*SRCID.*failed/          { errors++ }
                END { printf "%d %d %d %d %d", proc+0, skip+0, fits+0, flagged+0, errors+0 }
            ' "${log_file}" 2>/dev/null
        )

        # Accumulate totals
        N_PROCESSED=$((N_PROCESSED + n_processed))
        N_SKIPPED=$((N_SKIPPED + n_skipped))
        N_ALL_FITS=$((N_ALL_FITS + n_fits))
        N_FLAGGED=$((N_FLAGGED + n_flagged))
        N_ERRORS=$((N_ERRORS + n_errors))

        # Collect error lines for "recent errors" display
        if [ "${n_errors}" -gt 0 ]; then
            grep "ERROR.*SRCID.*failed" "${log_file}" 2>/dev/null | tail -3 >> "${RECENT_ERRORS_FILE}"
        fi

        # Determine status using file mtime (much faster than lsof)
        # stat -c %Y is Linux; stat -f %m is macOS
        file_mtime=$(stat -c %Y "${log_file}" 2>/dev/null || stat -f %m "${log_file}" 2>/dev/null || echo 0)
        file_age=$(( NOW - file_mtime ))

        if grep -q "Batch processing complete" "${log_file}" 2>/dev/null; then
            status_color="${GREEN}"
            status_raw="DONE"
            N_CHUNKS_DONE=$((N_CHUNKS_DONE + 1))
        elif [ "${file_age}" -lt 300 ]; then
            status_color="${CYAN}"
            status_raw="RUNNING"
            N_CHUNKS_RUNNING=$((N_CHUNKS_RUNNING + 1))
        elif [ "${n_processed}" -eq 0 ] 2>/dev/null; then
            status_color=""
            status_raw="PENDING"
            N_CHUNKS_PENDING=$((N_CHUNKS_PENDING + 1))
        else
            # Log hasn't been written in >5 min but isn't done → probably stuck
            stall_min=$(( file_age / 60 ))
            status_color="${RED}"
            status_raw="STALL ${stall_min}m"
            N_CHUNKS_STALLED=$((N_CHUNKS_STALLED + 1))
        fi

        # Per-chunk progress percentage
        chunk_touched=$((n_processed + n_skipped))
        if [ "${n_in_chunk}" != "?" ] && [ "${n_in_chunk}" -gt 0 ] 2>/dev/null; then
            pct=$(( chunk_touched * 100 / n_in_chunk ))
        else
            pct="?"
        fi

        # Get last SRCID being processed (tac reads from end — fast on large logs)
        last_srcid=$(tac "${log_file}" 2>/dev/null | grep -m1 "Processing SRCID" | sed -n 's/.*SRCID \([0-9]*\).*/\1/p' || true)
        last_srcid=${last_srcid:-—}

        CHUNK_LINES="${CHUNK_LINES}$(printf "   %-12s  ${status_color}%-10s${RESET}  %3s%%  fits:%-5s  err:%-3s  last: %s\n" \
            "${chunk_name}" "[${status_raw}]" "${pct}" "${n_fits}" "${n_errors}" "${last_srcid}")\n"
    done
fi

# Derived counts
N_TOUCHED=$((N_PROCESSED + N_SKIPPED))
N_CLEAN=$((N_ALL_FITS - N_FLAGGED))

# Error rate (among sources actually processed, not skipped)
if [ "${N_PROCESSED}" -gt 0 ]; then
    ERROR_RATE=$(echo "scale=1; 100 * ${N_ERRORS} / ${N_PROCESSED}" | bc 2>/dev/null || echo "?")
    FIT_RATE_PCT=$(echo "scale=1; 100 * ${N_ALL_FITS} / ${N_PROCESSED}" | bc 2>/dev/null || echo "?")
else
    ERROR_RATE="—"
    FIT_RATE_PCT="—"
fi

# ---------- Progress bar ----------
if [ "${TOTAL}" != "?" ] && [ "${TOTAL}" -gt 0 ] 2>/dev/null; then
    PCT_OVERALL=$(( N_TOUCHED * 100 / TOTAL ))
    BAR_WIDTH=40
    FILLED=$(( PCT_OVERALL * BAR_WIDTH / 100 ))
    EMPTY=$(( BAR_WIDTH - FILLED ))
    BAR="${GREEN}$(printf "%${FILLED}s" | tr ' ' '#')${RESET}$(printf "%${EMPTY}s" | tr ' ' '-')"
    echo " ${BOLD}Progress: [${BAR}] ${PCT_OVERALL}%${RESET}  (${N_TOUCHED}/${TOTAL})"
    echo ""
fi

# ---------- Overall summary ----------
echo " ${BOLD}Overall:${RESET}"
echo "   Total SRCIDs:    ${TOTAL}"
echo "   Touched:         ${N_TOUCHED}  (processed: ${N_PROCESSED}, skipped: ${N_SKIPPED})"
echo "   Successful fits: ${GREEN}${N_ALL_FITS}${RESET}  (clean: ${N_CLEAN}, flagged: ${N_FLAGGED})  [${FIT_RATE_PCT}%]"
echo "   Errors:          ${N_ERRORS}  [${ERROR_RATE}%]"
echo ""

# ---------- Chunk summary ----------
N_TOTAL_CHUNKS=$((N_CHUNKS_DONE + N_CHUNKS_RUNNING + N_CHUNKS_STALLED + N_CHUNKS_PENDING))
STALL_MSG=""
if [ "${N_CHUNKS_STALLED}" -gt 0 ]; then
    STALL_MSG=", ${RED}${N_CHUNKS_STALLED} stalled${RESET}"
fi
echo " ${BOLD}Chunks:${RESET} ${N_TOTAL_CHUNKS} total — ${GREEN}${N_CHUNKS_DONE} done${RESET}, ${CYAN}${N_CHUNKS_RUNNING} running${RESET}, ${N_CHUNKS_PENDING} pending${STALL_MSG}"

# ---------- Per-chunk detail ----------
if [ -n "${CHUNK_LINES}" ]; then
    if [ "${VERBOSE}" = true ] || [ "${N_TOTAL_CHUNKS}" -le 20 ]; then
        # Show all chunks
        echo ""
        echo "   CHUNK         STATUS       PROG  FITS   ERR   LAST SRCID"
        echo "   ------------ ------------ ----  -----  ---   ----------------"
        printf "%b" "${CHUNK_LINES}"
    else
        # Too many chunks — show only non-DONE ones
        NON_DONE=$(printf "%b" "${CHUNK_LINES}" | grep -v "DONE" || true)
        if [ -n "${NON_DONE}" ]; then
            echo ""
            echo "   Active chunks (${N_CHUNKS_RUNNING} running; use --verbose for all ${N_TOTAL_CHUNKS}):"
            echo "   CHUNK         STATUS       PROG  FITS   ERR   LAST SRCID"
            echo "   ------------ ------------ ----  -----  ---   ----------------"
            echo "${NON_DONE}"
        fi
    fi
    echo ""
fi

# ---------- Running processes ----------
N_RUNNING=$(pgrep -cf "automated_fits.py" 2>/dev/null || echo "0")
echo " Workers alive: ${N_RUNNING} python processes"

# ---------- Result FITS files ----------
N_RESULT_FILES=$(ls "${OUTPUT_DIR}"/fit_results_chunk_*.fits 2>/dev/null | wc -l | tr -d ' ')
MERGED_EXISTS=false
if [ -f "${OUTPUT_DIR}/fit_results_all.fits" ]; then
    MERGED_EXISTS=true
fi

if [ "${N_RESULT_FILES}" -gt 0 ]; then
    # Count total rows using Python (via env var to avoid shell injection)
    VENV_PYTHON="/home/mcoriat/XMM/5XMM/automated_fits/.venv/bin/python3"
    [ -x "${VENV_PYTHON}" ] || VENV_PYTHON="python3"
    TOTAL_ROWS=$(OUTPUT_DIR_VAR="${OUTPUT_DIR}" "${VENV_PYTHON}" -c "
import os, sys
try:
    from astropy.table import Table
    total = 0
    d = os.environ['OUTPUT_DIR_VAR']
    for f in sorted(os.listdir(d)):
        if f.startswith('fit_results_chunk_') and f.endswith('.fits'):
            try: total += len(Table.read(os.path.join(d, f)))
            except: pass
    print(total)
except: print('?')
" 2>/dev/null || echo "?")
    echo " Result files: ${N_RESULT_FILES} chunks, ${BOLD}${TOTAL_ROWS} rows${RESET} written"
    if [ "${MERGED_EXISTS}" = true ]; then
        echo " Merged file:  ${OUTPUT_DIR}/fit_results_all.fits"
    fi
else
    echo " Result files: none yet"
fi

# ---------- Disk usage ----------
# Use df on the mount point (instant) instead of du -sh (walks entire tree)
DISK_USAGE=$(df -h "${OUTPUT_DIR}" 2>/dev/null | tail -1 | awk '{printf "%s used / %s total (%s)", $3, $2, $5}')
MOUNT_POINT=$(df "${OUTPUT_DIR}" 2>/dev/null | tail -1 | awk '{print $NF}')
echo " Disk:         ${DISK_USAGE} (${MOUNT_POINT})"

# ---------- Timing and ETA ----------
if [ "${N_TOUCHED}" -gt 0 ] && [ "${TOTAL}" != "?" ]; then
    START_FILE="${OUTPUT_DIR}/.pipeline_start_time"
    PIPELINE_START=""
    if [ -f "${START_FILE}" ]; then
        PIPELINE_START=$(tr -d '[:space:]' < "${START_FILE}")
    fi
    if [ -n "${PIPELINE_START}" ]; then
        ELAPSED=$((NOW - PIPELINE_START))

        # Only show rates after at least 5 minutes of data
        if [ "${ELAPSED}" -gt 300 ] 2>/dev/null; then
            ELAPSED_H=$(printf "%.1f" "$(echo "scale=2; ${ELAPSED} / 3600" | bc 2>/dev/null)" 2>/dev/null || echo "?")

            PROC_RATE=$(echo "scale=0; ${N_TOUCHED} * 3600 / ${ELAPSED}" | bc 2>/dev/null || echo "?")
            FIT_RATE=$(printf "%.1f" "$(echo "scale=2; ${N_ALL_FITS} * 3600 / ${ELAPSED}" | bc 2>/dev/null)" 2>/dev/null || echo "?")

            REMAINING=$(echo "${TOTAL} - ${N_TOUCHED}" | bc 2>/dev/null || echo 0)
            ETA_H=$(printf "%.1f" "$(echo "scale=2; ${REMAINING} * ${ELAPSED} / ${N_TOUCHED} / 3600" | bc 2>/dev/null)" 2>/dev/null || echo "?")

            # ETA as wall-clock time
            if command -v bc &>/dev/null; then
                ETA_SECS=$(echo "scale=0; ${REMAINING} * ${ELAPSED} / ${N_TOUCHED}" | bc 2>/dev/null || echo 0)
                ETA_TIME=$(date -d "+${ETA_SECS} seconds" "+%a %H:%M" 2>/dev/null || true)
            else
                ETA_TIME=""
            fi

            EXPECTED_FITS="?"
            if [ "${N_PROCESSED}" -gt 0 ]; then
                EXPECTED_FITS=$(echo "scale=0; ${TOTAL} * ${N_ALL_FITS} / ${N_PROCESSED}" | bc 2>/dev/null || echo "?")
            fi

            echo ""
            echo " ${BOLD}Timing${RESET} (elapsed: ${ELAPSED_H}h):"
            echo "   Processing rate: ~${PROC_RATE} sources/hour"
            echo "   Fit rate:        ~${FIT_RATE} fits/hour"
            echo "   Expected fits:   ~${EXPECTED_FITS} (of ${TOTAL} total)"
            if [ -n "${ETA_TIME}" ]; then
                echo "   Est. remaining:  ~${ETA_H} hours (finish ~${ETA_TIME})"
            else
                echo "   Est. remaining:  ~${ETA_H} hours"
            fi
        else
            ELAPSED_MIN=$(( ELAPSED / 60 ))
            echo ""
            echo " Timing: ${ELAPSED_MIN}min elapsed (rates shown after 5min)"
        fi
    fi
fi

# ---------- Recent errors ----------
if [ -s "${RECENT_ERRORS_FILE}" ]; then
    N_ERROR_LINES=$(wc -l < "${RECENT_ERRORS_FILE}" | tr -d ' ')
    echo ""
    echo " ${BOLD}${RED}Recent errors${RESET} (last few from each chunk, ${N_ERROR_LINES} total):"
    # Show at most 8 most recent error lines
    tail -8 "${RECENT_ERRORS_FILE}" | while IFS= read -r line; do
        # Truncate long lines for readability
        if [ ${#line} -gt 100 ]; then
            echo "   ${line:0:97}..."
        else
            echo "   ${line}"
        fi
    done
fi
rm -f "${RECENT_ERRORS_FILE}" 2>/dev/null

echo ""
echo "=============================================="
