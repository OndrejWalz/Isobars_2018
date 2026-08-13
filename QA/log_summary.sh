#!/bin/bash

# Usage: ./parse_logs.sh <log_dir_or_files>
# Examples:
#   ./parse_logs.sh /path/to/logs/
#   ./parse_logs.sh *.out

OUTPUT="job_summary.csv"

# Collect input files
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <log_file(s) or directory>"
    exit 1
fi

# Build list of log files
LOGFILES=()
for arg in "$@"; do
    if [[ -d "$arg" ]]; then
        while IFS= read -r f; do
            LOGFILES+=("$f")
        done < <(find "$arg" -maxdepth 1 -type f -name "*.out" | sort)
    elif [[ -f "$arg" ]]; then
        LOGFILES+=("$arg")
    else
        echo "Warning: '$arg' not found, skipping."
    fi
done

if [[ ${#LOGFILES[@]} -eq 0 ]]; then
    echo "No log files found."
    exit 1
fi

TOTAL=${#LOGFILES[@]}
CURRENT=0

draw_progress() {
    local current=$1
    local total=$2
    local width=40
    local filled=$(( current * width / total ))
    local empty=$(( width - filled ))
    local bar=""
    for ((i=0; i<filled; i++)); do bar+="█"; done
    for ((i=0; i<empty;  i++)); do bar+="░"; done
    local pct=$(( current * 100 / total ))
    local remaining=$(( total - current ))
    printf "\r  [%s] %3d%%  (%d done, %d remaining)" "$bar" "$pct" "$current" "$remaining"
}

echo "Found $TOTAL log file(s) to process."
echo ""
echo "filename,duration_sec,num_files,num_events" > "$OUTPUT"

for logfile in "${LOGFILES[@]}"; do
    fname=$(basename "$logfile")

    # Extract start and end timestamps
    start_str=$(grep -m1 "We BEGIN on" "$logfile" | sed 's/.*We BEGIN on \[//; s/\]//')
    end_str=$(grep -m1 "We END on"   "$logfile" | sed 's/.*We END on \[//;   s/\]//')

    # Convert to epoch seconds for duration calculation
    if [[ -n "$start_str" && -n "$end_str" ]]; then
        start_epoch=$(date -d "$start_str" +%s 2>/dev/null)
        end_epoch=$(date -d "$end_str"   +%s 2>/dev/null)
        if [[ -n "$start_epoch" && -n "$end_epoch" ]]; then
            duration=$(( end_epoch - start_epoch ))
        else
            duration="N/A"
        fi
    else
        duration="N/A"
    fi

    # Number of input files (from INPUTFILECOUNT line)
    num_files=$(grep -m1 "INPUTFILECOUNT=" "$logfile" | cut -d= -f2 | tr -d '[:space:]')
    [[ -z "$num_files" ]] && num_files="N/A"

    # Number of analyzed events (from "total number of events" line)
    num_events=$(grep "total number of events" "$logfile" | tail -1 | awk '{print $NF}')
    [[ -z "$num_events" ]] && num_events="N/A"

    echo "$fname,$duration,$num_files,$num_events" >> "$OUTPUT"

    CURRENT=$(( CURRENT + 1 ))
    draw_progress "$CURRENT" "$TOTAL"
done

# Print newline after progress bar, then summary
echo ""
echo ""
echo "Done! Results saved to $OUTPUT"

# Top 20 longest jobs (skip header, sort by duration_sec descending, ignore N/A)
echo ""
echo "=========================================="
echo "  Top 20 jobs by duration"
echo "=========================================="
printf "  %-40s  %12s  %10s  %12s\n" "Filename" "Duration (s)" "Inp. Files" "Events"
echo "  ------------------------------------------  ------------  ----------  ------------"
tail -n +2 "$OUTPUT" \
    | grep -v ",N/A," \
    | sort -t',' -k2 -rn \
    | head -20 \
    | while IFS=',' read -r fname dur nfiles nevents; do
        printf "  %-40s  %12s  %10s  %12s\n" "$fname" "$dur" "$nfiles" "$nevents"
    done
echo "=========================================="