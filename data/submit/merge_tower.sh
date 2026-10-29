#!/bin/bash
# usage: ./submit/merge_tower.sh [submit/2026-10-02/job_isobar_2018] [merged.root]   (no argument: newest job folder)
# Merges only jobs whose log shows a clean finish and whose .root file exists.
HERE=$(cd "$(dirname "$0")" && pwd)     # the submit/ folder, so it works from any directory
JOB=${1%/}
if [ -z "$JOB" ]; then
  JOB=$(ls -dt "$HERE"/*/job_* 2>/dev/null | head -1)   # default: newest job folder
  [ -n "$JOB" ] && echo "No job folder given, using newest: $JOB"
fi
if [ -z "$JOB" ] || [ ! -d "$JOB/log" ]; then
  echo "usage: merge_tower.sh [job_folder] [merged.root]   (job folder must contain log/ and production/)"
  echo "       e.g. merge_tower.sh submit/2026-10-02/job_ru_2018"
  exit 1
fi

#output file location 
#OUT=${2:-$JOB/merged.root}
OUT=${2:-../Tower_QA/Isobars_1.root}
NAME=$(basename "$JOB"); NAME=${NAME#job_}
DONE_MSG="writing data and closing output file"

NJOBS=$(ls "$JOB"/csh/*.csh 2>/dev/null | wc -l)
good=(); bad=()
for log in "$JOB"/log/*.out; do
  [ -e "$log" ] || continue
  idx=$(basename "$log" .out); idx=${idx##*_}
  root="$JOB/production/${NAME}_${idx}.root"
  if grep -q "$DONE_MSG" "$log" && [ -f "$root" ]; then good+=("$root"); else bad+=("$idx"); fi
done
# Input files ROOT could not open are skipped silently by the job -> their events are MISSING.
RETRY="$JOB/retry_files.list"
cat "$JOB"/log/*.out "$JOB"/err/*.err 2>/dev/null | grep "open attempt failed on" \
  | sed 's/.*open attempt failed on //' | sort -u > "$RETRY"
NMISS=$(grep -c . "$RETRY")
if [ "$NMISS" -gt 0 ]; then
  echo "WARNING: $NMISS input file(s) could not be opened -> their events are NOT in the output."
  echo "         Listed in $RETRY. Resubmit them, e.g.:  ./submit/submit_tower.sh ${RETRY#./}  (path relative to project root)"
  echo "         then merge the retry job's output together with this one."
else
  rm -f "$RETRY"
fi

echo "jobs submitted: $NJOBS   finished cleanly: ${#good[@]}   incomplete/failed: ${#bad[@]}"
[ ${#bad[@]} -gt 0 ] && echo "Incomplete job indices (check err/ and log/): ${bad[*]}"
[ ${#good[@]} -eq 0 ] && { echo "nothing to merge"; exit 1; }
[ ${#good[@]} -lt "$NJOBS" ] && echo "WARNING: merging only ${#good[@]} of $NJOBS jobs"
hadd -f "$OUT" "${good[@]}" && echo "merged -> $OUT"