#!/bin/bash

# Stop on errors
set -e

# Directory where *this* script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default file names (can be overridden via arguments)
SYSTEM="${1:-Ru}"

if [[ "$SYSTEM" == "Ru" ]]; then
 IN_BASENAME="${2:-data_merged_Ru.root}"
 OUT_BASENAME="${3:-hists_data_Ru.root}"
else
    IN_BASENAME="${2:-data_merged_Zr.root}"
    OUT_BASENAME="${3:-hists_data_Zr.root}"
fi


# Paths relative to the script location
INPUT="${SCRIPT_DIR}/../../trees/${IN_BASENAME}"
OUTPUT="${SCRIPT_DIR}/${OUT_BASENAME}"
MACRO="${SCRIPT_DIR}/make_hists_data.C"

echo "----------------------------------------"
echo "Running DATA histogram production"
echo "Script dir : $SCRIPT_DIR"
echo "Input      : $INPUT"
echo "Output     : $OUTPUT"
echo "Macro      : $MACRO"
echo "----------------------------------------"

root -l -b -q "$MACRO(\"$SYSTEM\",\"$INPUT\",\"$OUTPUT\")"

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"
