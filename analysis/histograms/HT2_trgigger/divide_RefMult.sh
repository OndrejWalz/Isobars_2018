#!/bin/bash
#
# run_divide.sh
# Runs the divide_hrefmult.C ROOT macro on two input .root files.
#
# Usage:
#   ./run_divide.sh file1.root file2.root [output_base_name]
#
# Example:
#   ./run_divide.sh data_merged_Ru.root data_merged_Zr.root ratio_Ru_Zr

set -e




MACRO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO="${MACRO_DIR}/divide_hrefmult.C"

FILE1="${MACRO_DIR}/../../../trees/data_merged_Ru_HT2.root"
FILE2="${MACRO_DIR}/../../../trees/data_merged_Zr_HT2.root"
OUTBASE="${3:-hRefMultCorr_ratio}"

if [ ! -f "$FILE1" ]; then
    echo "ERROR: file not found: $FILE1"
    exit 1
fi

if [ ! -f "$FILE2" ]; then
    echo "ERROR: file not found: $FILE2"
    exit 1
fi

if [ ! -f "$MACRO" ]; then
    echo "ERROR: macro not found: $MACRO"
    exit 1
fi

echo "Running divide_hrefmult.C on:"
echo "  File 1: $FILE1"
echo "  File 2: $FILE2"
echo "  Output base name: $OUTBASE"
echo ""

root4star -l -b -q "${MACRO}(\"${FILE1}\", \"${FILE2}\", \"${OUTBASE}\")"