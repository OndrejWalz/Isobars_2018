#!/bin/bash

# ──────────────────────────────────────────────────────────────
# run_divide.sh
# Usage:  ./run_divide.sh
#
# At runtime you'll be asked whether to use the original binning
# (no rebinning) or your custom binning (hardcoded below).
# ──────────────────────────────────────────────────────────────

# Your custom bin edges. Must exactly match existing bin boundaries in
# the source histograms (ROOT can only merge whole bins together, not
# split one).

#CUSTOM_BIN_EDGES="0,2,4,6,8,10,12,14,16,18,20,25,30,35,40,60"

CUSTOM_BIN_EDGES="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,35,45,60"

echo "======================================"
echo "   Ru/Zr Jet Spectra Ratio Macro"
echo "======================================"
echo ""
echo "Select operation:"
echo "  1) Divide Ru / Zr      —  MB data"
echo "  2) Divide Ru / Zr      —  HT2 data"
echo "  3) Divide MB / HT2 —  Ru"
echo "  4) Divide MB / HT2 —  Zr"
echo ""
read -p "Enter choice [1-4]: " choice

echo ""
echo "Select binning:"
echo "  1) Original (no rebinning)"
echo "  2) Custom ($CUSTOM_BIN_EDGES)"
read -p "Enter choice [1-2]: " binChoice

case $binChoice in
    1) BIN_EDGES="" ;;
    2) BIN_EDGES="$CUSTOM_BIN_EDGES" ;;
    *) echo "Invalid choice. Exiting."; exit 1 ;;
esac

case $choice in
    1)
        DATA_TYPE="Ru/Zr MB"
        FILE_NUM="hists_data_Ru_MB.root"
        FILE_DEN="hists_data_Zr_MB.root"
        OUTPUT="ratio_Ru_over_Zr_MB.root"
        LABEL_NUM="Ru"
        LABEL_DEN="Zr"
        FILE_EVENTS_NUM="../../../trees/data_merged_Ru.root"
        FILE_EVENTS_DEN="../../../trees/data_merged_Zr.root"
        ;;
    2)
        DATA_TYPE="Ru/Zr HT2"
        FILE_NUM="hists_data_Ru_HT2.root"
        FILE_DEN="hists_data_Zr_HT2.root"
        OUTPUT="ratio_Ru_over_Zr_HT2.root"
        LABEL_NUM="Ru"
        LABEL_DEN="Zr"
        FILE_EVENTS_NUM="../../../trees/data_merged_Ru.root"
        FILE_EVENTS_DEN="../../../trees/data_merged_Zr.root"
        ;;
    3)
        DATA_TYPE="MB / HT2 Ru"
        FILE_NUM="hists_data_Ru_MB.root"
        FILE_DEN="hists_data_Ru_HT2.root"
        OUTPUT="ratio_MB_over_HT2_Ru.root"
        LABEL_NUM="MB"
        LABEL_DEN="HT2"
        FILE_EVENTS_NUM="../../../trees/data_merged_Ru.root"
        FILE_EVENTS_DEN="../../../trees/data_merged_Ru.root"
        ;;
    4)
        DATA_TYPE="MB / HT2 Zr"
        FILE_NUM="hists_data_Zr_MB.root"
        FILE_DEN="hists_data_Zr_HT2.root"
        OUTPUT="ratio_MB_over_HT2_Zr.root"
        LABEL_NUM="MB"
        LABEL_DEN="HT2"
        FILE_EVENTS_NUM="../../../trees/data_merged_Zr.root"
        FILE_EVENTS_DEN="../../../trees/data_merged_Zr.root"
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

# Derive folder name from output root filename (strip .root)
OUTPUT_DIR="${OUTPUT%.root}"

echo ""
echo "--------------------------------------"
echo " Operation    : $DATA_TYPE"
echo " Numerator    : $FILE_NUM"
echo " Denominator  : $FILE_DEN"
echo " Output file  : $OUTPUT"
echo " PDF folder   : $OUTPUT_DIR/"
echo "--------------------------------------"
echo ""

# Check input files exist
for f in "$FILE_NUM" "$FILE_DEN" "$FILE_EVENTS_NUM" "$FILE_EVENTS_DEN"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done

if [ -n "$BIN_EDGES" ]; then
    echo "Rebinning : custom edges: $BIN_EDGES"
else
    echo "Rebinning : none"
fi
echo ""

# Create output folder for PDFs
mkdir -p "$OUTPUT_DIR"
echo "Created folder: $OUTPUT_DIR/"
echo ""

# Run the ROOT macro.
# NOTE: this passes the macro+call as a single command-line argument
# (root4star -q 'file.C(args)'), NOT a heredoc piping ".L"/call/".q" via
# stdin. On this system's root4star, the heredoc/stdin pattern was
# silently doing nothing at all -- no load, no call, no errors, just a
# clean exit -- which is why the macro's own std::cout output never
# appeared and the script still reported SUCCESS. This argument-based
# form is the pattern confirmed to work here (see run_qa_pdfs.sh).
MACRO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root4star -l -b -q "${MACRO_DIR}/divideJetSpectra.C(\"$FILE_NUM\", \"$FILE_DEN\", \"$OUTPUT\", \"$LABEL_NUM\", \"$LABEL_DEN\", \"$OUTPUT_DIR\", \"$FILE_EVENTS_NUM\", \"$FILE_EVENTS_DEN\", \"$BIN_EDGES\")"

if [ $? -eq 0 ]; then
    echo ""
    echo "SUCCESS: Output written to $OUTPUT"
    echo "         PDFs saved in $OUTPUT_DIR/"
else
    echo ""
    echo "ERROR: ROOT macro failed. Check output above."
    exit 1
fi