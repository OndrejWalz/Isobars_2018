#!/bin/bash
# usage (from the project root):  ./submit/submit_tower.sh [filelist]   (default: see fileList below)
# The list must contain xrootd URLs (root://...), one picoDst per line.

# >>> EDIT THIS: default file list (relative to the project root); an argument overrides it <<<
fileList="${1:-filelists/isobars_2018_1.list}"
baseFolder="$PWD"
case "$fileList" in /*) listOfFiles="$fileList" ;; *) listOfFiles="$PWD/${fileList}" ;; esac   # absolute or project-relative
if [[ ! -f "$listOfFiles" ]]; then
    echo "File list not found (or not a regular file): $listOfFiles"
    exit 1
fi
filelist_name=$(basename "${listOfFiles}" .list)
# macro lives in StRoot/macros/, which is already shipped in the sandbox
macroPath="StRoot/macros/RunPicoTowerTest.C"
rootMacro=$(basename "${macroPath}")
productionId=$(date +%F)
# MUST match the library the macro checks for (SL_version in RunPicoTowerTest.C is SL20c for 2018)
starVersion="SL20c"
maxNFiles=20
year=2018

templateXml="template_tower.xml"
jobFolder="${baseFolder}/submit/${productionId}/job_${filelist_name}"
mkdir -p "${jobFolder}"
cd "${jobFolder}"
mkdir -p report err log list csh production

check=(
    ".sl73_x8664_gcc485"
    "StRoot"
    "${macroPath}"
    "submit/${templateXml}"
)

printf "Checking project … "
for item in "${check[@]}"; do
    path="${baseFolder}/${item}"
    [[ -e $path ]] || { echo "$item missing"; exit 1; }
    [[ -d $path ]] && ln -sf "$path"
done
echo "ok"

generatedXml="generated_${filelist_name}.xml"

cat <<XML >"${generatedXml}"
<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE note [
<!ENTITY rootMacro "${rootMacro}">
<!ENTITY baseFolder "${baseFolder}">
<!ENTITY jobFolder "${jobFolder}">
<!ENTITY listOfFiles "${listOfFiles}">
<!ENTITY starVersion "${starVersion}">
<!ENTITY maxNFiles "${maxNFiles}">
<!ENTITY filelist_name "${filelist_name}">
<!ENTITY year "${year}">
]>

XML
tail -n +2 "${baseFolder}/submit/${templateXml}" >>"${generatedXml}"

star-submit-beta -u ie "${generatedXml}"