#!/usr/bin/env bash
# SEXmer-dump.sh - Generate a filtered k-mer dump file for a single sample.
# Author: Dede Kurniawan

set -euo pipefail
export LC_ALL=C

# defaults
KMER_SIZE=21
MIN_COUNT=3
MEM="16G"
TRIGGER_SEQ=""
THREADS=4
PREFIX=""
TMPDIR_BASE="$(pwd)"
KMC_BIN=""
READS=()

# log helpers
# All log output goes to stderr to keep stdout clean for potential piping.
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*" >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"  >&2; }

# usage
usage() {
    cat >&2 <<EOF

SEXmer-dump.sh - Generate a filtered k-mer dump file for a single sample.
Requires KMC (kmc + kmc_dump) available on PATH or via --kmc-bin.

Usage: SEXmer-dump.sh --prefix <sample> <reads_1.fq.gz> [reads_2.fq.gz] [OPTIONS]

Mandatory:
  --prefix             Output filename prefix (output: <prefix>.dump.gz)
  <reads>              One or two FASTQ files (.fq, .fastq, .fq.gz, .fastq.gz)

Optional:
  -k, --kmer-size      K-mer size (1-63)                              [default: ${KMER_SIZE}]
  --min-count          Minimum k-mer count (KMC -ci)                  [default: ${MIN_COUNT}]
  --mem                RAM budget for KMC (e.g. 16G, 512M)            [default: ${MEM}]
  --trigger-seq        Retain only k-mers starting with this seq      [default: off]
  -t, --threads        CPU threads                                    [default: ${THREADS}]
  --tmpdir             Parent directory for the temporary work folder [default: current dir]
  --kmc-bin            Directory containing kmc and kmc_dump binaries [default: PATH]
  -h, --help           Show this help and exit

EOF
    exit 1
}

[[ $# -eq 0 ]] && usage

# argument parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)         PREFIX="$2";        shift 2 ;;
        -k|--kmer-size)   KMER_SIZE="$2";    shift 2 ;;
        --min-count)      MIN_COUNT="$2";    shift 2 ;;
        --mem)            MEM="$2";          shift 2 ;;
        --trigger-seq)    TRIGGER_SEQ="$2";  shift 2 ;;
        -t|--threads)     THREADS="$2";      shift 2 ;;
        --tmpdir)         TMPDIR_BASE="$2";  shift 2 ;;
        --kmc-bin)        KMC_BIN="$2";      shift 2 ;;
        -h|--help)        usage ;;
        -*) error "Unknown option '$1'"; usage ;;
        *)  READS+=("$1"); shift ;;
    esac
done

# validate mandatory arguments
[[ -z "$PREFIX" ]] && { error "--prefix is required."; usage; }

[[ ${#READS[@]} -eq 0 ]] && { error "At least one input read file is required."; usage; }

[[ ${#READS[@]} -gt 2 ]] && {
    error "At most two read files (paired-end) are accepted. Got: ${READS[*]}"; exit 1; }

# validate options
[[ "$KMER_SIZE" =~ ^[1-9][0-9]*$ ]] && [[ "$KMER_SIZE" -le 63 ]] || {
    error "--kmer-size must be an integer between 1 and 63."; exit 1; }

[[ "$MIN_COUNT" =~ ^[1-9][0-9]*$ ]] || {
    error "--min-count must be a positive integer."; exit 1; }

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    error "--threads must be a positive integer."; exit 1; }

[[ "$MEM" =~ ^[0-9]+[GgMm]$ ]] || {
    error "--mem must be a number followed by G or M (e.g. 16G, 512M)."; exit 1; }

if [[ -n "$TRIGGER_SEQ" ]]; then
    [[ "$TRIGGER_SEQ" =~ ^[ACGTacgt]+$ ]] || {
        error "--trigger-seq must contain only ACGT characters. Got: '${TRIGGER_SEQ}'"; exit 1; }
    TRIGGER_SEQ="${TRIGGER_SEQ^^}"
fi

# validate --kmc-bin directory if provided
if [[ -n "$KMC_BIN" ]]; then
    [[ -d "$KMC_BIN" ]] || { error "--kmc-bin directory does not exist: ${KMC_BIN}"; exit 1; }
    export PATH="${KMC_BIN}:${PATH}"
fi

# validate directories
[[ -d "$TMPDIR_BASE" ]] || { error "Temporary parent directory does not exist: ${TMPDIR_BASE}"; exit 1; }
[[ -w "$TMPDIR_BASE" ]] || { error "Temporary parent directory is not writable: ${TMPDIR_BASE}"; exit 1; }

# validate output directory derived from prefix
# dirname returns "." for bare prefixes — resolve to the actual current directory.
OUT_DIR="$(dirname "$PREFIX")"
[[ "$OUT_DIR" == "." ]] && OUT_DIR="$(pwd)"
[[ -d "$OUT_DIR" ]] || { error "Output directory does not exist: ${OUT_DIR}"; exit 1; }
[[ -w "$OUT_DIR" ]] || { error "Output directory is not writable: ${OUT_DIR}"; exit 1; }

# validate input files
for f in "${READS[@]}"; do
    [[ -r "$f" ]] || { error "Cannot read input file: $f"; exit 1; }
done

# check KMC binaries
for bin in kmc kmc_dump; do
    command -v "$bin" &>/dev/null || {
        error "'${bin}' not found on PATH. Install KMC: conda install bioconda::kmc"
        exit 1; }
done

# set up temp directory
# TMPDIR_BASE : user-specified parent directory (--tmpdir)
# DUMP_TMPDIR : actual working directory created for this run
DUMP_TMPDIR="${TMPDIR_BASE}/sexmer_dump_tmp_$$"
mkdir -p "$DUMP_TMPDIR"
cleanup() { rm -rf "$DUMP_TMPDIR"; }
trap cleanup EXIT

# mem_to_kmc_gb
# KMC -m flag expects an integer in GB. Convert the user-supplied MEM string
# (e.g. 16G, 512M) to a ceiling GB value with a minimum of 1.
mem_to_kmc_gb() {
    local mem="$1"
    local num unit
    num="${mem//[GgMm]/}"
    unit="${mem//[0-9]/}"
    unit="${unit^^}"
    if [[ "$unit" == "G" ]]; then
        echo "$num"
    else
        local gb=$(( (num + 1023) / 1024 ))
        [[ "$gb" -lt 1 ]] && gb=1
        echo "$gb"
    fi
}
KMC_MEM_GB=$(mem_to_kmc_gb "$MEM")

# derive output path from prefix
OUTPUT="${PREFIX}.dump.gz"

# log run parameters
info "SEXmer-dump starting"
info "Parameters : kmer-size=${KMER_SIZE}, min-count=${MIN_COUNT}, mem=${MEM}, threads=${THREADS}"
info "Trigger-seq: ${TRIGGER_SEQ:-off}"
info "Input reads: ${READS[*]}"
info "Output     : ${OUTPUT}"
info "Temp dir   : ${DUMP_TMPDIR}"

# STEP 1: Build KMC file-of-files
# KMC requires a plain text file listing input FASTQ paths, one per line.
info "Building KMC input file list..."

FILES_LIST="${DUMP_TMPDIR}/input_files.lst"
printf '%s\n' "${READS[@]}" > "$FILES_LIST"

info "  ${#READS[@]} file(s) listed."

# STEP 2: Run KMC
# -k   k-mer size
# -ci  minimum count (discard k-mers below this threshold at counting time)
# -m   RAM in GB
# -t   threads
# Format flag (-fm/-fq) is intentionally omitted: when input is supplied via
# a @file-of-files, KMC auto-detects format from each file's extension
# (.fq, .fastq, .fq.gz, .fastq.gz). Passing -fm in this mode overrides
# auto-detection and causes "Wrong input file" errors for some extensions.
# KMC writes two files: <db_prefix>.kmc_pre and <db_prefix>.kmc_suf
KMC_DB="${DUMP_TMPDIR}/kmc_db"
KMC_TMP="${DUMP_TMPDIR}/kmc_tmp"
mkdir -p "$KMC_TMP"

info "Running KMC (k=${KMER_SIZE}, min-count=${MIN_COUNT}, mem=${KMC_MEM_GB}G, threads=${THREADS})..."

kmc \
    -k"${KMER_SIZE}" \
    -ci"${MIN_COUNT}" \
    -m"${KMC_MEM_GB}" \
    -t"${THREADS}" \
    @"${FILES_LIST}" \
    "${KMC_DB}" \
    "${KMC_TMP}" \
    2>&1 | while IFS= read -r line; do info "  [kmc] $line"; done

[[ -f "${KMC_DB}.kmc_pre" && -f "${KMC_DB}.kmc_suf" ]] || {
    error "KMC did not produce expected database files. Check log above."; exit 1; }

info "KMC database built successfully."

# STEP 3: Dump k-mers to text
# kmc_dump writes space-separated KMER COUNT, one per line, sorted by k-mer.
# -ci is repeated here as a safety net in case KMC -ci behaved unexpectedly.
RAW_DUMP="${DUMP_TMPDIR}/raw_dump.txt"

info "Dumping k-mers from KMC database..."

kmc_dump \
    -ci"${MIN_COUNT}" \
    "${KMC_DB}" \
    "${RAW_DUMP}" \
    2>&1 | while IFS= read -r line; do info "  [kmc_dump] $line"; done

[[ -f "$RAW_DUMP" ]] || {
    error "kmc_dump did not produce output. Check log above."; exit 1; }

RAW_COUNT=$(wc -l < "$RAW_DUMP")
info "  ${RAW_COUNT} k-mers in raw dump."

# STEP 4: Apply trigger-seq filter (optional)
# Streaming awk pass that keeps only k-mers whose sequence starts with the
# required prefix. No intermediate file is created when the filter is off.
if [[ -n "$TRIGGER_SEQ" ]]; then
    info "Applying trigger-seq filter: keeping k-mers starting with '${TRIGGER_SEQ}'..."
    FILTERED_DUMP="${DUMP_TMPDIR}/filtered_dump.txt"
    awk -v prefix="$TRIGGER_SEQ" '$1 ~ ("^" prefix)' "$RAW_DUMP" > "$FILTERED_DUMP"
    FILTERED_COUNT=$(wc -l < "$FILTERED_DUMP")
    DROPPED=$(( RAW_COUNT - FILTERED_COUNT ))
    info "  ${FILTERED_COUNT} k-mers retained, ${DROPPED} removed by trigger-seq filter."
    FINAL_DUMP="$FILTERED_DUMP"
else
    info "No trigger-seq filter applied."
    FINAL_DUMP="$RAW_DUMP"
    FILTERED_COUNT=$RAW_COUNT
fi

# STEP 5: Write gzip-compressed output
# gzip -1 (fastest compression) keeps the file small while minimising wall time.
# The output is compatible with gzip -dc used in SEXmer-scan.sh open_file().
info "Writing gzip-compressed output..."
gzip -1 -c "$FINAL_DUMP" > "$OUTPUT"

output "Dump file written to: ${OUTPUT}"

info "SEXmer-dump complete."
info "  Output     : ${OUTPUT}"
info "  K-mer size : ${KMER_SIZE}"
info "  Min count  : ${MIN_COUNT}"
info "  Trigger-seq: ${TRIGGER_SEQ:-off}"
info "  Total k-mers: ${FILTERED_COUNT}"
