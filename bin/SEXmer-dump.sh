#!/usr/bin/env bash
# SEXmer-dump.sh - Generate a filtered k-mer dump file for a single sample.
# Author: Dede Kurniawan

set -euo pipefail
export LC_ALL=C

KMER_SIZE=21
MIN_COUNT=3
TRIGGER_SEQ=""
THREADS=4
PREFIX=""
OUTDIR="."
TMPDIR_BASE="$(pwd)"
KMC_BIN=""
READS=()

info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*"  >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"   >&2; }

usage() {
    cat <<EOF

SEXmer-dump.sh - Generate a filtered k-mer dump file for a single sample.

Usage: SEXmer-dump.sh --prefix <sample> <reads_1.fq.gz> [reads_2.fq.gz] [OPTIONS]

Mandatory:
  --prefix             Output filename prefix (output: <prefix>.dump.gz)
  <reads>              One or two FASTQ files (.fq, .fastq, .fq.gz, .fastq.gz)

Optional:
  -k, --kmer-size      K-mer size (1-63)                              [default: ${KMER_SIZE}]
  --min-count          Minimum k-mer count (KMC -ci)                  [default: ${MIN_COUNT}]
  --trigger-seq        Retain only k-mers starting with this seq      [default: off]
  -t, --threads        CPU threads                                    [default: ${THREADS}]
  -o, --outdir         Output directory                               [default: current dir]
  --tmpdir             Parent directory for the temporary work folder [default: current dir]
  --kmc-bin            Directory containing kmc and kmc_dump binaries [default: PATH]
  -h, --help           Show this help and exit

EOF
}

[[ $# -eq 0 ]] && { usage >&2; exit 1; }

# argument parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)         PREFIX="$2";        shift 2 ;;
        -k|--kmer-size)   KMER_SIZE="$2";    shift 2 ;;
        --min-count)      MIN_COUNT="$2";    shift 2 ;;
        --trigger-seq)    TRIGGER_SEQ="$2";  shift 2 ;;
        -t|--threads)     THREADS="$2";      shift 2 ;;
        -o|--outdir)      OUTDIR="$2";       shift 2 ;;
        --tmpdir)         TMPDIR_BASE="$2";  shift 2 ;;
        --kmc-bin)        KMC_BIN="$2";      shift 2 ;;
        -h|--help)        usage; exit 0 ;;
        -*) error "Unknown option '$1'"; usage >&2; exit 1 ;;
        *)  READS+=("$1"); shift ;;
    esac
done

[[ -z "$PREFIX" ]] && { error "--prefix is required."; usage >&2; exit 1; }

[[ ${#READS[@]} -eq 0 ]] && { error "At least one input read file is required."; usage >&2; exit 1; }

[[ ${#READS[@]} -gt 2 ]] && {
    error "At most two read files (paired-end) are accepted. Got: ${READS[*]}"; exit 1; }

[[ "$KMER_SIZE" =~ ^[1-9][0-9]*$ ]] && [[ "$KMER_SIZE" -le 63 ]] || {
    error "--kmer-size must be an integer between 1 and 63."; exit 1; }

[[ "$MIN_COUNT" =~ ^[1-9][0-9]*$ ]] || {
    error "--min-count must be a positive integer."; exit 1; }

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    error "--threads must be a positive integer."; exit 1; }

if [[ -n "$TRIGGER_SEQ" ]]; then
    [[ "$TRIGGER_SEQ" =~ ^[ACGTacgt]+$ ]] || {
        error "--trigger-seq must contain only ACGT characters. Got: '${TRIGGER_SEQ}'"; exit 1; }
    TRIGGER_SEQ="${TRIGGER_SEQ^^}"
fi

if [[ -n "$KMC_BIN" ]]; then
    [[ -d "$KMC_BIN" ]] || { error "--kmc-bin directory does not exist: ${KMC_BIN}"; exit 1; }
    export PATH="${KMC_BIN}:${PATH}"
fi

[[ -d "$TMPDIR_BASE" ]] || { error "Temporary parent directory does not exist: ${TMPDIR_BASE}"; exit 1; }
[[ -w "$TMPDIR_BASE" ]] || { error "Temporary parent directory is not writable: ${TMPDIR_BASE}"; exit 1; }

OUTPUT="${OUTDIR%/}/${PREFIX}.dump.gz"
OUT_DIR="$(dirname "$OUTPUT")"
mkdir -p "$OUT_DIR" || { error "Could not create output directory: ${OUT_DIR}"; exit 1; }
[[ -d "$OUT_DIR" ]] || { error "Output path is not a directory: ${OUT_DIR}"; exit 1; }
[[ -w "$OUT_DIR" ]] || { error "Output directory is not writable: ${OUT_DIR}"; exit 1; }

for f in "${READS[@]}"; do
    [[ -r "$f" ]] || { error "Cannot read input file: $f"; exit 1; }
done

for bin in kmc kmc_dump; do
    command -v "$bin" &>/dev/null || {
        error "'${bin}' not found on PATH. Install KMC: conda install bioconda::kmc"
        exit 1; }
done

DUMP_TMPDIR="${TMPDIR_BASE}/sexmer_dump_tmp_$$"
mkdir -p "$DUMP_TMPDIR"
cleanup() { rm -rf "$DUMP_TMPDIR"; }
trap cleanup EXIT

info "SEXmer-dump starting"
info "Parameters : kmer-size=${KMER_SIZE}, min-count=${MIN_COUNT}, threads=${THREADS}"
info "Trigger-seq: ${TRIGGER_SEQ:-off}"
info "Input reads: ${READS[*]}"
info "Output     : ${OUTPUT}"
info "Temp dir   : ${DUMP_TMPDIR}"

# STEP 1: Build KMC file-of-files
info "Building KMC input file list..."

FILES_LIST="${DUMP_TMPDIR}/input_files.lst"
printf '%s\n' "${READS[@]}" > "$FILES_LIST"

info "  ${#READS[@]} file(s) listed."

# STEP 2: Run KMC
KMC_DB="${DUMP_TMPDIR}/kmc_db"
KMC_TMP="${DUMP_TMPDIR}/kmc_tmp"
mkdir -p "$KMC_TMP"

info "Running KMC (k=${KMER_SIZE}, min-count=${MIN_COUNT}, threads=${THREADS})..."

kmc \
    -k"${KMER_SIZE}" \
    -ci"${MIN_COUNT}" \
    -t"${THREADS}" \
    @"${FILES_LIST}" \
    "${KMC_DB}" \
    "${KMC_TMP}" \
    2>&1 | while IFS= read -r line; do info "  [kmc] $line"; done

[[ -f "${KMC_DB}.kmc_pre" && -f "${KMC_DB}.kmc_suf" ]] || {
    error "KMC did not produce expected database files. Check log above."; exit 1; }

info "KMC database built successfully."

# STEP 3: Dump k-mers to text
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
info "Writing gzip-compressed output..."
gzip -1 -c "$FINAL_DUMP" > "$OUTPUT"

output "Dump file written to: ${OUTPUT}"

info "SEXmer-dump complete."
info "  Output     : ${OUTPUT}"
info "  K-mer size : ${KMER_SIZE}"
info "  Min count  : ${MIN_COUNT}"
info "  Trigger-seq: ${TRIGGER_SEQ:-off}"
info "  Total k-mers: ${FILTERED_COUNT}"
