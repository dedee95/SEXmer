#!/usr/bin/env bash
# SEXmer reads - Extract specific reads based on specific kmer sequence.
# Author: Dede Kurniawan

set -euo pipefail
export LC_ALL=C

# defaults
KMER_SIZE=21
MIN_HIT=3
THREADS=8
PREFIX=""
OUTDIR="."
TMPDIR_BASE="$(pwd)"
BBDUK_BIN=""
MARKERS=""
READS_INPUT=""
READS=()

# log helpers
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*"  >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"   >&2; }

usage() {
    cat <<EOF_USAGE

SEXmer reads - Extract specific reads based on specific kmer sequence.

Usage: SEXmer reads <markers.fa> -r <reads> --prefix <prefix> [OPTIONS]

Mandatory:
  <markers.fa>         Marker k-mer FASTA, e.g. MSK.fa or FSK.fa (.gz accepted)
  -r, --reads          Comma-separated FASTQ file(s)
                       Example: -r reads_1.fq.gz,reads_2.fq.gz or 
                                -r long_reads.fq.gz
  --prefix             Output filename prefix

Optionals:
  --hit                Minimum exact marker k-mer hits per read        [default: ${MIN_HIT}]
  -k, --kmer-size      K-mer size used by BBDuk (1-63)                 [default: ${KMER_SIZE}]
  -t, --threads        CPU threads                                     [default: ${THREADS}]
  -o, --outdir         Output directory                                [default: current dir]
  --tmpdir             Parent directory for temporary work folder      [default: current dir]
  --bbduk-bin          Directory containing bbduk.sh                   [default: PATH]
  -h, --help           Show this help and exit

EOF_USAGE
}
[[ $# -eq 0 ]] && { usage >&2; exit 1; }

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--reads)        READS_INPUT="$2"; shift 2 ;;
        --prefix)          PREFIX="$2";      shift 2 ;;
        -k|--kmer-size)    KMER_SIZE="$2";  shift 2 ;;
        --hit)             MIN_HIT="$2";     shift 2 ;;
        -t|--threads)      THREADS="$2";     shift 2 ;;
        -o|--outdir)       OUTDIR="$2";      shift 2 ;;
        --tmpdir)          TMPDIR_BASE="$2"; shift 2 ;;
        --bbduk-bin)       BBDUK_BIN="$2";   shift 2 ;;
        -h|--help)         usage; exit 0 ;;
        -*) error "Unknown option '$1'"; usage >&2; exit 1 ;;
        *)  POSITIONAL+=("$1"); shift ;;
    esac
done

[[ ${#POSITIONAL[@]} -eq 1 ]] || {
    error "Exactly one positional argument is required: <markers.fa>."; usage >&2; exit 1; }

MARKERS="${POSITIONAL[0]}"

[[ -z "$READS_INPUT" ]] && { error "-r/--reads is required."; usage >&2; exit 1; }
[[ -z "$PREFIX" ]] && { error "--prefix is required."; usage >&2; exit 1; }

[[ -r "$MARKERS" ]] || { error "Cannot read marker FASTA file: $MARKERS"; exit 1; }

[[ "$KMER_SIZE" =~ ^[1-9][0-9]*$ ]] && [[ "$KMER_SIZE" -le 63 ]] || {
    error "--kmer-size must be an integer between 1 and 63."; exit 1; }

[[ "$MIN_HIT" =~ ^[1-9][0-9]*$ ]] || {
    error "--hit must be a positive integer."; exit 1; }

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    error "--threads must be a positive integer."; exit 1; }

if [[ -n "$BBDUK_BIN" ]]; then
    [[ -d "$BBDUK_BIN" ]] || { error "--bbduk-bin directory does not exist: ${BBDUK_BIN}"; exit 1; }
    export PATH="${BBDUK_BIN}:${PATH}"
fi

[[ -d "$TMPDIR_BASE" ]] || { error "Temporary parent directory does not exist: ${TMPDIR_BASE}"; exit 1; }
[[ -w "$TMPDIR_BASE" ]] || { error "Temporary parent directory is not writable: ${TMPDIR_BASE}"; exit 1; }

mkdir -p "$OUTDIR" || { error "Cannot create output directory: ${OUTDIR}"; exit 1; }
[[ -d "$OUTDIR" ]] || { error "Output path is not a directory: ${OUTDIR}"; exit 1; }
[[ -w "$OUTDIR" ]] || { error "Output directory is not writable: ${OUTDIR}"; exit 1; }
if [[ "$OUTDIR" != "." ]]; then
    PREFIX="${OUTDIR%/}/${PREFIX}"
fi
OUT_DIR="$(dirname "$PREFIX")"
mkdir -p "$OUT_DIR" || { error "Cannot create output directory: ${OUT_DIR}"; exit 1; }
[[ -w "$OUT_DIR" ]] || { error "Output directory is not writable: ${OUT_DIR}"; exit 1; }

command -v bbduk.sh &>/dev/null || {
    error "bbduk.sh not found on PATH. Install BBTools/BBMap or provide --bbduk-bin."
    exit 1; }
command -v gzip &>/dev/null || { error "gzip not found on PATH."; exit 1; }
command -v awk &>/dev/null || { error "awk not found on PATH."; exit 1; }

trim_space() {
    local s="$1"
    s="${s#${s%%[![:space:]]*}}"
    s="${s%${s##*[![:space:]]}}"
    echo "$s"
}

IFS=',' read -ra RAW_READS <<< "$READS_INPUT"
for item in "${RAW_READS[@]}"; do
    item="$(trim_space "$item")"
    [[ -n "$item" ]] || { error "Empty entry found in -r/--reads comma-separated list."; exit 1; }
    READS+=("$item")
done

[[ ${#READS[@]} -gt 0 ]] || { error "At least one input read file is required."; exit 1; }

for f in "${READS[@]}"; do
    [[ -r "$f" ]] || { error "Cannot read input FASTQ file: $f"; exit 1; }
done

READS_TMPDIR="${TMPDIR_BASE}/sexmer_reads_tmp_$$"
mkdir -p "$READS_TMPDIR"
cleanup() { rm -rf "$READS_TMPDIR"; }
trap cleanup EXIT

open_file() {
    local f="$1"
    if [[ "$f" == *.gz ]]; then
        gzip -dc -- "$f"
    else
        cat -- "$f"
    fi
}

stream_files() {
    local f
    for f in "$@"; do
        open_file "$f"
    done
}

strip_fastq_suffix() {
    local b="$1"
    case "$b" in
        *.fastq.gz) b="${b%.fastq.gz}" ;;
        *.fq.gz)    b="${b%.fq.gz}" ;;
        *.fastq)    b="${b%.fastq}" ;;
        *.fq)       b="${b%.fq}" ;;
        *.gz)       b="${b%.gz}" ;;
    esac
    echo "$b"
}

filename_pair_score() {
    local p1="$1" p2="$2"
    local b1 b2 stem1 stem2 expected
    b1="$(basename "$p1")"
    b2="$(basename "$p2")"
    stem1="$(strip_fastq_suffix "$b1")"
    stem2="$(strip_fastq_suffix "$b2")"

    if [[ "$stem1" =~ ^(.+)_R1_001$ ]]; then
        expected="${BASH_REMATCH[1]}_R2_001"
        [[ "$stem2" == "$expected" ]] && return 0
    fi
    if [[ "$stem1" =~ ^(.+)_R1$ ]]; then
        expected="${BASH_REMATCH[1]}_R2"
        [[ "$stem2" == "$expected" ]] && return 0
    fi
    if [[ "$stem1" =~ ^(.+)_1$ ]]; then
        expected="${BASH_REMATCH[1]}_2"
        [[ "$stem2" == "$expected" ]] && return 0
    fi
    if [[ "$stem1" =~ ^(.+)\.1$ ]]; then
        expected="${BASH_REMATCH[1]}.2"
        [[ "$stem2" == "$expected" ]] && return 0
    fi
    if [[ "$stem1" =~ ^(.+)-1$ ]]; then
        expected="${BASH_REMATCH[1]}-2"
        [[ "$stem2" == "$expected" ]] && return 0
    fi
    return 1
}

normalize_read_id() {
    local h="$1"
    h="${h%%[[:space:]]*}"
    h="${h#@}"
    h="${h%/1}"
    h="${h%/2}"
    echo "$h"
}

header_pair_score() {
    local r1="$1" r2="$2" max_records="${3:-1000}"
    local h1 s1 p1 q1 h2 s2 p2 q2 id1 id2 total same
    total=0
    same=0

    exec 3< <(open_file "$r1")
    exec 4< <(open_file "$r2")

    while [[ "$total" -lt "$max_records" ]]; do
        if ! IFS= read -r h1 <&3; then
            break
        fi
        IFS= read -r s1 <&3 || break
        IFS= read -r p1 <&3 || break
        IFS= read -r q1 <&3 || break

        if ! IFS= read -r h2 <&4; then
            break
        fi
        IFS= read -r s2 <&4 || break
        IFS= read -r p2 <&4 || break
        IFS= read -r q2 <&4 || break

        total=$(( total + 1 ))
        id1="$(normalize_read_id "$h1")"
        id2="$(normalize_read_id "$h2")"
        [[ "$id1" == "$id2" ]] && same=$(( same + 1 ))
    done

    exec 3<&-
    exec 4<&-

    HEADER_SAME="$same"
    HEADER_TOTAL="$total"

    [[ "$total" -gt 0 ]] || return 1
    awk -v s="$same" -v t="$total" 'BEGIN { exit ((s / t) >= 0.90 ? 0 : 1) }'
}

validate_marker_fasta() {
    local markers="$1" k="$2" stats_out="$3"

    open_file "$markers" | awk -v k="$k" '
        BEGIN {
            seq = ""
            total = 0
            bad = 0
            wrong = 0
            minlen = 0
            maxlen = 0
        }
        function check_record(s, len) {
            if (s == "") return
            s = toupper(s)
            len = length(s)
            total++
            if (s !~ /^[ACGT]+$/) bad++
            if (len != k) wrong++
            if (minlen == 0 || len < minlen) minlen = len
            if (len > maxlen) maxlen = len
        }
        /^>/ {
            check_record(seq)
            seq = ""
            next
        }
        {
            gsub(/[[:space:]]/, "")
            seq = seq $0
        }
        END {
            check_record(seq)
            print "markers\t" total
            print "bad_bases\t" bad
            print "wrong_length\t" wrong
            print "min_length\t" minlen
            print "max_length\t" maxlen
            if (total == 0) exit 2
            if (bad > 0) exit 3
            if (wrong > 0) exit 4
        }
    ' > "$stats_out"
}

run_bbduk_logged() {
    local log_file="$1"
    shift

    "$@" 2> "$log_file"
    while IFS= read -r line; do
        info "  [bbduk] $line"
    done < "$log_file"
}

info "SEXmer reads starting"
info "Parameters : kmer-size=${KMER_SIZE}, hit=${MIN_HIT}, threads=${THREADS}"
info "Marker file: ${MARKERS}"
info "Input reads: ${READS[*]}"
info "Backend    : BBTools BBDuk"
info "Temp dir   : ${READS_TMPDIR}"

MARKER_STATS="${READS_TMPDIR}/marker_stats.txt"
set +e
validate_marker_fasta "$MARKERS" "$KMER_SIZE" "$MARKER_STATS"
MARKER_STATUS=$?
set -e
if [[ "$MARKER_STATUS" -ne 0 ]]; then
    case "$MARKER_STATUS" in
        2) error "Marker FASTA contains no sequences: ${MARKERS}" ;;
        3) error "Marker FASTA contains non-ACGT bases: ${MARKERS}" ;;
        4) error "Marker FASTA contains sequence(s) with length different from --kmer-size ${KMER_SIZE}." ;;
        *) error "Failed to validate marker FASTA: ${MARKERS}" ;;
    esac
    exit 1
fi

MARKER_COUNT=$(awk -F'\t' '$1 == "markers" { print $2 }' "$MARKER_STATS")
BAD_BASES=$(awk -F'\t' '$1 == "bad_bases" { print $2 }' "$MARKER_STATS")
WRONG_LENGTH=$(awk -F'\t' '$1 == "wrong_length" { print $2 }' "$MARKER_STATS")
MIN_LEN=$(awk -F'\t' '$1 == "min_length" { print $2 }' "$MARKER_STATS")
MAX_LEN=$(awk -F'\t' '$1 == "max_length" { print $2 }' "$MARKER_STATS")

info "Marker FASTA validation complete."
info "  Marker sequences: ${MARKER_COUNT}"
info "  Marker length   : ${MIN_LEN}-${MAX_LEN} bp"
info "  Non-ACGT records: ${BAD_BASES}"
info "  Wrong length    : ${WRONG_LENGTH}"

MODE="single"
FILENAME_PAIR="no"
HEADER_PAIR="no"
HEADER_SAME=0
HEADER_TOTAL=0

if [[ ${#READS[@]} -eq 2 ]]; then
    if filename_pair_score "${READS[0]}" "${READS[1]}"; then
        FILENAME_PAIR="yes"
    fi
    if header_pair_score "${READS[0]}" "${READS[1]}" 1000; then
        HEADER_PAIR="yes"
    fi
    if [[ "$FILENAME_PAIR" == "yes" || "$HEADER_PAIR" == "yes" ]]; then
        MODE="paired"
    fi
fi

BBDUK_LOG="${READS_TMPDIR}/bbduk.log"

if [[ "$MODE" == "paired" ]]; then
    info "Read mode detected: paired-end"
    info "  Filename pair pattern: ${FILENAME_PAIR}"
    info "  Header ID agreement  : ${HEADER_SAME}/${HEADER_TOTAL} sampled read pairs"
    info "  Pair rule            : write both mates if either mate has >= ${MIN_HIT} exact marker hit(s)"

    TMP_OUT1="${READS_TMPDIR}/matched_1.fq"
    TMP_OUT2="${READS_TMPDIR}/matched_2.fq"
    OUT1="${PREFIX}.sexmer_1.fq.gz"
    OUT2="${PREFIX}.sexmer_2.fq.gz"

    info "Running BBDuk paired-end extraction..."
    run_bbduk_logged "$BBDUK_LOG" \
        bbduk.sh \
            in="${READS[0]}" \
            in2="${READS[1]}" \
            outm="$TMP_OUT1" \
            outm2="$TMP_OUT2" \
            ref="$MARKERS" \
            k="$KMER_SIZE" \
            hdist=0 \
            minhits="$MIN_HIT" \
            threads="$THREADS"

    [[ -f "$TMP_OUT1" ]] || { error "BBDuk did not produce R1 matched output."; exit 1; }
    [[ -f "$TMP_OUT2" ]] || { error "BBDuk did not produce R2 matched output."; exit 1; }

    info "Writing gzip-compressed output with gzip -1..."
    gzip -1 -c "$TMP_OUT1" > "$OUT1"
    gzip -1 -c "$TMP_OUT2" > "$OUT2"

    output "R1 reads written to: ${OUT1}"
    output "R2 reads written to: ${OUT2}"
    info "SEXmer reads complete."
    info "  Mode     : paired-end"
    info "  Output R1: ${OUT1}"
    info "  Output R2: ${OUT2}"
else
    OUT="${PREFIX}.sexmer.fq.gz"
    info "Read mode detected: single/long-read"
    if [[ ${#READS[@]} -eq 2 ]]; then
        info "  Filename pair pattern: ${FILENAME_PAIR}"
        info "  Header ID agreement  : ${HEADER_SAME}/${HEADER_TOTAL} sampled records"
        error "Two input files were provided but did not look like a paired-end R1/R2 set."
        error "For long-read or unpaired mode, provide only one FASTQ file per SEXmer reads run."
        exit 1
    fi

    if [[ ${#READS[@]} -gt 1 ]]; then
        error "Multiple unpaired/long-read files are not supported in one SEXmer reads run."
        error "Please run SEXmer reads separately for each long-read FASTQ file."
        exit 1
    fi

    info "  Input files: ${#READS[@]}"
    info "  Read rule  : write read if it has >= ${MIN_HIT} exact marker hit(s)"

    TMP_OUT="${READS_TMPDIR}/matched.fq"
    info "Running BBDuk single-file extraction..."
    run_bbduk_logged "$BBDUK_LOG" \
        bbduk.sh \
            in="${READS[0]}" \
            outm="$TMP_OUT" \
            ref="$MARKERS" \
            k="$KMER_SIZE" \
            hdist=0 \
            minhits="$MIN_HIT" \
            threads="$THREADS"

    [[ -f "$TMP_OUT" ]] || { error "BBDuk did not produce matched read output."; exit 1; }

    info "Writing gzip-compressed output with gzip -1..."
    gzip -1 -c "$TMP_OUT" > "$OUT"

    [[ -f "$OUT" ]] || { error "Expected output file was not created: ${OUT}"; exit 1; }

    output "Reads written to: ${OUT}"
    info "SEXmer reads complete."
    info "  Mode  : single/long-read"
    info "  Output: ${OUT}"
fi
