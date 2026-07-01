#!/usr/bin/env bash
# SEXmer-assign.sh - Assign sex from SEXmer dump files using calibrated sex-specific marker k-mers.
# Author: Dede Kurniawan

set -euo pipefail
export LC_ALL=C

# defaults
KMER_SIZE=21
MEM="8G"
THREADS=8
TMPDIR_BASE="$(pwd)"
OUTPUT_PREFIX="sexmer"
MARKERS=""
UNKNOWN_INPUT=""
SAMPLE_INPUT=""
MALE_INPUT=""
FEMALE_INPUT=""
TYPE=""

# log helpers
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*" >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"  >&2; }

usage() {
    cat >&2 <<EOF

SEXmer-assign.sh - Assign sex from SEXmer dump files using calibrated marker k-mers.

Usage: SEXmer-assign.sh <markers.fa> -i <unknown_dumps> -s <samples> -m <male_dumps> -f <female_dumps> --type <XY|ZW> [OPTIONS]

Mandatory:
  <markers.fa>          Marker FASTA from SEXmer-scan or marker/gene sequence FASTA.
                        For XY systems, provide MSK marker FASTA.
                        For ZW systems, provide FSK marker FASTA.
  -i, --input           Comma-separated list of unknown SEXmer dump files (.dump or .dump.gz)
  -s, --sample          Comma-separated list of unknown sample names, same order as --input
  -m, --male            Comma-separated list of known male SEXmer dump files
  -f, --female          Comma-separated list of known female SEXmer dump files
  --type                Sex chromosome system: XY or ZW

Optional:
  -k, --kmer-size       K-mer size used for marker parsing              [default: ${KMER_SIZE}]
  --mem                 Memory budget for sort operations (e.g. 8G)     [default: ${MEM}]
  -t, --threads         CPU threads for sort operations                 [default: ${THREADS}]
  --tmpdir              Parent directory for temporary work folder      [default: current dir]
  -h, --help            Show this help and exit

Output:
  sexmer.assign.txt     Complete assignment report

EOF
    exit 1
}

[[ $# -eq 0 ]] && usage

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)        UNKNOWN_INPUT="$2"; shift 2 ;;
        -s|--sample)       SAMPLE_INPUT="$2";  shift 2 ;;
        -m|--male)         MALE_INPUT="$2";    shift 2 ;;
        -f|--female)       FEMALE_INPUT="$2";  shift 2 ;;
        --type)            TYPE="$2";          shift 2 ;;
        -k|--kmer-size)    KMER_SIZE="$2";     shift 2 ;;
        --mem)             MEM="$2";           shift 2 ;;
        -t|--threads)      THREADS="$2";       shift 2 ;;
        --tmpdir)          TMPDIR_BASE="$2";   shift 2 ;;
        -h|--help)         usage ;;
        -*) error "Unknown option '$1'"; usage ;;
        *)  POSITIONAL+=("$1"); shift ;;
    esac
done

[[ ${#POSITIONAL[@]} -eq 1 ]] || {
    error "Exactly one positional argument is required: <markers.fa>."; usage; }

MARKERS="${POSITIONAL[0]}"

[[ -z "$UNKNOWN_INPUT" ]] && { error "-i/--input files not specified."; usage; }
[[ -z "$SAMPLE_INPUT" ]]  && { error "-s/--sample names not specified."; usage; }
[[ -z "$MALE_INPUT" ]]    && { error "-m/--male files not specified."; usage; }
[[ -z "$FEMALE_INPUT" ]]  && { error "-f/--female files not specified."; usage; }
[[ -z "$TYPE" ]]          && { error "--type not specified."; usage; }

[[ -r "$MARKERS" ]] || { error "Cannot read marker FASTA file: $MARKERS"; exit 1; }

[[ "$KMER_SIZE" =~ ^[1-9][0-9]*$ ]] && [[ "$KMER_SIZE" -le 63 ]] || {
    error "--kmer-size must be an integer between 1 and 63."; exit 1; }

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    error "--threads must be a positive integer."; exit 1; }

case "$TYPE" in
    XY|xy) TYPE="XY" ;;
    ZW|zw) TYPE="ZW" ;;
    *) error "--type must be one of: XY, ZW."; exit 1 ;;
esac

[[ -d "$TMPDIR_BASE" ]] || { error "Temporary parent directory does not exist: $TMPDIR_BASE"; exit 1; }
[[ -w "$TMPDIR_BASE" ]] || { error "Temporary parent directory is not writable: $TMPDIR_BASE"; exit 1; }

for bin in gzip sort join awk python3; do
    command -v "$bin" &>/dev/null || { error "'${bin}' not found on PATH."; exit 1; }
done

if ! sort --version 2>&1 | grep -q 'GNU'; then
    error "GNU sort is required but not found. On macOS: brew install coreutils"
    exit 1
fi

trim_space() {
    local s="$1"
    s="${s#${s%%[![:space:]]*}}"
    s="${s%${s##*[![:space:]]}}"
    echo "$s"
}

split_csv_lines() {
    local input="$1"
    local item
    local RAW_ITEMS=()
    IFS=',' read -ra RAW_ITEMS <<< "$input"
    for item in "${RAW_ITEMS[@]}"; do
        item="$(trim_space "$item")"
        [[ -n "$item" ]] || { error "Empty entry found in comma-separated input: ${input}"; exit 1; }
        printf '%s\n' "$item"
    done
}

strip_dump_suffix() {
    local b="$1"
    b="$(basename "$b")"
    case "$b" in
        *.dump.gz) b="${b%.dump.gz}" ;;
        *.dump)    b="${b%.dump}" ;;
        *.txt.gz)  b="${b%.txt.gz}" ;;
        *.txt)     b="${b%.txt}" ;;
        *.gz)      b="${b%.gz}" ;;
    esac
    echo "$b"
}

open_file() {
    local f="$1"
    if [[ "$f" == *.gz ]]; then
        gzip -dc -- "$f"
    else
        cat -- "$f"
    fi
}

UNKNOWN_FILES=()
UNKNOWN_SAMPLES=()
MALE_FILES=()
FEMALE_FILES=()

while IFS= read -r item; do UNKNOWN_FILES+=("$item"); done < <(split_csv_lines "$UNKNOWN_INPUT")
while IFS= read -r item; do UNKNOWN_SAMPLES+=("$item"); done < <(split_csv_lines "$SAMPLE_INPUT")
while IFS= read -r item; do MALE_FILES+=("$item"); done < <(split_csv_lines "$MALE_INPUT")
while IFS= read -r item; do FEMALE_FILES+=("$item"); done < <(split_csv_lines "$FEMALE_INPUT")

[[ ${#UNKNOWN_FILES[@]} -eq ${#UNKNOWN_SAMPLES[@]} ]] || {
    error "Number of --input files (${#UNKNOWN_FILES[@]}) must match number of --sample names (${#UNKNOWN_SAMPLES[@]})."; exit 1; }

for f in "${UNKNOWN_FILES[@]}" "${MALE_FILES[@]}" "${FEMALE_FILES[@]}"; do
    [[ -r "$f" ]] || { error "Cannot read dump file: $f"; exit 1; }
done

ASSIGN_TMPDIR="${TMPDIR_BASE}/sexmer_assign_tmp_$$"
mkdir -p "$ASSIGN_TMPDIR" "$ASSIGN_TMPDIR/markers" "$ASSIGN_TMPDIR/counts" "$ASSIGN_TMPDIR/samples"
cleanup() { rm -rf "$ASSIGN_TMPDIR"; }
trap cleanup EXIT

REPORT_OUT="${OUTPUT_PREFIX}.assign.txt"

info "SEXmer-assign starting"
info "Parameters: kmer-size=${KMER_SIZE}, type=${TYPE}"
info "Settings  : mem=${MEM}, threads=${THREADS}"
info "Marker file: ${MARKERS}"
info "Temp dir   : ${ASSIGN_TMPDIR}"
info "Unknown samples (${#UNKNOWN_FILES[@]}): ${UNKNOWN_SAMPLES[*]}"
info "Known male files   (${#MALE_FILES[@]}): ${MALE_FILES[*]}"
info "Known female files (${#FEMALE_FILES[@]}): ${FEMALE_FILES[*]}"

MARKER_RAW="${ASSIGN_TMPDIR}/markers/marker_kmers.raw.txt"
MARKER_SORTED="${ASSIGN_TMPDIR}/markers/marker_kmers.sorted.txt"
MARKER_STATS="${ASSIGN_TMPDIR}/markers/marker_stats.txt"

info "Parsing marker FASTA and building marker k-mer index..."

python3 - "$MARKERS" "$KMER_SIZE" "$MARKER_RAW" "$MARKER_STATS" <<'PY'
import gzip
import sys
from pathlib import Path

markers = sys.argv[1]
k = int(sys.argv[2])
out_path = sys.argv[3]
stats_path = sys.argv[4]

def open_text(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt')

def parse_fasta(path):
    name = None
    chunks = []
    with open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, ''.join(chunks)

records = 0
short_records = 0
bad_windows = 0
exact_records = 0
long_records = 0
raw_kmers = 0
min_len = 0
max_len = 0

with open(out_path, 'w') as out:
    for name, seq in parse_fasta(markers):
        records += 1
        seq = ''.join(seq.split()).upper()
        length = len(seq)
        if min_len == 0 or length < min_len:
            min_len = length
        if length > max_len:
            max_len = length
        if length < k:
            short_records += 1
            continue
        if length == k:
            exact_records += 1
        else:
            long_records += 1
        for i in range(0, length - k + 1):
            kmer = seq[i:i+k]
            if set(kmer) <= {'A', 'C', 'G', 'T'}:
                out.write(kmer + '\n')
                raw_kmers += 1
            else:
                bad_windows += 1

with open(stats_path, 'w') as stats:
    stats.write(f"records\t{records}\n")
    stats.write(f"exact_kmer_records\t{exact_records}\n")
    stats.write(f"long_records\t{long_records}\n")
    stats.write(f"short_records\t{short_records}\n")
    stats.write(f"raw_marker_kmers\t{raw_kmers}\n")
    stats.write(f"bad_windows\t{bad_windows}\n")
    stats.write(f"min_length\t{min_len}\n")
    stats.write(f"max_length\t{max_len}\n")

if records == 0:
    print("[Error] Marker FASTA contains no records.", file=sys.stderr)
    sys.exit(2)
if raw_kmers == 0:
    print("[Error] No valid marker k-mers were generated from marker FASTA.", file=sys.stderr)
    sys.exit(3)
PY

sort -u -S "$MEM" --parallel="$THREADS" -T "$ASSIGN_TMPDIR" "$MARKER_RAW" > "$MARKER_SORTED"
MARKER_TOTAL=$(wc -l < "$MARKER_SORTED")

MARKER_RECORDS=$(awk -F'\t' '$1 == "records" { print $2 }' "$MARKER_STATS")
MARKER_EXACT=$(awk -F'\t' '$1 == "exact_kmer_records" { print $2 }' "$MARKER_STATS")
MARKER_LONG=$(awk -F'\t' '$1 == "long_records" { print $2 }' "$MARKER_STATS")
MARKER_SHORT=$(awk -F'\t' '$1 == "short_records" { print $2 }' "$MARKER_STATS")
MARKER_BAD=$(awk -F'\t' '$1 == "bad_windows" { print $2 }' "$MARKER_STATS")
MARKER_MINLEN=$(awk -F'\t' '$1 == "min_length" { print $2 }' "$MARKER_STATS")
MARKER_MAXLEN=$(awk -F'\t' '$1 == "max_length" { print $2 }' "$MARKER_STATS")

info "Marker k-mer index built successfully."
info "  FASTA records       : ${MARKER_RECORDS}"
info "  Exact k-mer records : ${MARKER_EXACT}"
info "  Long records        : ${MARKER_LONG}"
info "  Short records skipped: ${MARKER_SHORT}"
info "  Non-ACGT windows skipped: ${MARKER_BAD}"
info "  Marker length range : ${MARKER_MINLEN}-${MARKER_MAXLEN} bp"
info "  Unique marker k-mers: ${MARKER_TOTAL}"

if [[ "$MARKER_SHORT" -gt 0 ]]; then
    warn "${MARKER_SHORT} marker FASTA record(s) shorter than k=${KMER_SIZE} were skipped."
fi

COUNT_TSV="${ASSIGN_TMPDIR}/all_sample_marker_ratios.tsv"
printf "sample\tgroup\tdump_file\tsample_unique_kmers\tsample_total_kmer_count\tdetected_marker_kmers\ttotal_marker_kmers\tmarker_ratio\tmarker_kmer_hits\tmarker_hit_rate\n" > "$COUNT_TSV"

count_marker_ratio() {
    local sample="$1"
    local group="$2"
    local dump_file="$3"
    local safe idx sorted joined sample_unique sample_total_count detected marker_hits ratio hit_rate

    safe="${sample//[^A-Za-z0-9_.-]/_}_${group}"
    sorted="${ASSIGN_TMPDIR}/samples/${safe}.sorted_counts.txt"
    joined="${ASSIGN_TMPDIR}/counts/${safe}.marker_join.txt"

    info "Counting marker k-mers in ${group} sample '${sample}'..."

    open_file "$dump_file" \
    | awk '
        NF >= 1 {
            k = toupper($1)
            c = (NF >= 2 && $2 ~ /^[0-9]+$/ ? $2 : 1)
            if (k ~ /^[ACGT]+$/) print k, c
        }' \
    | sort -k1,1 -S "$MEM" --parallel="$THREADS" -T "$ASSIGN_TMPDIR" \
    | awk '
        BEGIN { prev=""; sum=0 }
        {
            if ($1 == prev) { sum += $2 }
            else {
                if (prev != "") print prev, sum
                prev = $1; sum = $2
            }
        }
        END { if (prev != "") print prev, sum }
    ' > "$sorted"

    sample_unique=$(wc -l < "$sorted")
    sample_total_count=$(awk '{ s += $2 } END { print s+0 }' "$sorted")

    join "$MARKER_SORTED" "$sorted" > "$joined" || true

    read -r detected marker_hits < <(awk '
        BEGIN { n=0; h=0 }
        { n++; h += $2 }
        END { print n+0, h+0 }
    ' "$joined")

    ratio=$(awk -v d="$detected" -v t="$MARKER_TOTAL" 'BEGIN { if (t > 0) printf "%.10f", (d / t) * 100.0; else printf "0" }')
    hit_rate=$(awk -v h="$marker_hits" -v t="$sample_total_count" 'BEGIN { if (t > 0) printf "%.10g", h / t; else printf "0" }')

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" "$group" "$dump_file" "$sample_unique" "$sample_total_count" "$detected" "$MARKER_TOTAL" "$ratio" "$marker_hits" "$hit_rate" \
        >> "$COUNT_TSV"

    info "  ${sample}: detected=${detected}/${MARKER_TOTAL}, marker_ratio=${ratio}%"
}

for f in "${MALE_FILES[@]}"; do
    count_marker_ratio "$(strip_dump_suffix "$f")" "known_male" "$f"
done

for f in "${FEMALE_FILES[@]}"; do
    count_marker_ratio "$(strip_dump_suffix "$f")" "known_female" "$f"
done

for idx in "${!UNKNOWN_FILES[@]}"; do
    count_marker_ratio "${UNKNOWN_SAMPLES[$idx]}" "unknown" "${UNKNOWN_FILES[$idx]}"
done

info "Generating calibrated assignment report..."

python3 - "$COUNT_TSV" "$REPORT_OUT" "$MARKERS" "$TYPE" "$KMER_SIZE" "$MARKER_TOTAL" \
    "$MARKER_RECORDS" "$MARKER_EXACT" "$MARKER_LONG" "$MARKER_SHORT" "$MARKER_BAD" <<'PY'
import csv
import statistics
import sys
from datetime import datetime

count_tsv, report_out, markers, sex_type, kmer_size, marker_total = sys.argv[1:7]
marker_records, marker_exact, marker_long, marker_short, marker_bad = sys.argv[7:12]

rows = []
with open(count_tsv, newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    for row in reader:
        row['marker_ratio_float'] = float(row['marker_ratio'])
        rows.append(row)

male_rows = [r for r in rows if r['group'] == 'known_male']
female_rows = [r for r in rows if r['group'] == 'known_female']
unknown_rows = [r for r in rows if r['group'] == 'unknown']

if not male_rows or not female_rows or not unknown_rows:
    print('[Error] Missing known male, known female, or unknown rows for assignment.', file=sys.stderr)
    sys.exit(1)

male_median = statistics.median([r['marker_ratio_float'] for r in male_rows])
female_median = statistics.median([r['marker_ratio_float'] for r in female_rows])
threshold = (male_median + female_median) / 2.0
separation = abs(male_median - female_median)
margin = separation * 0.10
amb_low = min(threshold - margin, threshold + margin)
amb_high = max(threshold - margin, threshold + margin)

expected_direction_pass = True
expected_msg = 'PASS'
if sex_type == 'XY':
    if male_median <= female_median:
        expected_direction_pass = False
        expected_msg = 'WARNING: XY expects known male marker ratio > known female marker ratio. Check whether markers.fa is MSK.fa.'
else:
    if female_median <= male_median:
        expected_direction_pass = False
        expected_msg = 'WARNING: ZW expects known female marker ratio > known male marker ratio. Check whether markers.fa is FSK.fa.'

warnings = []
if not expected_direction_pass:
    warnings.append(expected_msg)
if separation < 20.0:
    warnings.append('Known male/female marker-ratio separation is < 20 percentage points; assignments may be weak.')

def assign_ratio(ratio):
    if amb_low <= ratio <= amb_high:
        return 'ambiguous', 'low', 'marker ratio falls inside calibrated ambiguous range'
    if sex_type == 'XY':
        if ratio > threshold:
            sex = 'male'
            reason = 'marker ratio is above calibrated threshold for XY/MSK expectation'
        else:
            sex = 'female'
            reason = 'marker ratio is below calibrated threshold for XY/MSK expectation'
    else:
        if ratio > threshold:
            sex = 'female'
            reason = 'marker ratio is above calibrated threshold for ZW/FSK expectation'
        else:
            sex = 'male'
            reason = 'marker ratio is below calibrated threshold for ZW/FSK expectation'
    confidence = 'high' if expected_direction_pass and separation >= 20.0 else 'low'
    return sex, confidence, reason

for row in rows:
    sex, confidence, reason = assign_ratio(row['marker_ratio_float'])
    row['assignment'] = sex
    row['confidence'] = confidence
    row['reason'] = reason

def fmt(x):
    return f'{float(x):.4f}'

with open(report_out, 'w') as out:
    out.write('SEXmer-assign result\n')
    out.write('\n')
    out.write('Run information\n')
    out.write(f'Marker file          : {markers}\n')
    out.write(f'System type          : {sex_type}\n')
    out.write('Expected marker      : MSK for XY, FSK for ZW\n')
    out.write(f'K-mer size           : {kmer_size}\n')
    out.write(f'Unknown samples      : {len(unknown_rows)}\n')
    out.write(f'Known male samples   : {len(male_rows)}\n')
    out.write(f'Known female samples : {len(female_rows)}\n')
    out.write('\n')

    out.write('Marker statistics\n')
    out.write(f'FASTA records                : {marker_records}\n')
    out.write(f'Exact k-mer records          : {marker_exact}\n')
    out.write(f'Long marker/gene records     : {marker_long}\n')
    out.write(f'Short records skipped        : {marker_short}\n')
    out.write(f'Non-ACGT marker windows skipped : {marker_bad}\n')
    out.write(f'Total unique marker k-mers   : {marker_total}\n')
    out.write('\n')

    out.write('Calibration\n')
    out.write(f'Known male median ratio      : {fmt(male_median)}\n')
    out.write(f'Known female median ratio    : {fmt(female_median)}\n')
    out.write(f'Male/female separation       : {fmt(separation)}\n')
    out.write(f'Decision threshold           : {fmt(threshold)}\n')
    out.write(f'Ambiguous lower bound        : {fmt(amb_low)}\n')
    out.write(f'Ambiguous upper bound        : {fmt(amb_high)}\n')
    out.write(f'Direction check              : {expected_msg}\n')
    out.write('\n')

    if warnings:
        out.write('Warnings\n')
        for warning in warnings:
            out.write(f'- {warning}\n')
        out.write('\n')

    out.write('Known sample ratios\n')
    out.write('sample\tgroup\tsample_unique_kmers\tsample_total_kmer_count\tdetected_marker_kmers\ttotal_marker_kmers\tmarker_ratio\tmarker_kmer_hits\tmarker_hit_rate\tassignment\tconfidence\n')
    for row in male_rows + female_rows:
        out.write('\t'.join([
            row['sample'], row['group'], row['sample_unique_kmers'], row['sample_total_kmer_count'],
            row['detected_marker_kmers'], row['total_marker_kmers'], fmt(row['marker_ratio_float']), row['marker_kmer_hits'],
            row['marker_hit_rate'], row['assignment'], row['confidence']
        ]) + '\n')
    out.write('\n')

    out.write('Unknown sample assignment\n')
    out.write('sample\tsample_unique_kmers\tsample_total_kmer_count\tdetected_marker_kmers\ttotal_marker_kmers\tmarker_ratio\tmarker_kmer_hits\tmarker_hit_rate\tassignment\tconfidence\treason\n')
    for row in unknown_rows:
        out.write('\t'.join([
            row['sample'], row['sample_unique_kmers'], row['sample_total_kmer_count'], row['detected_marker_kmers'],
            row['total_marker_kmers'], fmt(row['marker_ratio_float']), row['marker_kmer_hits'],
            row['marker_hit_rate'], row['assignment'], row['confidence'], row['reason']
        ]) + '\n')
PY

if grep -q '^Direction check              : WARNING' "$REPORT_OUT"; then
    warn "Marker-ratio direction check failed. See ${REPORT_OUT}."
fi

output "Assignment report written to: ${REPORT_OUT}"

info "SEXmer-assign complete."
info "  Output: ${REPORT_OUT}"
info "  Marker k-mers: ${MARKER_TOTAL}"
info "  Unknown samples: ${#UNKNOWN_FILES[@]}"
