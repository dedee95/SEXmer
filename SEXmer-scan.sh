#!/usr/bin/env bash
# SEXmer-scan.sh - K-mer sex-specificity classifier from KMC or Jellyfish dump files.
# Author: Dede Kurniawan 

set -euo pipefail
export LC_ALL=C

# defaults
MIN_COUNT=10
MAX_COUNT=500
FOLD_THRESHOLD=5
PREFIX="output"
MEM="8G"
THREADS=8
MALE_INPUT=""
FEMALE_INPUT=""
NEUTRAL_MAX=100000
TMPDIR_BASE="$(pwd)"
SEED=42

# log helpers
# All log output goes to stderr to keep stdout clean for potential piping.
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*" >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"  >&2; }

# usage
usage() {
    cat >&2 <<EOF

SEXmer-scan.sh - K-mer sex-specificity classifier from KMC or Jellyfish dump files.

Usage: SEXmer-scan.sh -m <male_files> -f <female_files> [OPTIONS]

Mandatory:
  -m, --male           Comma-separated list of male dump files
  -f, --female         Comma-separated list of female dump files

Optional:
  --prefix             Output filename prefix                          [default: output]
  --mem                Memory budget for sort operations (e.g. 8G)     [default: ${MEM}]
  -t, --threads        CPU threads                                     [default: ${THREADS}]
  --neutral-max        Maximum neutral k-mers to retain, 0=keep all    [default: ${NEUTRAL_MAX}]
  --tmpdir             Parent directory for the temporary work folder  [default: current dir]
  --min-count          Minimum k-mer count to retain (all categories)  [default: ${MIN_COUNT}]
  --max-count          Maximum k-mer count to retain (all categories)  [default: ${MAX_COUNT}]
  --fold-threshold     Fold-change cutoff for MBK/FBK                  [default: ${FOLD_THRESHOLD}]
  --seed               Random seed for neutral k-mer sub-sampling      [default: ${SEED}]
  -h, --help           Show this help and exit

Categories in output TSV:
  MSK      Consistent across all males, pooled count >= min-count, absent in all females
  FSK      Consistent across all females, pooled count >= min-count, absent in all males
  MBK      Consistent across all males, male count >= fold-threshold x female count (NOT MSK)
  FBK      Consistent across all females, female count >= fold-threshold x male count (NOT FSK)
  neutral  Present in both sexes (mc > 0 AND fc > 0) and within fold-threshold in both
           directions (mc < fold-threshold x fc AND fc < fold-threshold x mc); not MSK/FSK/MBK/FBK

EOF
    exit 1
}

[[ $# -eq 0 ]] && usage

# argument parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--male)           MALE_INPUT="$2";    shift 2 ;;
        -f|--female)         FEMALE_INPUT="$2";  shift 2 ;;
        --prefix)            PREFIX="$2";        shift 2 ;;
        --mem)               MEM="$2";           shift 2 ;;
        -t|--threads)        THREADS="$2";       shift 2 ;;
        --neutral-max)       NEUTRAL_MAX="$2";   shift 2 ;;
        --tmpdir)            TMPDIR_BASE="$2";   shift 2 ;;
        --min-count)         MIN_COUNT="$2";     shift 2 ;;
        --max-count)         MAX_COUNT="$2";     shift 2 ;;
        --fold-threshold)    FOLD_THRESHOLD="$2"; shift 2 ;;
        --seed)              SEED="$2";          shift 2 ;;
        -h|--help)           usage ;;
        *) error "Unknown option '$1'"; usage ;;
    esac
done

# validate mandatory arguments
[[ -z "$MALE_INPUT" ]]   && { error "-m/--male files not specified.";   usage; }
[[ -z "$FEMALE_INPUT" ]] && { error "-f/--female files not specified."; usage; }

# validate options
[[ "$MIN_COUNT" =~ ^[1-9][0-9]*$ ]] || {
    error "--min-count must be a positive integer."; exit 1; }

[[ "$MAX_COUNT" =~ ^[1-9][0-9]*$ ]] || {
    error "--max-count must be a positive integer."; exit 1; }

[[ $(awk -v mn="$MIN_COUNT" -v mx="$MAX_COUNT" 'BEGIN{print (mn < mx)}') -eq 1 ]] || {
    error "--min-count must be less than --max-count."; exit 1; }

[[ "$FOLD_THRESHOLD" =~ ^[0-9]+(\.[0-9]+)?$ ]] && \
[[ $(awk -v v="$FOLD_THRESHOLD" 'BEGIN{print (v>0)}') -eq 1 ]] || {
    error "--fold-threshold must be a positive number."; exit 1; }

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    error "--threads must be a positive integer."; exit 1; }

[[ "$SEED" =~ ^[0-9]+$ ]] || {
    error "--seed must be a non-negative integer."; exit 1; }

# validate directories
[[ -d "$TMPDIR_BASE" ]] || { error "Temporary parent directory does not exist: $TMPDIR_BASE"; exit 1; }
[[ -w "$TMPDIR_BASE" ]] || { error "Temporary parent directory is not writable: $TMPDIR_BASE"; exit 1; }

# validate input files
IFS=',' read -ra MALE_FILES   <<< "$MALE_INPUT"
IFS=',' read -ra FEMALE_FILES <<< "$FEMALE_INPUT"

for f in "${MALE_FILES[@]}" "${FEMALE_FILES[@]}"; do
    [[ -r "$f" ]] || { error "Cannot read input file: $f"; exit 1; }
done

# check GNU sort
if ! sort --version 2>&1 | grep -q 'GNU'; then
    error "GNU sort is required but not found. On macOS: brew install coreutils"
    exit 1
fi

# set up temp directory
# TMPDIR_BASE : user-specified parent directory (--tmpdir)
# SCAN_TMPDIR : actual working directory created for this run
SCAN_TMPDIR="${TMPDIR_BASE}/sexmer_scan_tmp_$$"
mkdir -p "$SCAN_TMPDIR"
cleanup() { rm -rf "$SCAN_TMPDIR"; }
trap cleanup EXIT

# log run parameters
info "SEXmer-scan starting"
info "Parameters: min-count=${MIN_COUNT}, max-count=${MAX_COUNT}, fold-threshold=${FOLD_THRESHOLD}, neutral-max=${NEUTRAL_MAX}, seed=${SEED}"
info "Settings  : mem=${MEM}, threads=${THREADS}"
info "Temp dir  : ${SCAN_TMPDIR}"
info "Male files  (${#MALE_FILES[@]}): ${MALE_FILES[*]}"
info "Female files (${#FEMALE_FILES[@]}): ${FEMALE_FILES[*]}"

# divide_mem
# Split a memory string (e.g. 4G, 2048M) evenly across N parallel sort processes.
# Returns megabytes with a minimum floor of 256M.
divide_mem() {
    local mem="$1" n="$2"
    local num unit mb result
    if [[ "$mem" =~ ^([0-9]+)([GgMmKk]?)$ ]]; then
        num="${BASH_REMATCH[1]}"
        unit="${BASH_REMATCH[2]^^}"
        case "$unit" in
            G) mb=$(( num * 1024 )) ;;
            M) mb="$num" ;;
            K) mb=$(( (num + 1023) / 1024 )) ;;
            *) mb=1 ;;
        esac
        result=$(( mb / n ))
        [[ "$result" -lt 256 ]] && result=256
        echo "${result}M"
    else
        echo "$mem"
    fi
}

# open_file
# Print the contents of a dump file to stdout, transparently decompressing
# .gz files via gzip. No extra installation needed; gzip is standard on all
# Linux and macOS systems.
open_file() {
    local f="$1"
    if [[ "$f" == *.gz ]]; then
        gzip -dc "$f"
    else
        cat "$f"
    fi
}

# aggregate_counts
# Sort each input file independently in parallel, then merge with sort -m
# (O(k·log n) merge vs O(N·k·log(N·k)) full re-sort). Counts for duplicate
# k-mers are summed in a streaming awk pass.
# Output: sorted "KMER POOLED_COUNT" file.
aggregate_counts() {
    local label="$1" outfile="$2"
    shift 2
    local files=("$@")
    local n_files=${#files[@]}

    info "Aggregating ${label} counts from ${n_files} file(s)..."

    local per_sort_mem per_sort_threads
    per_sort_mem=$(divide_mem "$MEM" "$n_files")
    per_sort_threads=$(( THREADS / n_files ))
    [[ "$per_sort_threads" -lt 1 ]] && per_sort_threads=1

    local i
    for i in "${!files[@]}"; do
        open_file "${files[$i]}" \
        | sort -k1,1 -S "$per_sort_mem" --parallel="$per_sort_threads" \
            -T "$SCAN_TMPDIR" \
            > "$SCAN_TMPDIR/${label}_agg_sorted_${i}.txt" &
    done
    wait

    sort -m -k1,1 -S "$MEM" --parallel="$THREADS" -T "$SCAN_TMPDIR" \
        "$SCAN_TMPDIR"/${label}_agg_sorted_*.txt \
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
    ' > "$outfile"

    info "  ${label}: $(wc -l < "$outfile") unique k-mers after aggregation."
}

# build_consistency_index
# Identify k-mers present in every file of a sex group. Per-file sorts run in
# parallel; MEM and THREADS are divided across them to stay within budget.
# Pre-sorted files are merged with sort -m, then uniq -c selects k-mers whose
# occurrence count equals the total number of files.
# Output: sorted list of k-mers passing the consistency requirement.
build_consistency_index() {
    local label="$1" outfile="$2" n_files="$3"
    shift 3
    local files=("$@")

    info "Building ${label} consistency index (must appear in all ${n_files} file(s))..."

    local per_sort_mem per_sort_threads
    per_sort_mem=$(divide_mem "$MEM" "$n_files")
    per_sort_threads=$(( THREADS / n_files ))
    [[ "$per_sort_threads" -lt 1 ]] && per_sort_threads=1

    local i
    for i in "${!files[@]}"; do
        open_file "${files[$i]}" \
        | awk '$2 > 0 { print $1 }' \
        | sort -S "$per_sort_mem" --parallel="$per_sort_threads" -T "$SCAN_TMPDIR" \
        > "$SCAN_TMPDIR/${label}_presence_${i}.txt" &
    done
    wait

    sort -m -S "$MEM" --parallel="$THREADS" -T "$SCAN_TMPDIR" \
        "$SCAN_TMPDIR"/${label}_presence_*.txt \
    | uniq -c \
    | awk -v n="$n_files" '$1 == n { print $2 }' \
    > "$outfile"

    info "  ${label}: $(wc -l < "$outfile") k-mers consistent across all ${n_files} file(s)."
}

# STEP 1: Aggregate pooled counts per sex
aggregate_counts "male"   "$SCAN_TMPDIR/male_agg.txt"   "${MALE_FILES[@]}"
aggregate_counts "female" "$SCAN_TMPDIR/female_agg.txt" "${FEMALE_FILES[@]}"

# STEP 2: Build per-sex consistency indexes
N_MALE=${#MALE_FILES[@]}
N_FEMALE=${#FEMALE_FILES[@]}

build_consistency_index "male"   "$SCAN_TMPDIR/male_consistent.txt"   "$N_MALE"   "${MALE_FILES[@]}"
build_consistency_index "female" "$SCAN_TMPDIR/female_consistent.txt" "$N_FEMALE" "${FEMALE_FILES[@]}"

# STEP 3: Identify MSK candidates
# Consistent in all males AND pooled count >= MIN_COUNT.
info "Identifying MSK candidates (male-consistent, pooled count >= ${MIN_COUNT})..."

join "$SCAN_TMPDIR/male_consistent.txt" "$SCAN_TMPDIR/male_agg.txt" \
| awk -v mc="$MIN_COUNT" '$2 >= mc' \
> "$SCAN_TMPDIR/msk_candidates.txt"

info "  $(wc -l < "$SCAN_TMPDIR/msk_candidates.txt") MSK candidates before cross-sex filter."

# STEP 4: Identify FSK candidates
info "Identifying FSK candidates (female-consistent, pooled count >= ${MIN_COUNT})..."

join "$SCAN_TMPDIR/female_consistent.txt" "$SCAN_TMPDIR/female_agg.txt" \
| awk -v mc="$MIN_COUNT" '$2 >= mc' \
> "$SCAN_TMPDIR/fsk_candidates.txt"

info "  $(wc -l < "$SCAN_TMPDIR/fsk_candidates.txt") FSK candidates before cross-sex filter."

# STEP 5: Cross-sex exclusion
# Remove MSK candidates present in any female, and FSK candidates present in any male.
# All files are already sorted on k-mer sequence.
info "Applying cross-sex exclusion filter..."

awk '{print $1}' "$SCAN_TMPDIR/female_agg.txt" > "$SCAN_TMPDIR/female_kmers.txt"
awk '{print $1}' "$SCAN_TMPDIR/male_agg.txt"   > "$SCAN_TMPDIR/male_kmers.txt"

join -v 1 "$SCAN_TMPDIR/msk_candidates.txt" "$SCAN_TMPDIR/female_kmers.txt" > "$SCAN_TMPDIR/msk_final.txt"
join -v 1 "$SCAN_TMPDIR/fsk_candidates.txt" "$SCAN_TMPDIR/male_kmers.txt"   > "$SCAN_TMPDIR/fsk_final.txt"

MSK_COUNT=$(wc -l < "$SCAN_TMPDIR/msk_final.txt")
FSK_COUNT=$(wc -l < "$SCAN_TMPDIR/fsk_final.txt")
info "  ${MSK_COUNT} true MSK k-mers after cross-sex filter."
info "  ${FSK_COUNT} true FSK k-mers after cross-sex filter."

# STEP 6: Build the full union count table
# Field 0 (join key) is used as the kmer column so sex-specific k-mers on either
# side are never replaced by the -e fill value.
# The min/max-count filter is applied here to remove noise k-mers before
# annotation. A k-mer is retained only when at least one sex has a pooled count
# within [MIN_COUNT, MAX_COUNT]. This prevents low-count inconsistent k-mers
# from appearing as neutral in the MSK/FSK/MBK/FBK regions of the scatter plot.
info "Building full union count table (min-count=${MIN_COUNT}, max-count=${MAX_COUNT})..."

join -a 1 -a 2 -e 0 -o 0,1.2,2.2 \
    "$SCAN_TMPDIR/male_agg.txt" \
    "$SCAN_TMPDIR/female_agg.txt" \
| awk -v mn="$MIN_COUNT" -v mx="$MAX_COUNT" '
    {
        mc = ($2=="" ? 0 : $2)
        fc = ($3=="" ? 0 : $3)
        if ((mc >= mn || fc >= mn) && mc <= mx && fc <= mx)
            print $1, mc, fc
    }' \
> "$SCAN_TMPDIR/union_counts.txt"

info "  $(wc -l < "$SCAN_TMPDIR/union_counts.txt") k-mers in union table after count filtering."

awk '{print $1}' "$SCAN_TMPDIR/msk_final.txt" > "$SCAN_TMPDIR/msk_kmers.txt"
awk '{print $1}' "$SCAN_TMPDIR/fsk_final.txt" > "$SCAN_TMPDIR/fsk_kmers.txt"

# STEP 7: Pre-compute MBK and FBK k-mer lists using disk-based join.
# male_consistent.txt and female_consistent.txt can contain hundreds of millions
# of k-mers for large genomes, making it infeasible to load them into awk arrays.
# Instead, join intersects union_counts with each consistency index on disk (O(n),
# near-zero RAM), and awk applies the fold-threshold arithmetic on the result.
# The output mbk_kmers.txt and fbk_kmers.txt are small (same order as MSK/FSK),
# so they can be loaded safely as arrays in the annotation pass below.
#
# MBK: male-consistent in union, fc > 0, mc >= ft * fc, not already MSK.
# FBK: female-consistent in union, mc > 0, fc >= ft * mc, not already FSK.
# union_counts.txt fields: kmer mc fc (space-separated, sorted by kmer).
# male_consistent.txt / female_consistent.txt: single kmer column, sorted.
info "Pre-computing MBK/FBK candidates using disk-based join..."

join "$SCAN_TMPDIR/male_consistent.txt" "$SCAN_TMPDIR/union_counts.txt" \
| awk -v ft="$FOLD_THRESHOLD" '$3 > 0 && $2 >= ft * $3 { print $1 }' \
| join -v 1 - "$SCAN_TMPDIR/msk_kmers.txt" \
> "$SCAN_TMPDIR/mbk_kmers.txt"

join "$SCAN_TMPDIR/female_consistent.txt" "$SCAN_TMPDIR/union_counts.txt" \
| awk -v ft="$FOLD_THRESHOLD" '$2 > 0 && $3 >= ft * $2 { print $1 }' \
| join -v 1 - "$SCAN_TMPDIR/fsk_kmers.txt" \
> "$SCAN_TMPDIR/fbk_kmers.txt"

info "  $(wc -l < "$SCAN_TMPDIR/mbk_kmers.txt") MBK candidates."
info "  $(wc -l < "$SCAN_TMPDIR/fbk_kmers.txt") FBK candidates."

# STEP 8: Annotate categories in a single awk pass
# Only the four small k-mer sets (MSK, FSK, MBK, FBK) are loaded into arrays.
# union_counts.txt is streamed once. The neutral condition is pure arithmetic
# and requires no array lookup, identical to the original classification logic.
# Classification priority: MSK > FSK > MBK > FBK > neutral.
info "Annotating k-mers with categories..."

awk -v ft="$FOLD_THRESHOLD" \
    -v msk_f="$SCAN_TMPDIR/msk_kmers.txt" \
    -v fsk_f="$SCAN_TMPDIR/fsk_kmers.txt" \
    -v mbk_f="$SCAN_TMPDIR/mbk_kmers.txt" \
    -v fbk_f="$SCAN_TMPDIR/fbk_kmers.txt" '
    FILENAME == msk_f { msk[$1]=1; next }
    FILENAME == fsk_f { fsk[$1]=1; next }
    FILENAME == mbk_f { mbk[$1]=1; next }
    FILENAME == fbk_f { fbk[$1]=1; next }
    {
        kmer = $1
        mc   = $2
        fc   = $3

        if      (kmer in msk)                                        { cat = "MSK"     }
        else if (kmer in fsk)                                        { cat = "FSK"     }
        else if (kmer in mbk)                                        { cat = "MBK"     }
        else if (kmer in fbk)                                        { cat = "FBK"     }
        else if (mc > 0 && fc > 0 && mc < ft * fc && fc < ft * mc) { cat = "neutral" }
        else                                                         { next            }

        print kmer "\t" mc "\t" fc "\t" cat
    }
' "$SCAN_TMPDIR/msk_kmers.txt" "$SCAN_TMPDIR/fsk_kmers.txt" \
  "$SCAN_TMPDIR/mbk_kmers.txt" "$SCAN_TMPDIR/fbk_kmers.txt" \
  "$SCAN_TMPDIR/union_counts.txt" \
> "$SCAN_TMPDIR/annotated_all.txt"

info "  $(wc -l < "$SCAN_TMPDIR/annotated_all.txt") total annotated k-mers."

# STEP 9: Split categories and sub-sample neutral k-mers
# Non-neutral categories are never sub-sampled. Sub-sampling uses a streaming
# reservoir with a configurable seed for reproducibility.
info "Splitting categories and sub-sampling neutral k-mers if needed..."

# Pre-create both files so wc -l succeeds even if a category has zero k-mers.
: > "$SCAN_TMPDIR/neutral_all.txt"
: > "$SCAN_TMPDIR/non_neutral.txt"

awk -F'\t' '$4 == "neutral" { print > "'"$SCAN_TMPDIR/neutral_all.txt"'" }
            $4 != "neutral" { print > "'"$SCAN_TMPDIR/non_neutral.txt"'" }' \
    "$SCAN_TMPDIR/annotated_all.txt"

NEUTRAL_TOTAL=$(wc -l < "$SCAN_TMPDIR/neutral_all.txt")
NEUTRAL_FILE="$SCAN_TMPDIR/neutral_all.txt"
NEUTRAL_KEPT=$NEUTRAL_TOTAL
SUBSAMPLED="no"

info "  ${NEUTRAL_TOTAL} neutral k-mers before sub-sampling."

if [[ "$NEUTRAL_MAX" -gt 0 && "$NEUTRAL_TOTAL" -gt "$NEUTRAL_MAX" ]]; then
    info "  Sub-sampling neutral k-mers to ${NEUTRAL_MAX} (seed=${SEED})..."
    # Portable LCG (Numerical Recipes: a=1664525, c=1013904223, m=2^32).
    # Avoids srand()/rand() because mawk ignores the srand seed, making
    # results non-reproducible across runs. All intermediate values stay
    # below 2^53 so IEEE 754 double arithmetic is exact.
    awk -v max="$NEUTRAL_MAX" -v total="$NEUTRAL_TOTAL" -v seed="$SEED" '
    BEGIN { _a=1664525; _c=1013904223; _m=4294967296; _rng=seed+0 }
    function _rand() { _rng=(_a*_rng+_c)%_m; return _rng/_m }
    {
        remaining = total - NR + 1
        if (max > 0 && (remaining <= max || _rand() < max / remaining)) {
            print; max
        }
    }' "$SCAN_TMPDIR/neutral_all.txt" > "$SCAN_TMPDIR/neutral_sampled.txt"
    NEUTRAL_FILE="$SCAN_TMPDIR/neutral_sampled.txt"
    NEUTRAL_KEPT=$(wc -l < "$NEUTRAL_FILE")
    SUBSAMPLED="yes"
    info "  Kept ${NEUTRAL_KEPT} neutral k-mers after sub-sampling."
elif [[ "$NEUTRAL_MAX" -eq 0 ]]; then
    info "  neutral-max=0: retaining all ${NEUTRAL_TOTAL} neutral k-mers."
else
    info "  Neutral k-mer count (${NEUTRAL_TOTAL}) is within neutral-max; no sub-sampling."
fi

# STEP 10: Write main TSV output
# Concatenate non-neutral and neutral rows with header. 
# output order follows the union table order which is already sorted by k-mer.
TSV_OUT="${PREFIX}.kmers.tsv"
info "Writing main TSV output..."

{
    printf "kmer\tmale_count\tfemale_count\tcategory\n"
    cat "$SCAN_TMPDIR/non_neutral.txt" "$NEUTRAL_FILE"
} > "$TSV_OUT"

output "Main k-mer table written to: ${TSV_OUT}"

# STEP 11: Write FASTA outputs
# Extract MSK and FSK sequences directly from their final filtered files.
# msk_final.txt and fsk_final.txt contain "KMER COUNT" (space-separated),
# so count is available without re-reading the TSV.
MSK_FA="${PREFIX}.MSK.fa"
FSK_FA="${PREFIX}.FSK.fa"

info "Writing MSK FASTA output..."
awk 'BEGIN { n=0 }
     { n++; printf ">MSK_%d count=%s\n%s\n", n, $2, $1 }
' "$SCAN_TMPDIR/msk_final.txt" > "$MSK_FA"

info "Writing FSK FASTA output..."
awk 'BEGIN { n=0 }
     { n++; printf ">FSK_%d count=%s\n%s\n", n, $2, $1 }
' "$SCAN_TMPDIR/fsk_final.txt" > "$FSK_FA"

output "MSK sequences written to: ${MSK_FA}"
output "FSK sequences written to: ${FSK_FA}"

# Step 12: Tally category counts
# Single awk pass over the TSV to count per-category totals.
info "Tallying category counts..."

declare -A CAT_COUNTS
while IFS=$'\t' read -r cat cnt; do
    CAT_COUNTS["$cat"]="$cnt"
done < <(awk -F'\t' 'NR > 1 { c[$4]++ } END { for (k in c) print k "\t" c[k] }' "$TSV_OUT")

MSK_N=${CAT_COUNTS[MSK]:-0}
FSK_N=${CAT_COUNTS[FSK]:-0}
MBK_N=${CAT_COUNTS[MBK]:-0}
FBK_N=${CAT_COUNTS[FBK]:-0}


info "SEXmer scan complete."
info "  MSK     : ${MSK_N}"
info "  FSK     : ${FSK_N}"
info "  MBK     : ${MBK_N}"
info "  FBK     : ${FBK_N}"
info "  neutral : ${NEUTRAL_KEPT} (of ${NEUTRAL_TOTAL} total; sub-sampled=${SUBSAMPLED})"
