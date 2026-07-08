#!/usr/bin/env bash
# SEXmer assign - Assign sex from unknown samples using sex-specific k-mer (MSK or FSK).
# Author: Dede Kurniawan

set -euo pipefail
export LC_ALL=C

# defaults
KMER_SIZE=21
THREADS=8
TMPDIR_BASE="$(pwd)"
OUTDIR="."
OUTPUT_PREFIX="sexmer"
MARKERS=""
UNKNOWN_INPUT=""
SAMPLE_INPUT=""
MALE_INPUT=""
FEMALE_INPUT=""
TYPE=""
SIC_ENABLED=0

# log helpers
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*" >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"  >&2; }

usage() {
    cat <<EOF

SEXmer assign - Assign sex from unknown samples using sex-specific k-mer (MSK or FSK).

Usage:
  SEXmer assign <markers.fa> -i <dump_files> --type <XY|ZW> [OPTIONS]

Mandatory:
  <markers.fa>          Sex specific k-mer sequence, e.g. MSK.fa(.gz accepted)
                        For XY systems, provide MSK.
                        For ZW systems, provide FSK.
  -i, --input           K-mer dump file from unknown sample, separated by commas (.dump or .dump.gz)
  --type                Specify sex chromosome system: XY or ZW

SIC mode:
  --sic                 Enable SEXmer Iterative Classifier
  -m, --male            Known male dump files, separated by commas. Required only with --sic
  -f, --female          Known female dump files, separated by commas. Required only with --sic

Optional:
  -s, --sample          Specify each sample names, separated by commas.
                        [default: derived from dump filename by removing .dump.gz/.dump]
  -k, --kmer-size       Specify k-mer size used for marker parsing     [default: ${KMER_SIZE}]
  -t, --threads         Specify CPU threads for this task              [default: ${THREADS}]
  -o, --outdir          Specify output directory name                  [default: current dir]
  --tmpdir              Specify parent directory for the temp files    [default: current dir]
  -h, --help            Show this help message and exit

EOF
}

[[ $# -eq 0 ]] && { usage >&2; exit 1; }

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)        UNKNOWN_INPUT="$2"; shift 2 ;;
        -s|--sample)       SAMPLE_INPUT="$2";  shift 2 ;;
        -m|--male)         MALE_INPUT="$2";    shift 2 ;;
        -f|--female)       FEMALE_INPUT="$2";  shift 2 ;;
        --sic)             SIC_ENABLED=1;       shift ;;
        --type)            TYPE="$2";          shift 2 ;;
        -k|--kmer-size)    KMER_SIZE="$2";     shift 2 ;;
        -t|--threads)      THREADS="$2";       shift 2 ;;
        -o|--outdir)       OUTDIR="$2";        shift 2 ;;
        --tmpdir)          TMPDIR_BASE="$2";   shift 2 ;;
        -h|--help)         usage; exit 0 ;;
        -*) error "Unknown option '$1'"; usage >&2; exit 1 ;;
        *)  POSITIONAL+=("$1"); shift ;;
    esac
done

if [[ "$SIC_ENABLED" -eq 0 ]]; then
    [[ ${#POSITIONAL[@]} -eq 1 ]] || {
        error "Exactly one positional argument is required: <markers.fa>."; usage >&2; exit 1; }
    MARKERS="${POSITIONAL[0]}"
    [[ -r "$MARKERS" ]] || { error "Cannot read marker FASTA file: $MARKERS"; exit 1; }
else
    [[ ${#POSITIONAL[@]} -le 1 ]] || {
        error "At most one positional marker FASTA file is allowed."; usage >&2; exit 1; }
    if [[ ${#POSITIONAL[@]} -eq 1 ]]; then
        MARKERS="${POSITIONAL[0]}"
        [[ -r "$MARKERS" ]] || { error "Cannot read marker FASTA file: $MARKERS"; exit 1; }
    fi
fi

[[ -z "$UNKNOWN_INPUT" ]] && { error "-i/--input files not specified."; usage >&2; exit 1; }
[[ -z "$TYPE" ]]          && { error "--type not specified."; usage >&2; exit 1; }

if [[ "$SIC_ENABLED" -eq 1 ]]; then
    [[ -z "$MALE_INPUT" ]]   && { error "-m/--male files are required when --sic is enabled."; usage >&2; exit 1; }
    [[ -z "$FEMALE_INPUT" ]] && { error "-f/--female files are required when --sic is enabled."; usage >&2; exit 1; }
fi

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

mkdir -p "$OUTDIR" || { error "Cannot create output directory: ${OUTDIR}"; exit 1; }
[[ -d "$OUTDIR" ]] || { error "Output path is not a directory: ${OUTDIR}"; exit 1; }
[[ -w "$OUTDIR" ]] || { error "Output directory is not writable: ${OUTDIR}"; exit 1; }

command -v python3 &>/dev/null || { error "'python3' not found on PATH."; exit 1; }

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

derive_sample_name() {
    local path="$1"
    local name
    name="$(basename "$path")"
    case "$name" in
        *.dump.gz) name="${name%.dump.gz}" ;;
        *.dump)    name="${name%.dump}" ;;
        *.txt.gz)  name="${name%.txt.gz}" ;;
        *.gz)      name="${name%.gz}" ;;
    esac
    echo "$name"
}

UNKNOWN_FILES=()
UNKNOWN_SAMPLES=()
while IFS= read -r item; do UNKNOWN_FILES+=("$item"); done < <(split_csv_lines "$UNKNOWN_INPUT")

if [[ -n "$SAMPLE_INPUT" ]]; then
    while IFS= read -r item; do UNKNOWN_SAMPLES+=("$item"); done < <(split_csv_lines "$SAMPLE_INPUT")
    [[ ${#UNKNOWN_FILES[@]} -eq ${#UNKNOWN_SAMPLES[@]} ]] || {
        error "Number of --input files (${#UNKNOWN_FILES[@]}) must match number of --sample names (${#UNKNOWN_SAMPLES[@]})."; exit 1; }
else
    for f in "${UNKNOWN_FILES[@]}"; do
        UNKNOWN_SAMPLES+=("$(derive_sample_name "$f")")
    done
fi

for f in "${UNKNOWN_FILES[@]}"; do
    [[ -r "$f" ]] || { error "Cannot read dump file: $f"; exit 1; }
done

MALE_FILES=()
FEMALE_FILES=()
if [[ "$SIC_ENABLED" -eq 1 ]]; then
    while IFS= read -r item; do MALE_FILES+=("$item"); done < <(split_csv_lines "$MALE_INPUT")
    while IFS= read -r item; do FEMALE_FILES+=("$item"); done < <(split_csv_lines "$FEMALE_INPUT")
    for f in "${MALE_FILES[@]}" "${FEMALE_FILES[@]}"; do
        [[ -r "$f" ]] || { error "Cannot read dump file: $f"; exit 1; }
    done
fi

if [[ ${#UNKNOWN_SAMPLES[@]} -gt 1 ]]; then
    DUP_SAMPLE=$(printf '%s\n' "${UNKNOWN_SAMPLES[@]}" | sort | uniq -d | head -n 1 || true)
    if [[ -n "$DUP_SAMPLE" ]]; then
        warn "Duplicate sample name detected: ${DUP_SAMPLE}. Consider providing unique names with -s/--sample."
    fi
fi

ASSIGN_TMPDIR="${TMPDIR_BASE}/sexmer_assign_tmp_$$"
mkdir -p "$ASSIGN_TMPDIR"
cleanup() { rm -rf "$ASSIGN_TMPDIR"; }
trap cleanup EXIT

SAMPLE_TABLE="${ASSIGN_TMPDIR}/samples.tsv"
MALE_TABLE="${ASSIGN_TMPDIR}/known_male.tsv"
FEMALE_TABLE="${ASSIGN_TMPDIR}/known_female.tsv"
: > "$SAMPLE_TABLE"
: > "$MALE_TABLE"
: > "$FEMALE_TABLE"

for i in "${!UNKNOWN_FILES[@]}"; do
    printf '%s\t%s\n' "${UNKNOWN_SAMPLES[$i]}" "${UNKNOWN_FILES[$i]}" >> "$SAMPLE_TABLE"
done

if [[ "$SIC_ENABLED" -eq 1 ]]; then
    for f in "${MALE_FILES[@]}"; do
        printf '%s\t%s\n' "$(derive_sample_name "$f")" "$f" >> "$MALE_TABLE"
    done
    for f in "${FEMALE_FILES[@]}"; do
        printf '%s\t%s\n' "$(derive_sample_name "$f")" "$f" >> "$FEMALE_TABLE"
    done
fi

REPORT_OUT="${OUTDIR%/}/${OUTPUT_PREFIX}.assign.txt"
SIC_REPORT_OUT="${OUTDIR%/}/${OUTPUT_PREFIX}.SIC.report.txt"
ENGINE="${ASSIGN_TMPDIR}/sexmer_assign_engine.py"

info "SEXmer assign starting"
info "Parameters: kmer-size=${KMER_SIZE}, type=${TYPE}"
info "Settings  : threads=${THREADS}"
if [[ "$SIC_ENABLED" -eq 1 ]]; then
    info "SIC mode       : enabled"
    if [[ -n "$MARKERS" ]]; then
        info "Initial marker : ${MARKERS}"
    else
        info "Initial marker : generated from known samples"
    fi
    info "Known male     : ${#MALE_FILES[@]} sample(s)"
    info "Known female   : ${#FEMALE_FILES[@]} sample(s)"
else
    info "Marker file    : ${MARKERS}"
fi
if [[ -n "$SAMPLE_INPUT" ]]; then
    info "Sample names   : provided by -s/--sample"
else
    info "Sample names   : derived from dump filenames"
fi
info "Input samples  : ${#UNKNOWN_FILES[@]}"
info "Temp dir       : ${ASSIGN_TMPDIR}"
info "Output         : ${REPORT_OUT}"
[[ "$SIC_ENABLED" -eq 1 ]] && info "SIC report     : ${SIC_REPORT_OUT}"

cat > "$ENGINE" <<'PY'
#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import math
import multiprocessing as mp
import os
import re
import subprocess
import sys
from typing import Dict, Iterator, List, Optional, Sequence, Set, Tuple

BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}
_G_MARKERS: Optional[Set[int]] = None
MIN_COUNT = 3
MAX_COUNT = 1000
SIC_MIN_REF = 8
SIC_TARGET_REF = 10
SIC_MAX_ITER = 10
SIC_STRONG_GAP = 50.0
SIC_PATIENCE = 2

def log(kind: str, msg: str) -> None:
    print(f"[{kind}] {msg}", file=sys.stderr)

def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def parse_fasta(path: str) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0] or "sequence"
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks)

def encode_kmer(seq: str, k: int) -> Optional[int]:
    if len(seq) != k:
        return None
    code = 0
    base_code_get = BASE_CODE.get
    for base in seq:
        val = base_code_get(base)
        if val is None:
            return None
        code = (code << 2) | val
    return code

def iter_encoded_kmers_from_sequence(seq: str, k: int) -> Iterator[int]:
    if len(seq) < k:
        return
    mask = (1 << (2 * k)) - 1
    code = 0
    run = 0
    base_code_get = BASE_CODE.get
    for base in seq:
        val = base_code_get(base)
        if val is None:
            code = 0
            run = 0
            continue
        code = ((code << 2) | val) & mask
        run += 1
        if run >= k:
            yield code

def load_marker_kmers(path: str, k: int) -> Tuple[Set[int], Dict[str, int]]:
    kmers: Set[int] = set()
    stats = {
        "records": 0,
        "short_records": 0,
        "exact_kmer_records": 0,
        "long_records": 0,
        "bad_exact_records": 0,
        "raw_kmers": 0,
        "unique_kmers": 0,
    }

    for _name, seq in parse_fasta(path):
        stats["records"] += 1
        seq = re.sub(r"\s+", "", seq)
        if len(seq) < k:
            stats["short_records"] += 1
            continue
        if len(seq) == k:
            stats["exact_kmer_records"] += 1
            code = encode_kmer(seq, k)
            if code is not None:
                kmers.add(code)
                stats["raw_kmers"] += 1
            else:
                stats["bad_exact_records"] += 1
            continue
        stats["long_records"] += 1
        for code in iter_encoded_kmers_from_sequence(seq, k):
            kmers.add(code)
            stats["raw_kmers"] += 1

    stats["unique_kmers"] = len(kmers)
    if not kmers:
        raise ValueError(f"No valid marker k-mers were generated from: {path}")
    return kmers, stats

def load_marker_kmers_if_any(path: str, k: int) -> Tuple[Set[int], Dict[str, int]]:
    try:
        return load_marker_kmers(path, k)
    except ValueError as exc:
        if str(exc).startswith("No valid marker k-mers were generated"):
            return set(), {"records": 0, "short_records": 0, "exact_kmer_records": 0, "long_records": 0, "bad_exact_records": 0, "raw_kmers": 0, "unique_kmers": 0}
        raise

def read_sample_table(path: str) -> List[Tuple[int, str, str]]:
    rows: List[Tuple[int, str, str]] = []
    with open(path, "rt") as fh:
        for idx, raw in enumerate(fh):
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"Malformed sample table line {idx + 1}: {line}")
            sample, dump_path = parts
            rows.append((idx, sample, dump_path))
    return rows

def infer_marker_label(path: str, type_: str) -> Tuple[str, str]:
    base = os.path.basename(path).lower()
    if "msk" in base:
        return "MSK", "filename"
    if "fsk" in base:
        return "FSK", "filename"
    if type_ == "XY":
        return "MSK", "type"
    if type_ == "ZW":
        return "FSK", "type"
    return "marker", "unknown"

def iter_dump(path: str, k: int) -> Iterator[Tuple[int, int]]:
    with open_text(path) as fh:
        for raw in fh:
            if not raw:
                continue
            parts = raw.split()
            if not parts:
                continue
            code = encode_kmer(parts[0], k)
            if code is None:
                continue
            count = 1
            if len(parts) > 1:
                try:
                    count = int(float(parts[1]))
                except Exception:
                    count = 1
            if count > 0:
                yield code, count

def init_worker(markers: Set[int]) -> None:
    global _G_MARKERS
    _G_MARKERS = markers

def scan_dump_worker(item: Tuple[int, str, str, int]) -> Dict[str, object]:
    assert _G_MARKERS is not None
    order, sample, dump_path, k = item
    markers = _G_MARKERS
    sample_total_kmers = 0
    marker_detected = 0

    with open_text(dump_path) as fh:
        for raw in fh:
            if not raw:
                continue
            parts = raw.split()
            if not parts:
                continue
            code = encode_kmer(parts[0], k)
            if code is None:
                continue
            sample_total_kmers += 1
            if code in markers:
                marker_detected += 1

    return {
        "order": order,
        "sample": sample,
        "dump": dump_path,
        "sample_total_kmers": sample_total_kmers,
        "marker_detected": marker_detected,
    }

def scan_samples(samples: List[Tuple[int, str, str]], markers: Set[int], k: int, threads: int) -> List[Dict[str, object]]:
    items = [(order, sample, dump_path, k) for order, sample, dump_path in samples]
    workers = min(threads, len(items)) if items else 1
    if workers < 1:
        workers = 1
    if workers == 1:
        init_worker(markers)
        results = [scan_dump_worker(item) for item in items]
    else:
        ctx = mp.get_context("fork") if hasattr(os, "fork") else mp.get_context("spawn")
        with ctx.Pool(processes=workers, initializer=init_worker, initargs=(markers,)) as pool:
            results = list(pool.imap(scan_dump_worker, items, chunksize=1))
    marker_total = len(markers)
    for r in results:
        r["marker_ratio"] = pct(int(r["marker_detected"]), marker_total)
    return results

def write_scan_lite_script(path: str) -> None:
    script = r'''#!/usr/bin/env bash
set -euo pipefail
TMP="$1"
THREADS="$2"
K="$3"
MIN_COUNT="$4"
MAX_COUNT="$5"
shift 5
male_files=()
female_files=()
mode="male"
for f in "$@"; do
    if [[ "$f" == "--" ]]; then
        mode="female"
        continue
    fi
    if [[ "$mode" == "male" ]]; then
        male_files+=("$f")
    else
        female_files+=("$f")
    fi
done
open_file() {
    local f="$1"
    if [[ "$f" == *.gz ]]; then
        gzip -dc "$f"
    else
        cat "$f"
    fi
}
aggregate_counts() {
    local label="$1" outfile="$2"
    shift 2
    local files=("$@")
    local i
    for i in "${!files[@]}"; do
        open_file "${files[$i]}" \
        | awk -v k="$K" '$2 > 0 && length($1) == k && $1 ~ /^[ACGTacgt]+$/ { print toupper($1), $2 }' \
        | sort -k1,1 --parallel="$THREADS" -T "$TMP" \
        > "$TMP/${label}_agg_sorted_${i}.txt" &
    done
    wait
    sort -m -k1,1 --parallel="$THREADS" -T "$TMP" "$TMP"/${label}_agg_sorted_*.txt \
    | awk 'BEGIN { prev=""; sum=0 } { if ($1 == prev) { sum += $2 } else { if (prev != "") print prev, sum; prev = $1; sum = $2 } } END { if (prev != "") print prev, sum }' \
    > "$outfile"
}
build_consistency_index() {
    local label="$1" outfile="$2" n_files="$3"
    shift 3
    local files=("$@")
    local i
    for i in "${!files[@]}"; do
        open_file "${files[$i]}" \
        | awk -v k="$K" '$2 > 0 && length($1) == k && $1 ~ /^[ACGTacgt]+$/ { print toupper($1) }' \
        | sort -u --parallel="$THREADS" -T "$TMP" \
        > "$TMP/${label}_presence_${i}.txt" &
    done
    wait
    sort -m --parallel="$THREADS" -T "$TMP" "$TMP"/${label}_presence_*.txt \
    | uniq -c \
    | awk -v n="$n_files" '$1 == n { print $2 }' \
    > "$outfile"
}
aggregate_counts "male" "$TMP/male_agg.txt" "${male_files[@]}"
aggregate_counts "female" "$TMP/female_agg.txt" "${female_files[@]}"
build_consistency_index "male" "$TMP/male_consistent.txt" "${#male_files[@]}" "${male_files[@]}"
build_consistency_index "female" "$TMP/female_consistent.txt" "${#female_files[@]}" "${female_files[@]}"
join "$TMP/male_consistent.txt" "$TMP/male_agg.txt" | awk -v mn="$MIN_COUNT" -v mx="$MAX_COUNT" '$2 >= mn && $2 <= mx' > "$TMP/msk_candidates.txt"
join "$TMP/female_consistent.txt" "$TMP/female_agg.txt" | awk -v mn="$MIN_COUNT" -v mx="$MAX_COUNT" '$2 >= mn && $2 <= mx' > "$TMP/fsk_candidates.txt"
awk '{print $1}' "$TMP/female_agg.txt" > "$TMP/female_kmers.txt"
awk '{print $1}' "$TMP/male_agg.txt" > "$TMP/male_kmers.txt"
join -v 1 "$TMP/msk_candidates.txt" "$TMP/female_kmers.txt" > "$TMP/msk_final.txt"
join -v 1 "$TMP/fsk_candidates.txt" "$TMP/male_kmers.txt" > "$TMP/fsk_final.txt"
awk '{n++; printf ">MSK_%d count=%s\n%s\n", n, $2, $1}' "$TMP/msk_final.txt" > "$TMP/MSK.fa"
awk '{n++; printf ">FSK_%d count=%s\n%s\n", n, $2, $1}' "$TMP/fsk_final.txt" > "$TMP/FSK.fa"
'''
    with open(path, "wt") as out:
        out.write(script)
    os.chmod(path, 0o755)

def generate_sex_markers(male_samples: List[Tuple[int, str, str]], female_samples: List[Tuple[int, str, str]], k: int, tmp_parent: str, threads: int) -> Tuple[Set[int], Set[int], Dict[str, int]]:
    if not male_samples or not female_samples:
        raise ValueError("Scan-lite marker generation requires at least one male and one female reference sample")
    tmp = os.path.join(tmp_parent, f"scan_lite_{os.getpid()}_{len(male_samples)}_{len(female_samples)}")
    os.makedirs(tmp, exist_ok=True)
    script_path = os.path.join(tmp, "scan_lite.sh")
    write_scan_lite_script(script_path)
    male_paths = [path for _order, _sample, path in male_samples]
    female_paths = [path for _order, _sample, path in female_samples]
    cmd = ["bash", script_path, tmp, str(max(1, threads)), str(k), str(MIN_COUNT), str(MAX_COUNT), *male_paths, "--", *female_paths]
    try:
        proc = subprocess.run(cmd, text=True, capture_output=True)
        if proc.returncode != 0:
            msg = proc.stderr.strip() or proc.stdout.strip() or "scan-lite marker generation failed"
            raise RuntimeError(msg)
        msk, _msk_stats = load_marker_kmers_if_any(os.path.join(tmp, "MSK.fa"), k)
        fsk, _fsk_stats = load_marker_kmers_if_any(os.path.join(tmp, "FSK.fa"), k)
        stats = {
            "MSK": len(msk),
            "FSK": len(fsk),
            "male_count": len(male_samples),
            "female_count": len(female_samples),
        }
        return msk, fsk, stats
    finally:
        for root, dirs, files in os.walk(tmp, topdown=False):
            for name in files:
                try:
                    os.unlink(os.path.join(root, name))
                except OSError:
                    pass
            for name in dirs:
                try:
                    os.rmdir(os.path.join(root, name))
                except OSError:
                    pass
        try:
            os.rmdir(tmp)
        except OSError:
            pass

def build_classification_model(values: Sequence[float]) -> Dict[str, object]:
    if not values:
        raise ValueError("No marker-ratio values available for classification")

    sorted_vals = sorted(values)
    score_min = sorted_vals[0]
    score_max = sorted_vals[-1]

    if len(sorted_vals) < 2:
        return {
            "method": "fixed signal rule",
            "threshold": None,
            "largest_gap": 0.0,
            "gap_low": None,
            "gap_high": None,
            "separation": "one sample",
            "warning": "Only one sample was provided; largest-gap clustering is not possible. Fixed high/low marker signal rule was used.",
        }

    best_gap = -1.0
    best_low = sorted_vals[0]
    best_high = sorted_vals[1]
    for a, b in zip(sorted_vals, sorted_vals[1:]):
        gap = b - a
        if gap > best_gap:
            best_gap = gap
            best_low = a
            best_high = b

    threshold = (best_low + best_high) / 2.0
    has_low_signal = score_min <= 20.0
    has_high_signal = score_max >= 80.0
    if best_gap >= 20.0 and has_low_signal and has_high_signal:
        if best_gap >= 50.0:
            separation = "strong"
            warning = ""
        else:
            separation = "moderate"
            warning = "MSK/FSK-ratio separation is moderate; inspect the result manually."
        return {
            "method": "largest gap clustering",
            "threshold": threshold,
            "largest_gap": best_gap,
            "gap_low": best_low,
            "gap_high": best_high,
            "separation": separation,
            "warning": warning,
        }

    if score_min >= 80.0:
        separation = "single high-signal group"
        warning = "All samples have high MSK/FSK signal; fixed signal rule was used instead of forcing two clusters."
    elif score_max <= 20.0:
        separation = "single low-signal group"
        warning = "All samples have low MSK/FSK signal; fixed signal rule was used instead of forcing two clusters."
    else:
        separation = "no clear separation"
        warning = "MSK/FSK-ratio separation is weak or incomplete; fixed signal rule was used and intermediate samples are ambiguous."

    return {
        "method": "fixed signal rule",
        "threshold": None,
        "largest_gap": best_gap,
        "gap_low": best_low,
        "gap_high": best_high,
        "separation": separation,
        "warning": warning,
    }

def high_signal_sex(type_: str) -> str:
    return "male" if type_ == "XY" else "female"

def low_signal_sex(type_: str) -> str:
    return "female" if type_ == "XY" else "male"

def assign_sex(score: float, type_: str, model: Dict[str, object]) -> str:
    method = str(model["method"])
    if method == "largest gap clustering":
        threshold = float(model["threshold"])
        return high_signal_sex(type_) if score >= threshold else low_signal_sex(type_)
    if score >= 80.0:
        return high_signal_sex(type_)
    if score <= 20.0:
        return low_signal_sex(type_)
    return "ambiguous"

def confidence_for_ratio(score: float, model: Dict[str, object]) -> str:
    method = str(model["method"])
    separation = str(model["separation"])
    if method == "largest gap clustering":
        if separation == "strong":
            return "high"
        if separation == "moderate":
            return "medium"
        return "low"
    if 20.0 < score < 80.0:
        return "low"
    if separation == "one sample":
        return "medium"
    if separation in {"single high-signal group", "single low-signal group"}:
        return "high"
    return "medium"

def pct(num: int, den: int) -> float:
    if den <= 0:
        return 0.0
    return (num / den) * 100.0

def fmt_float(value: object, digits: int = 4, percent: bool = False) -> str:
    if value is None:
        return "NA"
    try:
        f = float(value)
    except Exception:
        return str(value)
    if math.isnan(f) or math.isinf(f):
        return "NA"
    text = f"{f:.{digits}f}"
    return f"{text}%" if percent else text

def marker_names(marker_label: str) -> Tuple[str, str, str, str, str]:
    if marker_label in {"MSK", "FSK"}:
        return (
            f"{marker_label}_detected",
            f"{marker_label}_total",
            f"{marker_label}_ratio(%)",
            f"{marker_label}-ratio",
            f"{marker_label} ratio (%)",
        )
    return "sex_marker_detected", "sex_marker_total", "sex_marker_ratio(%)", "Sex-marker-ratio", "sex marker ratio (%)"

def write_report(out_path: str, type_: str, marker_label: str, marker_label_source: str, marker_stats: Dict[str, int], results: List[Dict[str, object]], threshold_info: Dict[str, object]) -> None:
    marker_total = int(marker_stats["unique_kmers"])
    scores = [float(r["marker_ratio"]) for r in results]
    score_min = min(scores) if scores else 0.0
    score_max = max(scores) if scores else 0.0
    detected_col, total_col, ratio_col, model_prefix, score_label = marker_names(marker_label)

    with open(out_path, "wt") as out:
        out.write("SEXmer assign result\n\n")
        out.write("Classification model\n")
        out.write(f"Method : {threshold_info['method']}\n")
        out.write(f"Score used : {score_label}\n")
        out.write(f"{model_prefix} minimum : {fmt_float(score_min, 4, True)}\n")
        out.write(f"{model_prefix} maximum : {fmt_float(score_max, 4, True)}\n")
        out.write(f"Gap lower value : {fmt_float(threshold_info['gap_low'], 4, True)}\n")
        out.write(f"Gap upper value : {fmt_float(threshold_info['gap_high'], 4, True)}\n")
        out.write(f"Largest gap : {fmt_float(threshold_info['largest_gap'], 4, True)}\n")
        out.write(f"Separation : {threshold_info['separation']}\n")
        if threshold_info.get("warning"):
            out.write(f"Warning : {threshold_info['warning']}\n")
        out.write("\n")

        header = ["sample", "assignment", "confidence", detected_col, total_col, ratio_col]
        out.write("\t".join(header) + "\n")
        for r in sorted(results, key=lambda x: int(x["order"])):
            row = [
                str(r["sample"]),
                str(r["assignment"]),
                str(r["confidence"]),
                str(r["marker_detected"]),
                str(marker_total),
                fmt_float(r["marker_ratio"], 4),
            ]
            out.write("\t".join(row) + "\n")

        out.write("\nInterpretation notes\n")
        if type_ == "XY":
            out.write(f"High {marker_label} ratio indicates male; low {marker_label} ratio indicates female.\n")
        else:
            out.write(f"High {marker_label} ratio indicates female; low {marker_label} ratio indicates male.\n")
        out.write(f"Assignment is based on {score_label} only.\n")
        out.write("Largest-gap clustering is used only when both low-ratio and high-ratio groups are clear.\n")
        out.write("Fixed signal rule: high ratio >= 80%, low ratio <= 20%, intermediate = ambiguous.\n")
        if marker_label_source == "type":
            out.write(f"Marker label was inferred from --type {type_}; use MSK naming for XY or FSK naming for ZW to make it explicit.\n")

def select_candidate(pending: List[Tuple[int, str, str]], latest: Dict[str, float], sex: str, type_: str) -> Optional[Tuple[int, str, str]]:
    if not pending:
        return None
    high_sex = high_signal_sex(type_)
    if sex == high_sex:
        ranked = sorted(pending, key=lambda x: latest.get(x[1], -1.0), reverse=True)
        best = ranked[0]
        if latest.get(best[1], 0.0) <= 20.0:
            return None
        return best
    ranked = sorted(pending, key=lambda x: latest.get(x[1], 101.0))
    best = ranked[0]
    if latest.get(best[1], 100.0) >= 80.0:
        return None
    return best

def write_sic_report(out_path: str, marker_label: str, male_original: List[Tuple[int, str, str]], female_original: List[Tuple[int, str, str]], known_final: Dict[str, float], unknown_original: List[Tuple[int, str, str]], round_records: Dict[str, List[Dict[str, str]]], round_summaries: List[Dict[str, str]], final_results: List[Dict[str, object]]) -> None:
    final_sex = {str(r["sample"]): str(r["assignment"]) for r in final_results}
    max_round = max((len(v) for v in round_records.values()), default=0)
    with open(out_path, "wt") as out:
        out.write("Known sample\n")
        out.write(f"file_name\tsample\t{marker_label}_ratio(%)\tsex\n")
        for _order, sample, path in male_original:
            out.write(f"{path}\t{sample}\t{fmt_float(known_final.get(sample), 4)}\tmale\n")
        for _order, sample, path in female_original:
            out.write(f"{path}\t{sample}\t{fmt_float(known_final.get(sample), 4)}\tfemale\n")

        out.write("\nUnknown sample\n")
        header = ["file_name", "sample"]
        for i in range(1, max_round + 1):
            header.extend([f"round{i}_{marker_label}_ratio(%)", f"round{i}_pseudo_label", f"round{i}_action"])
        header.append("final_sex")
        out.write("\t".join(header) + "\n")
        for _order, sample, path in unknown_original:
            row = [path, sample]
            records = round_records.get(sample, [])
            for i in range(max_round):
                if i < len(records):
                    rec = records[i]
                    row.extend([rec.get("ratio", "NA"), rec.get("pseudo_label", "NA"), rec.get("action", "pending")])
                else:
                    row.extend(["NA", "NA", "NA"])
            row.append(final_sex.get(sample, "NA"))
            out.write("\t".join(row) + "\n")

        out.write("\nRound summary\n")
        out.write("round\tmale_count\tfemale_count\tmarker_type\tlargest_gap(%)\tadded_sample\tpseudo_label\tstop_reason\n")
        for rec in round_summaries:
            out.write("\t".join([
                rec.get("round", "NA"),
                rec.get("male_count", "NA"),
                rec.get("female_count", "NA"),
                rec.get("marker_type", marker_label),
                rec.get("largest_gap", "NA"),
                rec.get("added_sample", "NA"),
                rec.get("pseudo_label", "NA"),
                rec.get("stop_reason", "NA"),
            ]) + "\n")

def run_standard(args: argparse.Namespace) -> int:
    marker_label, marker_label_source = infer_marker_label(args.markers, args.type)
    if marker_label_source == "type":
        log("Warning", f"Marker type was not detected from file name; using expected {marker_label} label for {args.type} mode.")

    log("Info", "Building SEXmer marker k-mer index...")
    markers, marker_stats = load_marker_kmers(args.markers, args.kmer_size)
    log("Info", f"  Marker label          : {marker_label}")
    log("Info", f"  FASTA records         : {marker_stats['records']}")
    log("Info", f"  Exact k-mer records   : {marker_stats['exact_kmer_records']}")
    log("Info", f"  Long sequence records : {marker_stats['long_records']}")
    log("Info", f"  Skipped short records : {marker_stats['short_records']}")
    log("Info", f"  Unique {marker_label} k-mers  : {len(markers)}")
    log("Info", "  Scanner backend       : 2-bit integer exact matching")

    samples = read_sample_table(args.sample_table)
    log("Info", f"Scanning {len(samples)} sample dump file(s)...")
    results = scan_samples(samples, markers, args.kmer_size, args.threads)

    threshold_info = build_classification_model([float(r["marker_ratio"]) for r in results])
    for r in results:
        score = float(r["marker_ratio"])
        r["assignment"] = assign_sex(score, args.type, threshold_info)
        r["confidence"] = confidence_for_ratio(score, threshold_info)

    if threshold_info.get("warning"):
        log("Warning", str(threshold_info["warning"]))

    write_report(args.out, args.type, marker_label, marker_label_source, marker_stats, results, threshold_info)
    log("Info", "SEXmer assignment complete.")
    log("Info", f"  Method        : {threshold_info['method']}")
    log("Info", f"  Largest gap   : {float(threshold_info['largest_gap']):.4f}%")
    log("Info", f"  Gap lower     : {fmt_float(threshold_info['gap_low'], 4, True)}")
    log("Info", f"  Gap upper     : {fmt_float(threshold_info['gap_high'], 4, True)}")
    log("Info", f"  Separation    : {threshold_info['separation']}")
    return 0


def run_sic(args: argparse.Namespace) -> int:
    marker_label = "MSK" if args.type == "XY" else "FSK"
    male_original = read_sample_table(args.male_table)
    female_original = read_sample_table(args.female_table)
    unknown_original = read_sample_table(args.sample_table)
    male_current = list(male_original)
    female_current = list(female_original)
    pending = list(unknown_original)
    round_records: Dict[str, List[Dict[str, str]]] = {sample: [] for _order, sample, _path in unknown_original}
    round_summaries: List[Dict[str, str]] = []
    best_gap = -1.0
    no_improve = 0
    final_results: List[Dict[str, object]] = []
    final_model: Dict[str, object] = {}
    final_markers: Set[int] = set()
    stop_reason = "max_iteration_reached"
    decision_markers: Set[int] = set()
    msk_count = 0
    fsk_count = 0
    markers_need_refresh = True

    log("Info", "SIC started.")
    log("Info", f"  Initial male references   : {len(male_current)}")
    log("Info", f"  Initial female references : {len(female_current)}")
    log("Info", f"  Unknown candidates        : {len(pending)}")

    for iteration in range(1, SIC_MAX_ITER + 1):
        if iteration == 1 and args.markers:
            decision_markers, loaded_stats = load_marker_kmers(args.markers, args.kmer_size)
            msk_count = len(decision_markers) if marker_label == "MSK" else 0
            fsk_count = len(decision_markers) if marker_label == "FSK" else 0
            source = "provided FASTA"
            log("Info", f"SIC round {iteration}: using provided {marker_label} marker FASTA.")
            log("Info", f"  FASTA records       : {loaded_stats['records']}")
        elif markers_need_refresh:
            msk, fsk, marker_stats = generate_sex_markers(male_current, female_current, args.kmer_size, os.path.dirname(args.out) or ".", args.threads)
            decision_markers = msk if marker_label == "MSK" else fsk
            msk_count = int(marker_stats["MSK"])
            fsk_count = int(marker_stats["FSK"])
            source = "scan-lite regenerated"
            log("Info", f"SIC round {iteration}: generated markers with shell scan-lite.")
        else:
            source = "previous round markers"
            log("Info", f"SIC round {iteration}: reusing previous markers.")

        if not decision_markers:
            raise ValueError(f"No {marker_label} markers are available in SIC round {iteration}.")

        log("Info", f"  Marker source       : {source}")
        log("Info", f"  Male references     : {len(male_current)}")
        log("Info", f"  Female references   : {len(female_current)}")
        log("Info", f"  MSK k-mers          : {msk_count}")
        log("Info", f"  FSK k-mers          : {fsk_count}")
        log("Info", f"  Decision marker     : {marker_label}")
        log("Info", f"  Remaining unknowns  : {len(pending)}")
        log("Info", "  Scanner backend     : 2-bit integer exact matching")

        if not pending:
            stop_reason = "no_remaining_unknown"
            final_markers = decision_markers
            final_results = scan_samples(unknown_original, final_markers, args.kmer_size, args.threads)
            final_model = build_classification_model([float(r["marker_ratio"]) for r in final_results])
            round_summaries.append({
                "round": str(iteration),
                "male_count": str(len(male_current)),
                "female_count": str(len(female_current)),
                "marker_type": marker_label,
                "largest_gap": fmt_float(final_model["largest_gap"], 4),
                "added_sample": "NA",
                "pseudo_label": "NA",
                "stop_reason": stop_reason,
            })
            break

        scan_targets = pending
        scanned_results = scan_samples(scan_targets, decision_markers, args.kmer_size, args.threads)
        latest = {str(r["sample"]): float(r["marker_ratio"]) for r in scanned_results}
        model_targets = scanned_results if len(scanned_results) >= 2 else scan_samples(unknown_original, decision_markers, args.kmer_size, args.threads)
        model = build_classification_model([float(r["marker_ratio"]) for r in model_targets])
        largest_gap = float(model["largest_gap"])
        if largest_gap > best_gap:
            best_gap = largest_gap
            no_improve = 0
        else:
            no_improve += 1

        log("Info", f"  Largest gap         : {largest_gap:.4f}%")
        log("Info", f"  Separation          : {model['separation']}")

        added_samples: List[str] = []
        added_labels: List[str] = []
        scanned_names = {str(r["sample"]) for r in scanned_results}

        for r in scanned_results:
            sample = str(r["sample"])
            round_records[sample].append({
                "ratio": fmt_float(r["marker_ratio"], 4),
                "pseudo_label": "NA",
                "action": "pending",
            })
        for _order, sample, _path in unknown_original:
            if sample not in scanned_names and round_records.get(sample):
                round_records[sample].append({
                    "ratio": "NA",
                    "pseudo_label": "NA",
                    "action": "reference",
                })

        need_male = len(male_current) < SIC_MIN_REF
        need_female = len(female_current) < SIC_MIN_REF
        if not need_male and not need_female and largest_gap >= SIC_STRONG_GAP:
            stop_reason = "minimum_reference_and_strong_gap_reached"
            final_markers = decision_markers
            final_results = scan_samples(unknown_original, final_markers, args.kmer_size, args.threads)
            final_model = build_classification_model([float(r["marker_ratio"]) for r in final_results])
            round_summaries.append({
                "round": str(iteration),
                "male_count": str(len(male_current)),
                "female_count": str(len(female_current)),
                "marker_type": marker_label,
                "largest_gap": fmt_float(largest_gap, 4),
                "added_sample": "NA",
                "pseudo_label": "NA",
                "stop_reason": stop_reason,
            })
            break
        if len(male_current) >= SIC_TARGET_REF and len(female_current) >= SIC_TARGET_REF:
            stop_reason = "target_reference_reached"
            final_markers = decision_markers
            final_results = scan_samples(unknown_original, final_markers, args.kmer_size, args.threads)
            final_model = build_classification_model([float(r["marker_ratio"]) for r in final_results])
            round_summaries.append({
                "round": str(iteration),
                "male_count": str(len(male_current)),
                "female_count": str(len(female_current)),
                "marker_type": marker_label,
                "largest_gap": fmt_float(largest_gap, 4),
                "added_sample": "NA",
                "pseudo_label": "NA",
                "stop_reason": stop_reason,
            })
            break
        if no_improve >= SIC_PATIENCE:
            stop_reason = "no_largest_gap_improvement"
            final_markers = decision_markers
            final_results = scan_samples(unknown_original, final_markers, args.kmer_size, args.threads)
            final_model = build_classification_model([float(r["marker_ratio"]) for r in final_results])
            round_summaries.append({
                "round": str(iteration),
                "male_count": str(len(male_current)),
                "female_count": str(len(female_current)),
                "marker_type": marker_label,
                "largest_gap": fmt_float(largest_gap, 4),
                "added_sample": "NA",
                "pseudo_label": "NA",
                "stop_reason": stop_reason,
            })
            break

        targets: List[str] = []
        if need_male and need_female:
            targets = ["male", "female"]
        elif need_male:
            targets = ["male"]
        elif need_female:
            targets = ["female"]
        elif len(male_current) < len(female_current):
            targets = ["male"]
        elif len(female_current) < len(male_current):
            targets = ["female"]
        else:
            stop_reason = "reference_size_sufficient"
            final_markers = decision_markers
            final_results = scan_samples(unknown_original, final_markers, args.kmer_size, args.threads)
            final_model = build_classification_model([float(r["marker_ratio"]) for r in final_results])
            round_summaries.append({
                "round": str(iteration),
                "male_count": str(len(male_current)),
                "female_count": str(len(female_current)),
                "marker_type": marker_label,
                "largest_gap": fmt_float(largest_gap, 4),
                "added_sample": "NA",
                "pseudo_label": "NA",
                "stop_reason": stop_reason,
            })
            break

        for target in targets:
            candidate = select_candidate(pending, latest, target, args.type)
            if candidate is None:
                continue
            pending = [x for x in pending if x[1] != candidate[1]]
            if target == "male":
                male_current.append(candidate)
            else:
                female_current.append(candidate)
            added_samples.append(candidate[1])
            added_labels.append(target)
            if round_records.get(candidate[1]):
                round_records[candidate[1]][-1]["pseudo_label"] = target
                round_records[candidate[1]][-1]["action"] = "added_to_known"
            log("Info", f"  Added pseudo-label : {candidate[1]} -> {target} ({latest.get(candidate[1], 0.0):.4f}%)")

        if not added_samples:
            stop_reason = "no_confident_pseudo_label"
            final_markers = decision_markers
            final_results = scan_samples(unknown_original, final_markers, args.kmer_size, args.threads)
            final_model = build_classification_model([float(r["marker_ratio"]) for r in final_results])
            round_summaries.append({
                "round": str(iteration),
                "male_count": str(len(male_current)),
                "female_count": str(len(female_current)),
                "marker_type": marker_label,
                "largest_gap": fmt_float(largest_gap, 4),
                "added_sample": "NA",
                "pseudo_label": "NA",
                "stop_reason": stop_reason,
            })
            break

        round_summaries.append({
            "round": str(iteration),
            "male_count": str(len(male_current)),
            "female_count": str(len(female_current)),
            "marker_type": marker_label,
            "largest_gap": fmt_float(largest_gap, 4),
            "added_sample": ",".join(added_samples),
            "pseudo_label": ",".join(added_labels),
            "stop_reason": "continue",
        })
        markers_need_refresh = True
        final_markers = decision_markers
        final_results = scan_samples(unknown_original, final_markers, args.kmer_size, args.threads)
        final_model = build_classification_model([float(r["marker_ratio"]) for r in final_results])

    if not final_results:
        if not final_markers:
            msk, fsk, _marker_stats = generate_sex_markers(male_current, female_current, args.kmer_size, os.path.dirname(args.out) or ".", args.threads)
            final_markers = msk if marker_label == "MSK" else fsk
        final_results = scan_samples(unknown_original, final_markers, args.kmer_size, args.threads)
        final_model = build_classification_model([float(r["marker_ratio"]) for r in final_results])

    for r in final_results:
        score = float(r["marker_ratio"])
        r["assignment"] = assign_sex(score, args.type, final_model)
        r["confidence"] = confidence_for_ratio(score, final_model)

    marker_stats = {"unique_kmers": len(final_markers)}
    if final_model.get("warning"):
        log("Warning", str(final_model["warning"]))
    if stop_reason not in {"minimum_reference_and_strong_gap_reached", "target_reference_reached", "reference_size_sufficient", "no_remaining_unknown"}:
        log("Warning", f"SIC stopped with status: {stop_reason}.")

    write_report(args.out, args.type, marker_label, "SIC", marker_stats, final_results, final_model)

    known_all = [(i, s, p) for i, s, p in male_original + female_original]
    known_scan = scan_samples(known_all, final_markers, args.kmer_size, args.threads)
    known_final = {str(r["sample"]): float(r["marker_ratio"]) for r in known_scan}
    write_sic_report(args.sic_report, marker_label, male_original, female_original, known_final, unknown_original, round_records, round_summaries, final_results)

    log("Info", "SIC complete.")
    log("Info", f"  Stop reason   : {stop_reason}")
    log("Info", f"  Final {marker_label} k-mers : {len(final_markers)}")
    log("Info", f"  Largest gap   : {float(final_model['largest_gap']):.4f}%")
    log("Info", f"  Separation    : {final_model['separation']}")
    return 0

def main() -> int:
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--markers")
    ap.add_argument("--sample-table", required=True)
    ap.add_argument("--male-table")
    ap.add_argument("--female-table")
    ap.add_argument("--type", required=True, choices=["XY", "ZW"])
    ap.add_argument("--kmer-size", type=int, required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sic-report")
    ap.add_argument("--sic", action="store_true")
    args = ap.parse_args()

    try:
        if args.sic:
            if not args.male_table or not args.female_table or not args.sic_report:
                raise ValueError("SIC mode requires --male-table, --female-table, and --sic-report")
            return run_sic(args)
        if not args.markers:
            raise ValueError("Standard mode requires --markers")
        return run_standard(args)
    except Exception as exc:
        log("Error", str(exc))
        return 1

if __name__ == "__main__":
    sys.exit(main())
PY
chmod +x "$ENGINE"

PY_ARGS=(
    --sample-table "$SAMPLE_TABLE"
    --type "$TYPE"
    --kmer-size "$KMER_SIZE"
    --threads "$THREADS"
    --out "$REPORT_OUT"
)

if [[ "$SIC_ENABLED" -eq 1 ]]; then
    PY_ARGS+=(--sic --male-table "$MALE_TABLE" --female-table "$FEMALE_TABLE" --sic-report "$SIC_REPORT_OUT")
    [[ -n "$MARKERS" ]] && PY_ARGS+=(--markers "$MARKERS")
else
    PY_ARGS+=(--markers "$MARKERS")
fi

python3 "$ENGINE" "${PY_ARGS[@]}"

output "Assignment report written to: ${REPORT_OUT}"
if [[ "$SIC_ENABLED" -eq 1 ]]; then
    output "SIC report written to: ${SIC_REPORT_OUT}"
fi
info "SEXmer assign complete."