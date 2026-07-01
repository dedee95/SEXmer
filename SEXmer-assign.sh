#!/usr/bin/env bash
# SEXmer-assign.sh - Assign sex from unknown samples using sex-specific marker k-mers.
# Author: Dede Kurniawan

set -euo pipefail
export LC_ALL=C

# defaults
KMER_SIZE=21
THREADS=8
TMPDIR_BASE="$(pwd)"
OUTPUT_PREFIX="sexmer"
MARKERS=""
MARKER_SEQ=""
UNKNOWN_INPUT=""
SAMPLE_INPUT=""
TYPE=""

# log helpers
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*" >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"  >&2; }

usage() {
    cat >&2 <<EOF

SEXmer-assign.sh - Assign sex from unknown samples using sex-specific marker k-mers.

Usage: SEXmer-assign.sh <markers.fa> -i <dump_files> -s <samples> --type <XY|ZW> [OPTIONS]

Mandatory:
  <markers.fa>          SEXmer marker FASTA from SEXmer-scan.
                        For XY systems, provide MSK marker FASTA.
                        For ZW systems, provide FSK marker FASTA.
  -i, --input           Comma-separated list of unknown SEXmer dump files (.dump or .dump.gz)
  -s, --sample          Comma-separated list of sample names, same order as --input
  --type                Sex chromosome system: XY or ZW

Optional:
  --marker-seq          Additional gene/marker FASTA for marker-sequence k-mer coverage evidence
  -k, --kmer-size       K-mer size used for marker parsing             [default: ${KMER_SIZE}]
  -t, --threads         Number of samples processed in parallel        [default: ${THREADS}]
  --tmpdir              Parent directory for temporary work folder     [default: current dir]
  -h, --help            Show this help and exit

EOF
    exit 1
}

[[ $# -eq 0 ]] && usage

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)        UNKNOWN_INPUT="$2"; shift 2 ;;
        -s|--sample)       SAMPLE_INPUT="$2";  shift 2 ;;
        --type)            TYPE="$2";          shift 2 ;;
        --marker-seq)      MARKER_SEQ="$2";    shift 2 ;;
        -k|--kmer-size)    KMER_SIZE="$2";    shift 2 ;;
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
[[ -z "$TYPE" ]]          && { error "--type not specified."; usage; }

[[ -r "$MARKERS" ]] || { error "Cannot read marker FASTA file: $MARKERS"; exit 1; }
if [[ -n "$MARKER_SEQ" ]]; then
    [[ -r "$MARKER_SEQ" ]] || { error "Cannot read --marker-seq FASTA file: $MARKER_SEQ"; exit 1; }
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

for bin in gzip python3; do
    command -v "$bin" &>/dev/null || { error "'${bin}' not found on PATH."; exit 1; }
done

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

UNKNOWN_FILES=()
UNKNOWN_SAMPLES=()
while IFS= read -r item; do UNKNOWN_FILES+=("$item"); done < <(split_csv_lines "$UNKNOWN_INPUT")
while IFS= read -r item; do UNKNOWN_SAMPLES+=("$item"); done < <(split_csv_lines "$SAMPLE_INPUT")

[[ ${#UNKNOWN_FILES[@]} -eq ${#UNKNOWN_SAMPLES[@]} ]] || {
    error "Number of --input files (${#UNKNOWN_FILES[@]}) must match number of --sample names (${#UNKNOWN_SAMPLES[@]})."; exit 1; }

for f in "${UNKNOWN_FILES[@]}"; do
    [[ -r "$f" ]] || { error "Cannot read dump file: $f"; exit 1; }
done

ASSIGN_TMPDIR="${TMPDIR_BASE}/sexmer_assign_tmp_$$"
mkdir -p "$ASSIGN_TMPDIR"
cleanup() { rm -rf "$ASSIGN_TMPDIR"; }
trap cleanup EXIT

SAMPLE_TABLE="${ASSIGN_TMPDIR}/samples.tsv"
: > "$SAMPLE_TABLE"
for i in "${!UNKNOWN_FILES[@]}"; do
    printf '%s\t%s\n' "${UNKNOWN_SAMPLES[$i]}" "${UNKNOWN_FILES[$i]}" >> "$SAMPLE_TABLE"
done

REPORT_OUT="${OUTPUT_PREFIX}.assign.txt"
ENGINE="${ASSIGN_TMPDIR}/sexmer_assign_engine.py"

info "SEXmer-assign starting"
info "Parameters: kmer-size=${KMER_SIZE}, type=${TYPE}"
info "Settings  : threads=${THREADS}"
info "Marker file    : ${MARKERS}"
info "Marker-seq file: ${MARKER_SEQ:-off}"
info "Input samples  : ${#UNKNOWN_FILES[@]}"
info "Temp dir       : ${ASSIGN_TMPDIR}"
info "Output         : ${REPORT_OUT}"

cat > "$ENGINE" <<'PY'
#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import math
import multiprocessing as mp
import os
import sys
from statistics import median
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
_G_MARKERS: Optional[Set[int]] = None
_G_MARKER_SEQ: Optional[Set[int]] = None


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


def encode_kmer(seq: str) -> Optional[int]:
    code = 0
    for base in seq:
        val = BASE_CODE.get(base)
        if val is None:
            return None
        code = (code << 2) | val
    return code


def iter_kmer_codes_from_sequence(seq: str, k: int) -> Iterator[int]:
    seq = seq.upper()
    if len(seq) < k:
        return
    for i in range(0, len(seq) - k + 1):
        code = encode_kmer(seq[i:i + k])
        if code is not None:
            yield code


def load_marker_kmers(path: str, k: int, label: str) -> Tuple[Set[int], Dict[str, int]]:
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

    for name, seq in parse_fasta(path):
        stats["records"] += 1
        seq = seq.replace(" ", "").replace("\t", "").upper()
        if len(seq) < k:
            stats["short_records"] += 1
            continue
        if len(seq) == k:
            stats["exact_kmer_records"] += 1
            code = encode_kmer(seq)
            if code is not None:
                kmers.add(code)
                stats["raw_kmers"] += 1
            else:
                stats["bad_exact_records"] += 1
            continue
        stats["long_records"] += 1
        for code in iter_kmer_codes_from_sequence(seq, k):
            kmers.add(code)
            stats["raw_kmers"] += 1

    stats["unique_kmers"] = len(kmers)
    if not kmers:
        raise ValueError(f"No valid {label} k-mers were generated from: {path}")
    return kmers, stats


def read_sample_table(path: str) -> List[Tuple[str, str]]:
    rows: List[Tuple[str, str]] = []
    with open(path, "rt") as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"Malformed sample table line {line_no}: {line}")
            sample, dump_path = parts
            rows.append((sample, dump_path))
    if not rows:
        raise ValueError("No samples found in sample table")
    return rows


def init_worker(markers: Set[int], marker_seq: Optional[Set[int]]) -> None:
    global _G_MARKERS, _G_MARKER_SEQ
    _G_MARKERS = markers
    _G_MARKER_SEQ = marker_seq


def scan_dump_worker(item: Tuple[str, str]) -> Dict[str, object]:
    assert _G_MARKERS is not None
    sample, dump_path = item
    markers = _G_MARKERS
    marker_seq = _G_MARKER_SEQ

    sample_total_kmers = 0
    sample_total_counts = 0
    marker_detected = 0
    marker_count_sum = 0
    marker_seq_detected = 0
    marker_seq_count_sum = 0
    malformed_lines = 0

    with open_text(dump_path) as fh:
        for raw in fh:
            if not raw:
                continue
            parts = raw.split()
            if not parts:
                continue
            kmer = parts[0].upper()
            if not kmer:
                continue
            sample_total_kmers += 1
            code = encode_kmer(kmer)
            count = 1
            if len(parts) >= 2:
                try:
                    count = int(parts[1])
                except Exception:
                    malformed_lines += 1
                    count = 1
            sample_total_counts += count
            if code is not None and code in markers:
                marker_detected += 1
                marker_count_sum += count
            if marker_seq is not None and code is not None and code in marker_seq:
                marker_seq_detected += 1
                marker_seq_count_sum += count

    return {
        "sample": sample,
        "dump": dump_path,
        "sample_total_kmers": sample_total_kmers,
        "sample_total_counts": sample_total_counts,
        "marker_detected": marker_detected,
        "marker_count_sum": marker_count_sum,
        "marker_seq_detected": marker_seq_detected,
        "marker_seq_count_sum": marker_seq_count_sum,
        "malformed_lines": malformed_lines,
    }


def largest_gap_threshold(values: Sequence[float]) -> Dict[str, object]:
    if len(values) < 2:
        return {
            "method": "fixed_50_single_sample",
            "threshold": 50.0,
            "largest_gap": 0.0,
            "gap_low": None,
            "gap_high": None,
            "separation": "not_available",
            "warning": "Only one sample was provided; clustering is not possible. Fixed 50% threshold was used.",
        }

    sorted_vals = sorted(values)
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
    if best_gap >= 20.0:
        separation = "strong"
        warning = ""
    elif best_gap >= 10.0:
        separation = "moderate"
        warning = "Marker-ratio separation is moderate; inspect the result manually."
    else:
        separation = "weak"
        warning = "Marker-ratio separation is weak; assignment confidence is low."
    return {
        "method": "largest_gap_clustering",
        "threshold": threshold,
        "largest_gap": best_gap,
        "gap_low": best_low,
        "gap_high": best_high,
        "separation": separation,
        "warning": warning,
    }


def assign_sex(ratio: float, type_: str, threshold: float) -> str:
    high = ratio >= threshold
    if type_ == "XY":
        return "male" if high else "female"
    if type_ == "ZW":
        return "female" if high else "male"
    raise ValueError(f"Unsupported type: {type_}")


def confidence_for_ratio(ratio: float, threshold_info: Dict[str, object]) -> str:
    method = str(threshold_info["method"])
    separation = str(threshold_info["separation"])
    if method == "fixed_50_single_sample":
        if ratio >= 80.0 or ratio <= 20.0:
            return "medium"
        return "low"
    if separation == "strong":
        return "high"
    if separation == "moderate":
        return "medium"
    return "low"


def pct(num: int, den: int) -> float:
    if den <= 0:
        return 0.0
    return (num / den) * 100.0


def fmt_float(value: object, digits: int = 4) -> str:
    if value is None:
        return "NA"
    try:
        f = float(value)
    except Exception:
        return str(value)
    if math.isnan(f) or math.isinf(f):
        return "NA"
    return f"{f:.{digits}f}"


def write_report(
    out_path: str,
    markers_path: str,
    marker_seq_path: Optional[str],
    k: int,
    type_: str,
    threads: int,
    marker_stats: Dict[str, int],
    marker_seq_stats: Optional[Dict[str, int]],
    results: List[Dict[str, object]],
    threshold_info: Dict[str, object],
) -> None:
    marker_total = int(marker_stats["unique_kmers"])
    marker_seq_total = int(marker_seq_stats["unique_kmers"]) if marker_seq_stats is not None else 0
    threshold = float(threshold_info["threshold"])

    ratios = [float(r["marker_ratio"]) for r in results]
    ratio_min = min(ratios) if ratios else 0.0
    ratio_max = max(ratios) if ratios else 0.0
    ratio_median = median(ratios) if ratios else 0.0

    with open(out_path, "wt") as out:
        out.write("SEXmer-assign result\n\n")
        out.write("Run information\n")
        out.write(f"Marker file        : {markers_path}\n")
        out.write(f"Marker-seq file    : {marker_seq_path if marker_seq_path else 'off'}\n")
        out.write(f"System type        : {type_}\n")
        out.write("Expected marker    : MSK for XY, FSK for ZW\n")
        out.write(f"K-mer size         : {k}\n")
        out.write(f"Threads            : {threads}\n")
        out.write(f"Samples            : {len(results)}\n\n")

        out.write("Marker statistics\n")
        out.write(f"Marker FASTA records        : {marker_stats['records']}\n")
        out.write(f"Exact k-mer records         : {marker_stats['exact_kmer_records']}\n")
        out.write(f"Long sequence records       : {marker_stats['long_records']}\n")
        out.write(f"Skipped short records       : {marker_stats['short_records']}\n")
        out.write(f"Bad exact-kmer records      : {marker_stats['bad_exact_records']}\n")
        out.write(f"Raw marker k-mers           : {marker_stats['raw_kmers']}\n")
        out.write(f"Unique marker k-mers        : {marker_total}\n\n")

        if marker_seq_stats is not None:
            out.write("Marker-seq statistics\n")
            out.write(f"Marker-seq FASTA records    : {marker_seq_stats['records']}\n")
            out.write(f"Exact k-mer records         : {marker_seq_stats['exact_kmer_records']}\n")
            out.write(f"Long sequence records       : {marker_seq_stats['long_records']}\n")
            out.write(f"Skipped short records       : {marker_seq_stats['short_records']}\n")
            out.write(f"Bad exact-kmer records      : {marker_seq_stats['bad_exact_records']}\n")
            out.write(f"Raw marker-seq k-mers       : {marker_seq_stats['raw_kmers']}\n")
            out.write(f"Unique marker-seq k-mers    : {marker_seq_total}\n\n")

        out.write("Classification model\n")
        out.write(f"Method              : {threshold_info['method']}\n")
        out.write(f"Marker-ratio minimum: {fmt_float(ratio_min, 4)}\n")
        out.write(f"Marker-ratio median : {fmt_float(ratio_median, 4)}\n")
        out.write(f"Marker-ratio maximum: {fmt_float(ratio_max, 4)}\n")
        out.write(f"Largest gap         : {fmt_float(threshold_info['largest_gap'], 4)}\n")
        out.write(f"Gap lower value     : {fmt_float(threshold_info['gap_low'], 4)}\n")
        out.write(f"Gap upper value     : {fmt_float(threshold_info['gap_high'], 4)}\n")
        out.write(f"Decision threshold  : {fmt_float(threshold, 4)}\n")
        out.write(f"Separation          : {threshold_info['separation']}\n")
        if threshold_info.get("warning"):
            out.write(f"Warning             : {threshold_info['warning']}\n")
        out.write("\n")

        out.write("Sample results\n")
        header = [
            "sample",
            "assignment",
            "confidence",
            "marker_detected",
            "marker_total",
            "marker_ratio",
            "marker_count_sum",
            "sample_total_kmers",
            "sample_total_counts",
        ]
        if marker_seq_stats is not None:
            header.extend(["marker_seq_detected", "marker_seq_total", "marker_seq_cov", "marker_seq_count_sum"])
        header.extend(["malformed_lines", "dump_file"])
        out.write("\t".join(header) + "\n")
        for r in sorted(results, key=lambda x: str(x["sample"])):
            row = [
                str(r["sample"]),
                str(r["assignment"]),
                str(r["confidence"]),
                str(r["marker_detected"]),
                str(marker_total),
                fmt_float(r["marker_ratio"], 4),
                str(r["marker_count_sum"]),
                str(r["sample_total_kmers"]),
                str(r["sample_total_counts"]),
            ]
            if marker_seq_stats is not None:
                row.extend([
                    str(r["marker_seq_detected"]),
                    str(marker_seq_total),
                    fmt_float(r["marker_seq_cov"], 4),
                    str(r["marker_seq_count_sum"]),
                ])
            row.extend([str(r["malformed_lines"]), str(r["dump"])])
            out.write("\t".join(row) + "\n")

        out.write("\nInterpretation notes\n")
        if type_ == "XY":
            out.write("For XY mode, the marker FASTA should be MSK. Samples above the threshold are assigned male.\n")
            out.write("Samples below the threshold are assigned female.\n")
        else:
            out.write("For ZW mode, the marker FASTA should be FSK. Samples above the threshold are assigned female.\n")
            out.write("Samples below the threshold are assigned male.\n")
        out.write("The optional marker-seq result is supporting evidence calculated by k-mer coverage, not read mapping.\n")


def main() -> int:
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--marker-seq", default="")
    ap.add_argument("--sample-table", required=True)
    ap.add_argument("--type", required=True, choices=["XY", "ZW"])
    ap.add_argument("--kmer-size", type=int, required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    try:
        log("Info", "Building SEXmer marker k-mer index...")
        markers, marker_stats = load_marker_kmers(args.markers, args.kmer_size, "marker")
        log("Info", f"  Unique marker k-mers: {len(markers)} (2-bit encoded)")

        marker_seq: Optional[Set[int]] = None
        marker_seq_stats: Optional[Dict[str, int]] = None
        if args.marker_seq:
            log("Info", "Building optional marker-seq k-mer index...")
            marker_seq, marker_seq_stats = load_marker_kmers(args.marker_seq, args.kmer_size, "marker-seq")
            log("Info", f"  Unique marker-seq k-mers: {len(marker_seq)} (2-bit encoded)")

        samples = read_sample_table(args.sample_table)
        workers = min(args.threads, len(samples))
        if workers < 1:
            workers = 1
        log("Info", f"Scanning {len(samples)} sample dump file(s) with {workers} worker(s)...")

        if workers == 1:
            init_worker(markers, marker_seq)
            results = [scan_dump_worker(item) for item in samples]
        else:
            ctx = mp.get_context("fork") if hasattr(os, "fork") else mp.get_context("spawn")
            with ctx.Pool(processes=workers, initializer=init_worker, initargs=(markers, marker_seq)) as pool:
                results = list(pool.imap_unordered(scan_dump_worker, samples, chunksize=1))

        marker_total = len(markers)
        marker_seq_total = len(marker_seq) if marker_seq is not None else 0
        for r in results:
            r["marker_ratio"] = pct(int(r["marker_detected"]), marker_total)
            if marker_seq is not None:
                r["marker_seq_cov"] = pct(int(r["marker_seq_detected"]), marker_seq_total)
            else:
                r["marker_seq_cov"] = 0.0

        threshold_info = largest_gap_threshold([float(r["marker_ratio"]) for r in results])
        threshold = float(threshold_info["threshold"])
        for r in results:
            ratio = float(r["marker_ratio"])
            r["assignment"] = assign_sex(ratio, args.type, threshold)
            r["confidence"] = confidence_for_ratio(ratio, threshold_info)

        if threshold_info.get("warning"):
            log("Warning", str(threshold_info["warning"]))

        write_report(
            args.out,
            args.markers,
            args.marker_seq if args.marker_seq else None,
            args.kmer_size,
            args.type,
            workers,
            marker_stats,
            marker_seq_stats,
            results,
            threshold_info,
        )

        log("Info", "SEXmer assignment complete.")
        log("Info", f"  Decision threshold: {float(threshold_info['threshold']):.4f}")
        log("Info", f"  Separation        : {threshold_info['separation']}")
        return 0
    except Exception as exc:
        log("Error", str(exc))
        return 1


if __name__ == "__main__":
    sys.exit(main())
PY
chmod +x "$ENGINE"

python3 "$ENGINE" \
    --markers "$MARKERS" \
    --marker-seq "$MARKER_SEQ" \
    --sample-table "$SAMPLE_TABLE" \
    --type "$TYPE" \
    --kmer-size "$KMER_SIZE" \
    --threads "$THREADS" \
    --out "$REPORT_OUT"

output "Assignment report written to: ${REPORT_OUT}"
info "SEXmer-assign complete."
