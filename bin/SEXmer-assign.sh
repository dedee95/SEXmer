#!/usr/bin/env bash
# SEXmer-assign.sh - Assign sex from unknown samples using sex-specific marker k-mers.
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
TYPE=""

# log helpers
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*" >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"  >&2; }

usage() {
    cat <<EOF

SEXmer-assign.sh - Assign sex from unknown samples using sex-specific marker k-mers.

Usage: SEXmer-assign.sh <markers.fa> -i <dump_files> --type <XY|ZW> [OPTIONS]

Mandatory:
  <markers.fa>          SEXmer marker FASTA from SEXmer-scan.
                        For XY systems, provide MSK marker FASTA.
                        For ZW systems, provide FSK marker FASTA.
  -i, --input           Comma-separated list of unknown SEXmer dump files (.dump or .dump.gz)
  --type                Sex chromosome system: XY or ZW

Optional:
  -s, --sample          Comma-separated list of sample names
                        [default: derived from dump filename by removing .dump.gz/.dump]
  -k, --kmer-size       K-mer size used for marker parsing             [default: ${KMER_SIZE}]
  -t, --threads         Number of samples processed in parallel        [default: ${THREADS}]
  -o, --outdir         Output directory                                [default: current dir]
  --tmpdir              Parent directory for temporary work folder     [default: current dir]
  -h, --help            Show this help and exit

EOF
}

[[ $# -eq 0 ]] && { usage >&2; exit 1; }

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)        UNKNOWN_INPUT="$2"; shift 2 ;;
        -s|--sample)       SAMPLE_INPUT="$2";  shift 2 ;;
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

[[ ${#POSITIONAL[@]} -eq 1 ]] || {
    error "Exactly one positional argument is required: <markers.fa>."; usage >&2; exit 1; }

MARKERS="${POSITIONAL[0]}"

[[ -z "$UNKNOWN_INPUT" ]] && { error "-i/--input files not specified."; usage >&2; exit 1; }
[[ -z "$TYPE" ]]          && { error "--type not specified."; usage >&2; exit 1; }

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

mkdir -p "$OUTDIR" || { error "Cannot create output directory: ${OUTDIR}"; exit 1; }
[[ -d "$OUTDIR" ]] || { error "Output path is not a directory: ${OUTDIR}"; exit 1; }
[[ -w "$OUTDIR" ]] || { error "Output directory is not writable: ${OUTDIR}"; exit 1; }

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
: > "$SAMPLE_TABLE"
for i in "${!UNKNOWN_FILES[@]}"; do
    printf '%s\t%s\n' "${UNKNOWN_SAMPLES[$i]}" "${UNKNOWN_FILES[$i]}" >> "$SAMPLE_TABLE"
done

REPORT_OUT="${OUTDIR%/}/${OUTPUT_PREFIX}.assign.txt"
ENGINE="${ASSIGN_TMPDIR}/sexmer_assign_engine.py"

info "SEXmer-assign starting"
info "Parameters: kmer-size=${KMER_SIZE}, type=${TYPE}"
info "Settings  : threads=${THREADS}"
info "Marker file    : ${MARKERS}"
if [[ -n "$SAMPLE_INPUT" ]]; then
    info "Sample names   : provided by -s/--sample"
else
    info "Sample names   : derived from dump filenames"
fi
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
import re
import sys
from typing import Dict, Iterator, List, Optional, Sequence, Set, Tuple

BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
_G_MARKERS: Optional[Set[int]] = None


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
        seq = re.sub(r"\s+", "", seq).upper()
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
        raise ValueError(f"No valid marker k-mers were generated from: {path}")
    return kmers, stats


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
    if not rows:
        raise ValueError("No samples found in sample table")
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


def init_worker(markers: Set[int]) -> None:
    global _G_MARKERS
    _G_MARKERS = markers


def scan_dump_worker(item: Tuple[int, str, str]) -> Dict[str, object]:
    assert _G_MARKERS is not None
    order, sample, dump_path = item
    markers = _G_MARKERS

    sample_total_kmers = 0
    marker_detected = 0
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
            if code is None:
                malformed_lines += 1
                continue
            if code in markers:
                marker_detected += 1

    return {
        "order": order,
        "sample": sample,
        "dump": dump_path,
        "sample_total_kmers": sample_total_kmers,
        "marker_detected": marker_detected,
        "malformed_lines": malformed_lines,
    }


def build_classification_model(values: Sequence[float]) -> Dict[str, object]:
    """Build a safe classification model from MSK/FSK-ratio values.
    """
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
    if type_ == "XY":
        return "male"
    if type_ == "ZW":
        return "female"
    raise ValueError(f"Unsupported type: {type_}")


def low_signal_sex(type_: str) -> str:
    if type_ == "XY":
        return "female"
    if type_ == "ZW":
        return "male"
    raise ValueError(f"Unsupported type: {type_}")


def assign_sex(score: float, type_: str, model: Dict[str, object]) -> str:
    method = str(model["method"])
    if method == "largest gap clustering":
        threshold = float(model["threshold"])
        return high_signal_sex(type_) if score >= threshold else low_signal_sex(type_)

    # Fixed signal rule. Do not force classification for intermediate signal.
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

    # Fixed signal rule confidence.
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


def write_report(
    out_path: str,
    type_: str,
    marker_label: str,
    marker_label_source: str,
    marker_stats: Dict[str, int],
    results: List[Dict[str, object]],
    threshold_info: Dict[str, object],
) -> None:
    marker_total = int(marker_stats["unique_kmers"])
    scores = [float(r["marker_ratio"]) for r in results]
    score_min = min(scores) if scores else 0.0
    score_max = max(scores) if scores else 0.0

    if marker_label in {"MSK", "FSK"}:
        detected_col = f"{marker_label}_detected"
        total_col = f"{marker_label}_total"
        ratio_col = f"{marker_label}_ratio(%)"
        model_prefix = f"{marker_label}-ratio"
        score_label = f"{marker_label} ratio (%)"
    else:
        detected_col = "sex_marker_detected"
        total_col = "sex_marker_total"
        ratio_col = "sex_marker_ratio(%)"
        model_prefix = "Sex-marker-ratio"
        score_label = "sex marker ratio (%)"

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

        header = [
            "sample",
            "assignment",
            "confidence",
            detected_col,
            total_col,
            ratio_col,
        ]
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


def main() -> int:
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--sample-table", required=True)
    ap.add_argument("--type", required=True, choices=["XY", "ZW"])
    ap.add_argument("--kmer-size", type=int, required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    try:
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
        log("Info", f"  Unique {marker_label} k-mers  : {len(markers)} (2-bit encoded)")

        samples = read_sample_table(args.sample_table)
        workers = min(args.threads, len(samples))
        if workers < 1:
            workers = 1
        log("Info", f"Scanning {len(samples)} sample dump file(s) with {workers} worker(s)...")

        if workers == 1:
            init_worker(markers)
            results = [scan_dump_worker(item) for item in samples]
        else:
            ctx = mp.get_context("fork") if hasattr(os, "fork") else mp.get_context("spawn")
            with ctx.Pool(processes=workers, initializer=init_worker, initargs=(markers,)) as pool:
                # imap preserves input order; the report also sorts by original order as a safeguard.
                results = list(pool.imap(scan_dump_worker, samples, chunksize=1))

        marker_total = len(markers)
        total_malformed = 0
        for r in results:
            r["marker_ratio"] = pct(int(r["marker_detected"]), marker_total)
            total_malformed += int(r["malformed_lines"])

        if total_malformed > 0:
            log("Warning", f"Skipped {total_malformed} dump line(s) containing non-ACGT k-mers.")

        threshold_info = build_classification_model([float(r["marker_ratio"]) for r in results])
        for r in results:
            score = float(r["marker_ratio"])
            r["assignment"] = assign_sex(score, args.type, threshold_info)
            r["confidence"] = confidence_for_ratio(score, threshold_info)

        if threshold_info.get("warning"):
            log("Warning", str(threshold_info["warning"]))

        write_report(
            args.out,
            args.type,
            marker_label,
            marker_label_source,
            marker_stats,
            results,
            threshold_info,
        )

        log("Info", "SEXmer assignment complete.")
        log("Info", f"  Method        : {threshold_info['method']}")
        log("Info", f"  Largest gap   : {float(threshold_info['largest_gap']):.4f}%")
        log("Info", f"  Gap lower     : {fmt_float(threshold_info['gap_low'], 4, True)}")
        log("Info", f"  Gap upper     : {fmt_float(threshold_info['gap_high'], 4, True)}")
        log("Info", f"  Separation    : {threshold_info['separation']}")
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
    --sample-table "$SAMPLE_TABLE" \
    --type "$TYPE" \
    --kmer-size "$KMER_SIZE" \
    --threads "$THREADS" \
    --out "$REPORT_OUT"

output "Assignment report written to: ${REPORT_OUT}"
info "SEXmer-assign complete."
