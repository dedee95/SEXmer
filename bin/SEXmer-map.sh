#!/usr/bin/env bash
# SEXmer map - Map specific k-mer sequence or reads into the genome reference.
# Author: Dede Kurniawan

set -euo pipefail
export LC_ALL=C

# defaults
KMER_SIZE=21
WINDOW=10000
STEP=2500
THREADS=8
PREFIX=""
OUTDIR="."
TMPDIR_BASE="$(pwd)"
GENOME=""
MARKERS=""
READS_INPUT=""
SEQ_TYPE="short"

# log helpers
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*" >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"  >&2; }

usage() {
    cat <<EOF

SEXmer map - Map specific k-mer sequence or reads into the genome reference.

Usage: SEXmer map <genome.fa> <markers.fa> --prefix <prefix> [OPTIONS]

Mandatory:
  <genome.fa>          Specify reference genome in FASTA file (.gz is accepted)
  <markers.fa>         Specify sex specific k-mer sequence, e.g. MSK.fa (.gz is accepted)
  --prefix             The prefix used on generated files

Optional:
  -k, --kmer-size      Specify k-mer size (1-63)                       [default: ${KMER_SIZE}]
  -w, --window         Specify window size in bp                       [default: ${WINDOW}]
  -s, --step           Specify sliding step size in bp                 [default: ${STEP}]
  -r, --reads          Input extracted raw reads from SEXmer reads (comma-separated)
                       Examples: -r reads.fq.gz OR -r reads_1.fq.gz,reads_2.fq.gz
  --seq-type           Specify reads type: short or long               [default: ${SEQ_TYPE}]
  -t, --threads        Specify CPU threads for this task               [default: ${THREADS}]
  -o, --outdir         Specify output directory name                   [default: current dir]
  --tmpdir             Specify parent directory for the temp files     [default: current dir]
  -h, --help           Show this help message and exit

EOF
}

[[ $# -eq 0 ]] && { usage >&2; exit 1; }

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)          PREFIX="$2";       shift 2 ;;
        -k|--kmer-size)    KMER_SIZE="$2";   shift 2 ;;
        -w|--window)       WINDOW="$2";      shift 2 ;;
        -s|--step)         STEP="$2";        shift 2 ;;
        -r|--reads)        READS_INPUT="$2"; shift 2 ;;
        --seq-type)        SEQ_TYPE="$2";    shift 2 ;;
        -t|--threads)      THREADS="$2";     shift 2 ;;
        -o|--outdir)       OUTDIR="$2";      shift 2 ;;
        --tmpdir)          TMPDIR_BASE="$2"; shift 2 ;;
        -h|--help)         usage; exit 0 ;;
        -*) error "Unknown option '$1'"; usage >&2; exit 1 ;;
        *)  POSITIONAL+=("$1"); shift ;;
    esac
done

[[ ${#POSITIONAL[@]} -eq 2 ]] || {
    error "Exactly two positional arguments are required: <genome.fa> <markers.fa>."; usage >&2; exit 1; }

GENOME="${POSITIONAL[0]}"
MARKERS="${POSITIONAL[1]}"

[[ -z "$PREFIX" ]] && { error "--prefix is required."; usage >&2; exit 1; }

[[ -r "$GENOME" ]] || { error "Cannot read genome FASTA file: $GENOME"; exit 1; }
[[ -r "$MARKERS" ]] || { error "Cannot read marker FASTA file: $MARKERS"; exit 1; }

normalize_int() {
    local value="$1"
    echo "${value//_/}"
}

KMER_SIZE="$(normalize_int "$KMER_SIZE")"
WINDOW="$(normalize_int "$WINDOW")"
STEP="$(normalize_int "$STEP")"
THREADS="$(normalize_int "$THREADS")"
[[ "$KMER_SIZE" =~ ^[1-9][0-9]*$ ]] && [[ "$KMER_SIZE" -le 63 ]] || {
    error "--kmer-size must be an integer between 1 and 63."; exit 1; }

[[ "$WINDOW" =~ ^[1-9][0-9]*$ ]] || {
    error "--window must be a positive integer."; exit 1; }

[[ "$STEP" =~ ^[1-9][0-9]*$ ]] || {
    error "--step must be a positive integer."; exit 1; }

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    error "--threads must be a positive integer."; exit 1; }

case "$SEQ_TYPE" in
    short|long) ;;
    *) error "--seq-type must be one of: short, long."; exit 1 ;;
esac

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

command -v python3 &>/dev/null || { error "python3 not found on PATH."; exit 1; }

READ_FILES=()
if [[ -n "$READS_INPUT" ]]; then
    IFS=',' read -ra READ_FILES <<< "$READS_INPUT"
    [[ ${#READ_FILES[@]} -gt 0 ]] || { error "--reads was provided but no read file was parsed."; exit 1; }
    for f in "${READ_FILES[@]}"; do
        [[ -n "$f" ]] || { error "Empty entry found in --reads comma-separated list."; exit 1; }
        [[ -r "$f" ]] || { error "Cannot read FASTQ file from --reads: $f"; exit 1; }
    done
    command -v bbmap.sh &>/dev/null || { error "bbmap.sh not found on PATH. BBTools is required when -r/--reads is used."; exit 1; }
    if [[ "$SEQ_TYPE" == "long" ]]; then
        command -v reformat.sh &>/dev/null || { error "reformat.sh not found on PATH. BBTools reformat.sh is required for --seq-type long."; exit 1; }
    fi
fi

MAP_TMPDIR="${TMPDIR_BASE}/sexmer_map_tmp_$$"
mkdir -p "$MAP_TMPDIR"
cleanup() { rm -rf "$MAP_TMPDIR"; }
trap cleanup EXIT

SCANNER="${MAP_TMPDIR}/sexmer_map_engine.py"
cat > "$SCANNER" <<'PY'
#!/usr/bin/env python3
import argparse
import gzip
import math
import multiprocessing as mp
import os
import re
import sys
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

BASE_CODE: Dict[str, int] = {
    "A": 0, "C": 1, "G": 2, "T": 3,
    "a": 0, "c": 1, "g": 2, "t": 3,
}
RC_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUME = {"M", "D", "N", "=", "X"}
QUERY_ALIGNED = {"M", "=", "X"}

_G_MARKERS: Optional[set] = None
_G_K: Optional[int] = None
_G_WINDOW: Optional[int] = None
_G_STEP: Optional[int] = None

def log(kind: str, msg: str) -> None:
    print(f"[{kind}] {msg}", file=sys.stderr)

def open_text(path: str):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def revcomp(seq: str) -> str:
    return seq.translate(RC_TRANS)[::-1]

def encode_kmer(seq: str) -> Optional[int]:
    code = 0
    for base in seq:
        val = BASE_CODE.get(base)
        if val is None:
            return None
        code = (code << 2) | val
    return code

def parse_fasta(path: str) -> Iterator[Tuple[str, str]]:
    name = None
    chunks: List[str] = []
    with open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks)

def load_genome_sizes(path: str) -> Tuple[List[str], Dict[str, int]]:
    order: List[str] = []
    sizes: Dict[str, int] = {}
    for chrom, seq in parse_fasta(path):
        if chrom in sizes:
            raise ValueError(f"Duplicate FASTA record name in genome: {chrom}")
        order.append(chrom)
        sizes[chrom] = len(seq)
    if not order:
        raise ValueError(f"No sequence found in genome FASTA: {path}")
    return order, sizes

def load_markers(path: str, k: int) -> Tuple[set, int, int]:
    raw_count = 0
    indexed = set()
    for name, seq in parse_fasta(path):
        if not seq:
            raise ValueError(f"Empty marker sequence in FASTA record: {name}")
        seq = seq.strip()
        if len(seq) != k:
            raise ValueError(f"Marker length mismatch: expected {k}, found {len(seq)} in record '{name}'")
        code = encode_kmer(seq)
        if code is None:
            raise ValueError(f"Marker contains non-ACGT base in record '{name}'")
        rc_code = encode_kmer(revcomp(seq))
        if rc_code is None:
            raise ValueError(f"Marker contains non-ACGT base in record '{name}'")
        raw_count += 1
        indexed.add(code)
        indexed.add(rc_code)
    if raw_count == 0:
        raise ValueError(f"No marker sequence found in FASTA: {path}")
    return indexed, raw_count, len(indexed)

def window_starts(seq_len: int, k: int, step: int) -> List[int]:
    if seq_len < k:
        return []
    return list(range(0, seq_len - k + 1, step))

def count_valid_sites_for_window(q_start: int, q_end: int, intervals: Sequence[Tuple[int, int]], start_hint: int) -> Tuple[int, int]:
    if q_end < q_start:
        return 0, start_hint
    n = len(intervals)
    i = start_hint
    while i < n and intervals[i][1] < q_start:
        i += 1
    new_hint = i
    total = 0
    while i < n and intervals[i][0] <= q_end:
        a, b = intervals[i]
        left = q_start if q_start > a else a
        right = q_end if q_end < b else b
        if right >= left:
            total += right - left + 1
        i += 1
    return total, new_hint

def scan_chrom_int2bit(chrom: str, seq: str, markers: set, k: int, window: int, step: int) -> Tuple[str, int, int, int, int, str]:
    seq_len = len(seq)
    starts = window_starts(seq_len, k, step)
    nwin = len(starts)
    counts = [0] * nwin
    valid_kmers = 0
    total_hits = 0
    valid_intervals: List[Tuple[int, int]] = []

    if seq_len >= k and nwin > 0:
        mask = (1 << (2 * k)) - 1
        span = window - k
        code = 0
        run = 0
        run_start = 0
        base_code_get = BASE_CODE.get
        marker_set = markers
        counts_local = counts
        step_local = step
        nwin_last = nwin - 1
        k_local = k
        mask_local = mask
        span_local = span
        starts_local = starts
        valid_intervals_append = valid_intervals.append

        for pos, base in enumerate(seq):
            val = base_code_get(base)
            if val is None:
                if run >= k_local:
                    valid_intervals_append((run_start, pos - k_local))
                code = 0
                run = 0
                run_start = pos + 1
                continue
            if run == 0:
                run_start = pos
            code = ((code << 2) | val) & mask_local
            run += 1
            if run < k_local:
                continue
            kmer_start = pos - k_local + 1
            valid_kmers += 1
            if code not in marker_set:
                continue
            total_hits += 1
            left = kmer_start - span_local
            if left <= 0:
                i_min = 0
            else:
                i_min = (left + step_local - 1) // step_local
            i_max = kmer_start // step_local
            if i_max > nwin_last:
                i_max = nwin_last
            for idx in range(i_min, i_max + 1):
                start = starts_local[idx]
                end = start + window
                if end > seq_len:
                    end = seq_len
                if kmer_start <= end - k_local:
                    counts_local[idx] += 1
        if run >= k_local:
            valid_intervals_append((run_start, seq_len - k_local))

    total_possible = max(seq_len - k + 1, 0)
    skipped_non_acgt = total_possible - valid_kmers
    rows: List[str] = []
    interval_hint = 0
    for idx, start in enumerate(starts):
        end = start + window
        if end > seq_len:
            end = seq_len
        valid_sites, interval_hint = count_valid_sites_for_window(start, end - k, valid_intervals, interval_hint)
        hits = counts[idx]
        hits_per_10kb = (hits / valid_sites * 10000.0) if valid_sites > 0 else 0.0
        rows.append(f"{chrom}\t{start}\t{end}\t{hits}\t{hits_per_10kb:.10g}\n")
    return chrom, nwin, valid_kmers, skipped_non_acgt, total_hits, "".join(rows)

def init_worker(markers: set, k: int, window: int, step: int) -> None:
    global _G_MARKERS, _G_K, _G_WINDOW, _G_STEP
    _G_MARKERS = markers
    _G_K = k
    _G_WINDOW = window
    _G_STEP = step

def scan_chrom_worker(item: Tuple[str, str]) -> Tuple[str, int, int, int, int, str]:
    chrom, seq = item
    assert _G_MARKERS is not None and _G_K is not None and _G_WINDOW is not None and _G_STEP is not None
    return scan_chrom_int2bit(chrom, seq, _G_MARKERS, _G_K, _G_WINDOW, _G_STEP)

def run_kmer_scan(genome: str, markers: set, k: int, window: int, step: int, threads: int, out_path: str) -> Tuple[int, int, int, int, int, int]:
    chrom_count = 0
    total_bases = 0
    total_windows = 0
    total_valid_kmers = 0
    total_skipped_n = 0
    total_hits = 0
    with open(out_path, "w") as out:
        out.write("chr\tstart\tend\thits\thits_per_10kb\n")
        if threads <= 1:
            for chrom, seq in parse_fasta(genome):
                chrom_count += 1
                total_bases += len(seq)
                log("Info", f"Scanning {chrom} ({len(seq)} bp)...")
                _chrom, nwin, valid, skipped, hits, rows = scan_chrom_int2bit(chrom, seq, markers, k, window, step)
                out.write(rows)
                total_windows += nwin
                total_valid_kmers += valid
                total_skipped_n += skipped
                total_hits += hits
                log("Info", f"  {chrom}: windows={nwin}, valid_kmers={valid}, skipped_N/non-ACGT={skipped}, marker_hits={hits}")
        else:
            log("Info", f"Using chromosome/contig-level multiprocessing with {threads} worker(s).")
            ctx = mp.get_context("fork") if hasattr(os, "fork") else mp.get_context("spawn")
            with ctx.Pool(processes=threads, initializer=init_worker, initargs=(markers, k, window, step)) as pool:
                for chrom, nwin, valid, skipped, hits, rows in pool.imap(scan_chrom_worker, parse_fasta(genome), chunksize=1):
                    chrom_count += 1
                    out.write(rows)
                    total_windows += nwin
                    total_valid_kmers += valid
                    total_skipped_n += skipped
                    total_hits += hits
                    log("Info", f"  {chrom}: windows={nwin}, valid_kmers={valid}, skipped_N/non-ACGT={skipped}, marker_hits={hits}")
            for _chrom, seq in parse_fasta(genome):
                total_bases += len(seq)
    return chrom_count, total_bases, total_windows, total_valid_kmers, total_skipped_n, total_hits

def parse_cigar_blocks(pos0: int, cigar: str) -> Tuple[List[Tuple[int, int]], int, int]:
    """Return query-aligned reference blocks [start,end), reference span start/end."""
    ref = pos0
    span_start = pos0
    blocks: List[Tuple[int, int]] = []
    for length_s, op in CIGAR_RE.findall(cigar):
        length = int(length_s)
        if op in QUERY_ALIGNED:
            if length > 0:
                blocks.append((ref, ref + length))
            ref += length
        elif op in {"D", "N"}:
            ref += length
        elif op in {"I", "S", "H", "P"}:
            continue
        else:
            continue
    return blocks, span_start, ref

def load_chunk_map(path: Optional[str]) -> Dict[str, Tuple[str, int]]:
    chunk_map: Dict[str, Tuple[str, int]] = {}
    if not path:
        return chunk_map
    with open_text(path) as fh:
        next(fh, None)
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chunk, original, start_s, _end_s = parts[:4]
            try:
                chunk_map[chunk] = (original, int(start_s))
            except ValueError:
                continue
    return chunk_map


def add_coverage_event(events: List[Tuple[int, int]], chrom_len: int, start: int, end: int) -> None:
    if end <= start:
        return
    if start < 0:
        start = 0
    if end > chrom_len:
        end = chrom_len
    if end <= start:
        return
    events.append((start, 1))
    events.append((end, -1))


def build_coverage_prefix(events: List[Tuple[int, int]], chrom_len: int) -> Tuple[List[int], List[int], List[int]]:
    if chrom_len <= 0:
        return [0], [0], []
    if not events:
        return [0, chrom_len], [0, 0], [0]

    events.sort()
    points: List[int] = [0]
    prefix: List[int] = [0]
    coverage: List[int] = []
    cov = 0
    prev = 0
    i = 0
    n = len(events)

    while i < n:
        pos = events[i][0]
        if pos < 0:
            pos = 0
        elif pos > chrom_len:
            pos = chrom_len

        if pos > prev:
            coverage.append(cov)
            prefix.append(prefix[-1] + cov * (pos - prev))
            points.append(pos)
            prev = pos

        delta = 0
        while i < n:
            raw_pos, raw_delta = events[i]
            clipped_pos = raw_pos
            if clipped_pos < 0:
                clipped_pos = 0
            elif clipped_pos > chrom_len:
                clipped_pos = chrom_len
            if clipped_pos != pos:
                break
            delta += raw_delta
            i += 1
        cov += delta

        if prev >= chrom_len:
            break

    if chrom_len > prev:
        coverage.append(cov)
        prefix.append(prefix[-1] + cov * (chrom_len - prev))
        points.append(chrom_len)

    if len(points) == 1:
        points.append(chrom_len)
        prefix.append(0)
        coverage.append(0)
    return points, prefix, coverage


def coverage_integral(x: int, points: List[int], prefix: List[int], coverage: List[int]) -> int:
    from bisect import bisect_right
    if x <= points[0]:
        return prefix[0]
    if x >= points[-1]:
        return prefix[-1]
    idx = bisect_right(points, x) - 1
    if idx >= len(coverage):
        return prefix[-1]
    return prefix[idx] + coverage[idx] * (x - points[idx])


def run_sam_to_windows(genome: str, sam_path: str, window: int, step: int, out_path: str, chunk_map_path: Optional[str] = None) -> Tuple[int, int, int, int]:
    order, sizes = load_genome_sizes(genome)
    chunk_map = load_chunk_map(chunk_map_path)
    events_by_chr: Dict[str, List[Tuple[int, int]]] = {chrom: [] for chrom in order}
    total_windows = 0
    for chrom in order:
        length = sizes[chrom]
        total_windows += len(range(0, length, step)) if length > 0 else 0

    total_records = 0
    used_alignments = 0
    skipped_records = 0
    with open_text(sam_path) as fh:
        for raw in fh:
            if not raw or raw.startswith("@"):
                continue
            parts = raw.rstrip("\n").split("\t", 11)
            if len(parts) < 6:
                skipped_records += 1
                continue
            total_records += 1
            try:
                flag = int(parts[1])
                chrom = parts[2]
                pos = int(parts[3])
                cigar = parts[5]
            except Exception:
                skipped_records += 1
                continue
            if (flag & 4) or (flag & 256) or (flag & 2048):
                skipped_records += 1
                continue
            if chrom == "*" or cigar == "*" or pos <= 0:
                skipped_records += 1
                continue
            if chrom in chunk_map:
                original, offset = chunk_map[chrom]
                chrom = original
                pos += offset
            if chrom not in sizes:
                skipped_records += 1
                continue
            pos0 = pos - 1
            blocks, _span_start, _span_end = parse_cigar_blocks(pos0, cigar)
            if not blocks:
                skipped_records += 1
                continue
            chrom_len = sizes[chrom]
            added = False
            events = events_by_chr[chrom]
            for b_start, b_end in blocks:
                before = len(events)
                add_coverage_event(events, chrom_len, b_start, b_end)
                if len(events) > before:
                    added = True
            if not added:
                skipped_records += 1
                continue
            used_alignments += 1

    with open(out_path, "w") as out:
        out.write("chr\tstart\tend\tdepth\n")
        for chrom in order:
            chrom_len = sizes[chrom]
            points, prefix, coverage = build_coverage_prefix(events_by_chr[chrom], chrom_len)
            for start in range(0, chrom_len, step):
                end = start + window
                if end > chrom_len:
                    end = chrom_len
                denom = end - start
                bases = coverage_integral(end, points, prefix, coverage) - coverage_integral(start, points, prefix, coverage)
                depth = (bases / denom) if denom > 0 else 0.0
                out.write(f"{chrom}\t{start}\t{end}\t{depth:.10g}\n")
    return total_records, used_alignments, skipped_records, total_windows

def main() -> int:
    ap = argparse.ArgumentParser(add_help=False)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_kmer = sub.add_parser("kmer", add_help=False)
    p_kmer.add_argument("--genome", required=True)
    p_kmer.add_argument("--markers", required=True)
    p_kmer.add_argument("--kmer-size", type=int, required=True)
    p_kmer.add_argument("--window", type=int, required=True)
    p_kmer.add_argument("--step", type=int, required=True)
    p_kmer.add_argument("--threads", type=int, required=True)
    p_kmer.add_argument("--out", required=True)

    p_sam = sub.add_parser("sam-windows", add_help=False)
    p_sam.add_argument("--genome", required=True)
    p_sam.add_argument("--sam", required=True)
    p_sam.add_argument("--window", type=int, required=True)
    p_sam.add_argument("--step", type=int, required=True)
    p_sam.add_argument("--out", required=True)
    p_sam.add_argument("--chunk-map", default=None)

    args = ap.parse_args()
    try:
        if args.cmd == "kmer":
            if args.window < args.kmer_size:
                raise ValueError(f"Window size ({args.window}) must be >= k-mer size ({args.kmer_size})")
            markers, raw_count, indexed_count = load_markers(args.markers, args.kmer_size)
            log("Info", "2-bit integer marker index built successfully.")
            log("Info", f"  K-mer size      : {args.kmer_size}")
            log("Info", f"  Marker sequences: {raw_count}")
            log("Info", f"  Indexed k-mers   : {indexed_count} (forward + reverse-complement, deduplicated)")
            log("Info", "  Scanner backend : rolling 2-bit integer exact matching")
            chrom_count, total_bases, total_windows, total_valid_kmers, total_skipped_n, total_hits = run_kmer_scan(
                args.genome, markers, args.kmer_size, args.window, args.step, args.threads, args.out
            )
            log("Info", "SEXmer k-mer mapping complete.")
            log("Info", f"  Chromosomes/contigs : {chrom_count}")
            log("Info", f"  Genome bases        : {total_bases}")
            log("Info", f"  Windows written     : {total_windows}")
            log("Info", f"  Valid genome k-mers : {total_valid_kmers}")
            log("Info", f"  Skipped N/non-ACGT  : {total_skipped_n}")
            log("Info", f"  Marker hits         : {total_hits}")
            return 0
        if args.cmd == "sam-windows":
            total, used, skipped, windows = run_sam_to_windows(args.genome, args.sam, args.window, args.step, args.out, args.chunk_map)
            log("Info", "SEXmer read-window table complete.")
            log("Info", f"  SAM records read     : {total}")
            log("Info", f"  Alignments used      : {used}")
            log("Info", f"  Records skipped      : {skipped}")
            log("Info", f"  Windows written      : {windows}")
            return 0
        raise ValueError("Unknown command")
    except Exception as exc:
        log("Error", str(exc))
        return 1

if __name__ == "__main__":
    sys.exit(main())
PY
chmod +x "$SCANNER"

PLOTTER="${MAP_TMPDIR}/sexmer_manhattan_plot.py"
cat > "$PLOTTER" <<'PY'
#!/usr/bin/env python3

from __future__ import annotations

import argparse
import logging
import math
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

DPI = 300
COLORS = ["#0e616d", "#ffc348", "#1497a5", "#f76e68"]
DEFAULT_FORMATS = ("svg",)

class TaggedFormatter(logging.Formatter):
    LEVEL_TAGS = {
        logging.INFO: "[Info]",
        logging.WARNING: "[Warning]",
        logging.ERROR: "[Error]",
    }

    def format(self, record: logging.LogRecord) -> str:
        tag = getattr(record, "tag", None) or self.LEVEL_TAGS.get(record.levelno, "[Info]")
        return f"{tag} {record.getMessage()}"

def setup_logger() -> logging.Logger:
    logger = logging.getLogger("sexmer_manhattan")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(TaggedFormatter())
    logger.addHandler(handler)
    return logger

def log_output(logger: logging.Logger, message: str) -> None:
    logger.info(message, extra={"tag": "[Output]"})

def normalize_columns(columns: Iterable[str]) -> List[str]:
    return [
        c.strip()
        .lower()
        .lstrip("#")
        .replace(" ", "_")
        .replace("-", "_")
        for c in columns
    ]

def detect_input_type(df: pd.DataFrame) -> Tuple[str, str, str]:
    """Return input_type, y_column, y_label."""
    cols = set(df.columns)

    if {"chr", "start", "end", "hits"}.issubset(cols):
        return "kmer", "hits", "Kmer hits"

    if {"chr", "start", "end", "depth"}.issubset(cols):
        return "reads", "depth", "Reads depth"

    if {"chr", "start", "end", "mean_depth"}.issubset(cols):
        return "depth", "mean_depth", "depth"
    if {"chrom", "start", "end", "mean_depth"}.issubset(cols):
        df.rename(columns={"chrom": "chr"}, inplace=True)
        return "depth", "mean_depth", "depth"

    raise ValueError(
        "Input must contain SEXmer map columns: "
        "chr/start/end/hits for k-mer output, or chr/start/end/depth for reads output."
    )

def load_sexmer_table(path: Path) -> Tuple[pd.DataFrame, str, str]:
    df = pd.read_csv(path, sep="\t", comment="#")
    if df.empty:
        raise ValueError(f"Input table is empty: {path}")

    df.columns = normalize_columns(df.columns)
    input_type, y_column, y_label = detect_input_type(df)

    required = ["chr", "start", "end", y_column]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(missing)}")

    for col in ["start", "end", y_column]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["chr", "start", "end", y_column]).copy()
    if df.empty:
        raise ValueError(f"No valid numeric rows found in: {path}")

    df = df[df["end"] > df["start"]].copy()
    if df.empty:
        raise ValueError(f"No valid windows with end > start found in: {path}")

    df["chr"] = df["chr"].astype(str)
    df["pos"] = (df["start"] + df["end"]) / 2.0
    df["plot_value"] = df[y_column].clip(lower=0)
    return df, input_type, y_label

def strip_chr_prefix(chrom: str) -> str:
    chrom_s = str(chrom)
    return chrom_s[3:] if chrom_s.lower().startswith("chr") else chrom_s

def natural_sort_key(chrom: str) -> Tuple[int, float, str]:
    c = strip_chr_prefix(chrom).lower()
    aliases = {"x": 10_000_001, "y": 10_000_002, "m": 10_000_003, "mt": 10_000_003}
    if c.isdigit():
        return (0, int(c), chrom)
    if c in aliases:
        return (1, aliases[c], chrom)
    digits = "".join(ch for ch in c if ch.isdigit())
    if digits:
        return (2, int(digits), chrom)
    return (3, math.inf, chrom)

def prepare_manhattan_data(df: pd.DataFrame) -> Tuple[pd.DataFrame, List[dict]]:
    chroms = sorted(df["chr"].unique(), key=natural_sort_key)
    pieces = []
    chr_info = []
    offset = 0.0

    for idx, chrom in enumerate(chroms):
        chr_df = df[df["chr"] == chrom].copy().sort_values("pos")
        chrom_len = float(chr_df["end"].max())
        start_offset = offset
        end_offset = offset + chrom_len
        center = start_offset + chrom_len / 2.0

        chr_df["cum_pos"] = chr_df["pos"] + offset
        chr_df["chr_idx"] = idx
        pieces.append(chr_df)
        chr_info.append(
            {
                "chr": chrom,
                "start": start_offset,
                "end": end_offset,
                "center": center,
                "color": COLORS[idx % len(COLORS)],
            }
        )
        offset = end_offset

    if not pieces:
        raise ValueError("No chromosome/window data available for plotting")
    return pd.concat(pieces, ignore_index=True), chr_info

def nice_y_upper(values: pd.Series) -> float:
    """Return a rounded y-axis upper bound with a labeled top tick.
    """
    max_val = float(values.max()) if len(values) else 0.0
    if not math.isfinite(max_val) or max_val <= 0:
        return 1.0

    exponent = math.floor(math.log10(max_val))
    fraction = max_val / (10 ** exponent)

    for nice_fraction in (1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0):
        if fraction <= nice_fraction:
            return nice_fraction * (10 ** exponent)

    return 10.0 * (10 ** exponent)

def choose_figsize(chr_count: int, width: Optional[float], height: Optional[float]) -> Tuple[float, float]:
    if width is not None or height is not None:
        return (width if width is not None else 16.0, height if height is not None else 5.3)
    auto_width = min(max(14.0, chr_count * 0.55), 28.0)
    return (auto_width, 5.3)

def create_plot(
    df: pd.DataFrame,
    chr_info: Sequence[dict],
    output_file: Path,
    y_label: str,
    plot_type: str,
    width: Optional[float],
    height: Optional[float],
) -> None:
    figsize = choose_figsize(len(chr_info), width, height)
    fig, ax = plt.subplots(figsize=figsize, dpi=DPI)

    for item in chr_info:
        chrom = item["chr"]
        chr_df = df[df["chr"] == chrom].sort_values("cum_pos")
        if chr_df.empty:
            continue
        x = chr_df["cum_pos"] / 1e6
        y = chr_df["plot_value"]
        if plot_type == "line":
            ax.plot(x, y, color=item["color"], linewidth=0.9, alpha=0.9, rasterized=True)
        elif plot_type == "bar":
            step = chr_df["cum_pos"].diff().median()
            width_mb = (step if pd.notna(step) and step > 0 else 100_000) / 1e6
            ax.bar(x, y, color=item["color"], width=width_mb, alpha=0.9, edgecolor="none", rasterized=True)
        else:
            ax.scatter(x, y, c=item["color"], s=10, alpha=0.75, edgecolors="none", rasterized=True)

    tick_positions = [0.0] + [item["end"] / 1e6 for item in chr_info]
    tick_labels = [""] + [strip_chr_prefix(item["chr"]) for item in chr_info]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontweight="bold", rotation=0)
    ax.set_xlim(left=0, right=max(tick_positions) if tick_positions else None)

    ax.set_xlabel("Chromosome", fontsize=18, labelpad=10)
    ax.set_ylabel(y_label, fontsize=18, labelpad=10)

    y_upper = nice_y_upper(df["plot_value"])

    y_lower = -0.02 * y_upper
    ax.set_ylim(y_lower, y_upper)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, prune=None, steps=[1, 2, 2.5, 5, 10]))
    ticks = [tick for tick in ax.get_yticks() if 0 <= tick < y_upper]
    if 0.0 not in ticks:
        ticks.insert(0, 0.0)
    if not ticks or ticks[-1] != y_upper:
        ticks.append(y_upper)
    ax.set_yticks(ticks)

    ax.tick_params(axis="x", which="major", labelsize=12, length=6, width=1, direction="out", pad=8)
    ax.tick_params(axis="y", which="major", labelsize=13, length=6, width=1, direction="out", pad=8)

    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_position(("outward", 8))
    ax.spines["bottom"].set_position(("outward", 5))
    ax.spines["left"].set_bounds(0, y_upper)

    fig.tight_layout()
    fig.savefig(output_file, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def default_output_stem(input_path: Path, input_type: str) -> str:
    stem = input_path.name
    for suffix in (".kmer.windows.tsv", ".reads.windows.tsv", ".windows.tsv", ".tsv"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return f"{stem}.{input_type}.manhattan"

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create standard Manhattan plots for SEXmer map window TSV outputs."
    )
    parser.add_argument("inputs", nargs="+", help="SEXmer map TSV file(s)")
    return parser.parse_args()

def main() -> int:
    args = parse_args()
    logger = setup_logger()

    try:
        for input_name in args.inputs:
            input_path = Path(input_name)
            if not input_path.exists():
                raise FileNotFoundError(f"Input file not found: {input_path}")

            logger.info(f"Reading input data: {input_path}")
            raw_df, input_type, y_label = load_sexmer_table(input_path)
            logger.info(f"Detected input type: {input_type}; Y-axis: {y_label}")

            logger.info("Preparing Manhattan plot data")
            plot_df, chr_info = prepare_manhattan_data(raw_df)

            output_file = input_path.parent / f"{default_output_stem(input_path, input_type)}.svg"
            logger.info("Generating SVG Manhattan plot")
            create_plot(
                df=plot_df,
                chr_info=chr_info,
                output_file=output_file,
                y_label=y_label,
                plot_type="scatter",
                width=None,
                height=None,
            )
            log_output(logger, f"Plot file: {output_file}")

        logger.info("Done")
        return 0
    except Exception as exc:
        logger.error(str(exc))
        return 1

if __name__ == "__main__":
    sys.exit(main())

PY
chmod +x "$PLOTTER"

KMER_OUT="${PREFIX}.kmer.windows.tsv"
READS_OUT="${PREFIX}.reads.windows.tsv"

info "SEXmer map starting"
info "Parameters : kmer-size=${KMER_SIZE}, window=${WINDOW}, step=${STEP}"
info "Settings   : threads=${THREADS}"
info "Genome     : ${GENOME}"
info "Markers    : ${MARKERS}"
info "Kmer output: ${KMER_OUT}"
info "Temp dir   : ${MAP_TMPDIR}"
info "Coordinate : 0-based, half-open windows"
info "Strand mode: forward + reverse-complement marker lookup"
info "Normalize  : hits_per_10kb = hits / valid_ACGT_kmer_sites * 10000"

info "Running coordinate-aware k-mer window scanner..."
python3 "$SCANNER" kmer \
    --genome "$GENOME" \
    --markers "$MARKERS" \
    --kmer-size "$KMER_SIZE" \
    --window "$WINDOW" \
    --step "$STEP" \
    --threads "$THREADS" \
    --out "$KMER_OUT"

[[ -f "$KMER_OUT" ]] || { error "K-mer scanner did not produce output file: $KMER_OUT"; exit 1; }
output "K-mer window table written to: ${KMER_OUT}"

info "Generating Manhattan plot for k-mer window table..."
python3 "$PLOTTER" "$KMER_OUT"
KMER_PLOT="$(dirname "$KMER_OUT")/$(basename "$KMER_OUT" .kmer.windows.tsv).kmer.manhattan.svg"
[[ -f "$KMER_PLOT" ]] || { error "Manhattan plotter did not produce output file: $KMER_PLOT"; exit 1; }
output "K-mer Manhattan plot written to: ${KMER_PLOT}"

filename_pair_score() {
    local path1="$1" path2="$2"
    local b1 b2 stem1 stem2
    b1="$(basename "$path1")"
    b2="$(basename "$path2")"
    stem1="$b1"; stem2="$b2"
    local suffix
    for suffix in ".fastq.gz" ".fq.gz" ".fastq" ".fq" ".gz"; do
        [[ "$stem1" == *"$suffix" ]] && stem1="${stem1%$suffix}"
        [[ "$stem2" == *"$suffix" ]] && stem2="${stem2%$suffix}"
    done
    if [[ "$stem1" =~ ^(.+)_1$ && "${BASH_REMATCH[1]}_2" == "$stem2" ]]; then return 0; fi
    if [[ "$stem1" =~ ^(.+)_R1$ && "${BASH_REMATCH[1]}_R2" == "$stem2" ]]; then return 0; fi
    if [[ "$stem1" =~ ^(.+)_R1_001$ && "${BASH_REMATCH[1]}_R2_001" == "$stem2" ]]; then return 0; fi
    if [[ "$stem1" =~ ^(.+)\.1$ && "${BASH_REMATCH[1]}.2" == "$stem2" ]]; then return 0; fi
    if [[ "$stem1" =~ ^(.+)-1$ && "${BASH_REMATCH[1]}-2" == "$stem2" ]]; then return 0; fi
    return 1
}

if [[ ${#READ_FILES[@]} -gt 0 ]]; then
    info "Read validation requested."
    info "Read files : ${READ_FILES[*]}"
    info "Seq type   : ${SEQ_TYPE}"
    info "Read output: ${READS_OUT}"

    SAM_TMP="${MAP_TMPDIR}/bbmap_stream.sam"
    BBMAP_INDEX_DIR="${MAP_TMPDIR}/bbmap_ref"
    mkdir -p "$BBMAP_INDEX_DIR"
    : > "$SAM_TMP"

    BBMAP_REF_GENOME="${MAP_TMPDIR}/bbmap_reference.fa"
    BBMAP_REF_CHUNK_MAP="${MAP_TMPDIR}/bbmap_reference.chunk_map.tsv"
    BBMAP_CHUNK_SIZE=500000000
    BBMAP_CHUNK_OVERLAP=10000

    info "Preparing BBMap-compatible reference for read validation..."
    python3 - "$GENOME" "$BBMAP_REF_GENOME" "$BBMAP_REF_CHUNK_MAP" "$BBMAP_CHUNK_SIZE" "$BBMAP_CHUNK_OVERLAP" <<'PY'
import gzip
import sys

genome, out_fa, map_tsv, chunk_size_s, overlap_s = sys.argv[1:6]
chunk_size = int(chunk_size_s)
overlap = int(overlap_s)
if overlap < 0 or overlap >= chunk_size:
    raise SystemExit("Invalid BBMap chunk overlap")

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def wrap_write(out, seq, width=80):
    for i in range(0, len(seq), width):
        out.write(seq[i:i+width] + "\n")

def emit_record(out, mp, name, seq):
    length = len(seq)
    if length == 0:
        return 0, 0
    if length <= chunk_size:
        out.write(f">{name}\n")
        wrap_write(out, seq)
        return 1, 0
    chunks = 0
    step = chunk_size - overlap
    start = 0
    while start < length:
        end = min(start + chunk_size, length)
        chunk_name = f"{name}__SEXMERCHUNK__{start}__{end}"
        out.write(f">{chunk_name}\n")
        wrap_write(out, seq[start:end])
        mp.write(f"{chunk_name}\t{name}\t{start}\t{end}\n")
        chunks += 1
        if end == length:
            break
        start += step
    return chunks, 1

records = 0
chunks = 0
split_records = 0
name = None
parts = []
with open_text(genome) as fh, open(out_fa, "w") as out, open(map_tsv, "w") as mp:
    mp.write("chunk\toriginal\tstart\tend\n")
    for raw in fh:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                c, sp = emit_record(out, mp, name, "".join(parts))
                records += 1
                chunks += c
                split_records += sp
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)
    if name is not None:
        c, sp = emit_record(out, mp, name, "".join(parts))
        records += 1
        chunks += c
        split_records += sp

if records == 0 or chunks == 0:
    raise SystemExit(f"No sequence found while preparing BBMap reference: {genome}")
print(f"[Info] BBMap reference records: original={records}, bbmap_records={chunks}, split_large_records={split_records}", file=sys.stderr)
PY

    LONG_READ_CHUNK=6000
    MAPPER="bbmap.sh"
    if [[ "$SEQ_TYPE" == "long" ]]; then
        if command -v mapPacBio.sh &>/dev/null; then
            MAPPER="mapPacBio.sh"
            info "Long-read mode selected; using mapPacBio.sh after read splitting."
        else
            warn "mapPacBio.sh not found; falling back to bbmap.sh after read splitting."
        fi
    fi

    split_long_reads() {
        local input="$1"
        local output="$2"
        info "Splitting long reads into <=${LONG_READ_CHUNK} bp pieces: $(basename "$input")"
        reformat.sh in="$input" out="$output" breaklength="$LONG_READ_CHUNK" overwrite=t ziplevel=2 2>&1             | while IFS= read -r line; do info "  [reformat] $line"; done
        [[ -s "$output" ]] || { error "Long-read splitting did not produce output file: $output"; exit 1; }
    }

    run_bbmap_to_sam() {
        local r1="$1"
        local r2="${2:-}"
        local out_sam="$3"
        if [[ -n "$r2" ]]; then
            info "Running BBMap paired-end validation: $(basename "$r1") , $(basename "$r2")"
            "$MAPPER" ref="$BBMAP_REF_GENOME" path="$BBMAP_INDEX_DIR" in="$r1" in2="$r2" out="$out_sam" overwrite=t threads="$THREADS" 2>&1                 | while IFS= read -r line; do info "  [bbmap] $line"; done
        else
            info "Running BBMap single-read validation: $(basename "$r1")"
            "$MAPPER" ref="$BBMAP_REF_GENOME" path="$BBMAP_INDEX_DIR" in="$r1" out="$out_sam" overwrite=t threads="$THREADS" 2>&1                 | while IFS= read -r line; do info "  [bbmap] $line"; done
        fi
    }

    if [[ "$SEQ_TYPE" == "short" && ${#READ_FILES[@]} -eq 2 ]] && filename_pair_score "${READ_FILES[0]}" "${READ_FILES[1]}"; then
        info "Read mode detected: paired-end based on filename pattern."
        run_bbmap_to_sam "${READ_FILES[0]}" "${READ_FILES[1]}" "$SAM_TMP"
    else
        if [[ "$SEQ_TYPE" == "long" ]]; then
            info "Read mode detected: long-read file set (${#READ_FILES[@]} file(s)); mapping each split file independently."
        elif [[ ${#READ_FILES[@]} -eq 2 ]]; then
            info "Read mode detected: two single-read files; filename pattern is not paired-end."
        else
            info "Read mode detected: single-read file set (${#READ_FILES[@]} file(s))."
        fi
        FIRST=1
        for rf in "${READ_FILES[@]}"; do
            MAP_INPUT="$rf"
            if [[ "$SEQ_TYPE" == "long" ]]; then
                SPLIT_FASTQ="${MAP_TMPDIR}/long_reads_${FIRST}.split.fq.gz"
                split_long_reads "$rf" "$SPLIT_FASTQ"
                MAP_INPUT="$SPLIT_FASTQ"
            fi
            TMP_ONE="${MAP_TMPDIR}/bbmap_one_${FIRST}.sam"
            info "Running BBMap for read file ${FIRST}/${#READ_FILES[@]}: $(basename "$MAP_INPUT")"
            run_bbmap_to_sam "$MAP_INPUT" "" "$TMP_ONE"
            if [[ "$FIRST" -eq 1 ]]; then
                cat "$TMP_ONE" > "$SAM_TMP"
            else
                awk 'BEGIN{OFS="\t"} /^@/ {next} {print}' "$TMP_ONE" >> "$SAM_TMP"
            fi
            FIRST=$(( FIRST + 1 ))
        done
    fi

    [[ -s "$SAM_TMP" ]] || { error "BBMap did not produce SAM output for read validation."; exit 1; }

    SAM_FOR_WINDOWS="$SAM_TMP"

    info "Converting BBMap SAM output to sliding-window read depth table..."
    python3 "$SCANNER" sam-windows \
        --genome "$GENOME" \
        --sam "$SAM_FOR_WINDOWS" \
        --chunk-map "$BBMAP_REF_CHUNK_MAP" \
        --window "$WINDOW" \
        --step "$STEP" \
        --out "$READS_OUT"

    [[ -f "$READS_OUT" ]] || { error "Read-window parser did not produce output file: $READS_OUT"; exit 1; }
    output "Read validation window table written to: ${READS_OUT}"

    info "Generating Manhattan plot for read-window table..."
    python3 "$PLOTTER" "$READS_OUT"
    READS_PLOT="$(dirname "$READS_OUT")/$(basename "$READS_OUT" .reads.windows.tsv).reads.manhattan.svg"
    [[ -f "$READS_PLOT" ]] || { error "Manhattan plotter did not produce output file: $READS_PLOT"; exit 1; }
    output "Read Manhattan plot written to: ${READS_PLOT}"
fi

info "SEXmer map complete."