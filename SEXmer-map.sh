#!/usr/bin/env bash
# SEXmer-map.v2.sh - Map SEXmer marker k-mers and optional extracted reads to reference genome windows.
# Author: Dede Kurniawan

set -euo pipefail
export LC_ALL=C

# defaults
KMER_SIZE=21
WINDOW=10000
STEP=2500
THREADS=8
PREFIX=""
TMPDIR_BASE="$(pwd)"
GENOME=""
MARKERS=""
READS_INPUT=""
SEQ_TYPE="auto"

# log helpers
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*" >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"  >&2; }

usage() {
    cat >&2 <<EOF

SEXmer-map.sh - Map SEXmer marker k-mers and optional extracted reads to reference genome windows.

Usage: SEXmer-map.sh <genome.fa> <markers.fa> --prefix <prefix> [OPTIONS]

Mandatory:
  <genome.fa>          Reference genome FASTA (.fa, .fasta, optionally .gz)
  <markers.fa>         Marker k-mer FASTA, e.g. MSK.fa or FSK.fa (.gz accepted)
  --prefix             Output filename prefix

Optionals:
  -k, --kmer-size      K-mer size (1-63)                               [default: ${KMER_SIZE}]
  -w, --window         Window size in bp                               [default: ${WINDOW}]
  -s, --step           Sliding step size in bp                         [default: ${STEP}]
  -r, --reads          Comma-separated FASTQ file(s) for BBMap validation
                       Examples: -r reads.fq.gz OR -r reads_1.fq.gz,reads_2.fq.gz
  --seq-type           Read technology: auto, short, ONT, PacBio       [default: ${SEQ_TYPE}]
  -t, --threads        CPU threads                                     [default: ${THREADS}]
  --tmpdir             Parent directory for temporary work folder      [default: current dir]
  -h, --help           Show this help and exit

EOF
    exit 1
}

[[ $# -eq 0 ]] && usage

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
        --tmpdir)          TMPDIR_BASE="$2"; shift 2 ;;
        -h|--help)         usage ;;
        -*) error "Unknown option '$1'"; usage ;;
        *)  POSITIONAL+=("$1"); shift ;;
    esac
done

[[ ${#POSITIONAL[@]} -eq 2 ]] || {
    error "Exactly two positional arguments are required: <genome.fa> <markers.fa>."; usage; }

GENOME="${POSITIONAL[0]}"
MARKERS="${POSITIONAL[1]}"

[[ -z "$PREFIX" ]] && { error "--prefix is required."; usage; }

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
    auto|short|ONT|ont|PacBio|pacbio|PACBIO) ;;
    *) error "--seq-type must be one of: auto, short, ONT, PacBio."; exit 1 ;;
esac
case "$SEQ_TYPE" in
    ont) SEQ_TYPE="ONT" ;;
    pacbio|PACBIO) SEQ_TYPE="PacBio" ;;
esac

[[ -d "$TMPDIR_BASE" ]] || { error "Temporary parent directory does not exist: ${TMPDIR_BASE}"; exit 1; }
[[ -w "$TMPDIR_BASE" ]] || { error "Temporary parent directory is not writable: ${TMPDIR_BASE}"; exit 1; }

OUT_DIR="$(dirname "$PREFIX")"
[[ "$OUT_DIR" == "." ]] && OUT_DIR="$(pwd)"
[[ -d "$OUT_DIR" ]] || { error "Output directory does not exist: ${OUT_DIR}"; exit 1; }
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

def add_interval_bases(arr: List[int], chrom_len: int, start: int, end: int, window: int, step: int) -> None:
    if end <= start or not arr:
        return
    if start < 0:
        start = 0
    if end > chrom_len:
        end = chrom_len
    if end <= start:
        return
    # Windows overlap interval if win_start < end and win_end > start.
    i_min = 0 if start < window else (start - window) // step + 1
    i_max = (end - 1) // step
    if i_max >= len(arr):
        i_max = len(arr) - 1
    for idx in range(i_min, i_max + 1):
        w_start = idx * step
        w_end = w_start + window
        if w_end > chrom_len:
            w_end = chrom_len
        ov = min(end, w_end) - max(start, w_start)
        if ov > 0:
            arr[idx] += ov

def add_interval_hit(arr: List[int], chrom_len: int, start: int, end: int, window: int, step: int) -> None:
    if end <= start or not arr:
        return
    if start < 0:
        start = 0
    if end > chrom_len:
        end = chrom_len
    if end <= start:
        return
    i_min = 0 if start < window else (start - window) // step + 1
    i_max = (end - 1) // step
    if i_max >= len(arr):
        i_max = len(arr) - 1
    for idx in range(i_min, i_max + 1):
        w_start = idx * step
        w_end = w_start + window
        if w_end > chrom_len:
            w_end = chrom_len
        if min(end, w_end) > max(start, w_start):
            arr[idx] += 1

def run_sam_to_windows(genome: str, sam_path: str, window: int, step: int, out_path: str) -> Tuple[int, int, int, int]:
    order, sizes = load_genome_sizes(genome)
    nwin_by_chr: Dict[str, int] = {}
    bases_by_chr: Dict[str, List[int]] = {}
    hits_by_chr: Dict[str, List[int]] = {}
    total_windows = 0
    for chrom in order:
        length = sizes[chrom]
        if length <= 0:
            nwin = 0
        else:
            nwin = len(range(0, length, step))
        nwin_by_chr[chrom] = nwin
        bases_by_chr[chrom] = [0] * nwin
        hits_by_chr[chrom] = [0] * nwin
        total_windows += nwin

    total_records = 0
    used_alignments = 0
    skipped_records = 0
    with open_text(sam_path) as fh:
        for raw in fh:
            if not raw or raw.startswith("@"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 11:
                skipped_records += 1
                continue
            total_records += 1
            try:
                flag = int(parts[1])
                chrom = parts[2]
                pos = int(parts[3])
                mapq = int(parts[4])
                cigar = parts[5]
            except Exception:
                skipped_records += 1
                continue
            # Skip unmapped, secondary, and supplementary alignments to avoid double-counting.
            if (flag & 4) or (flag & 256) or (flag & 2048):
                skipped_records += 1
                continue
            if chrom == "*" or cigar == "*" or chrom not in sizes or pos <= 0:
                skipped_records += 1
                continue
            pos0 = pos - 1
            blocks, span_start, span_end = parse_cigar_blocks(pos0, cigar)
            if span_end <= span_start or not blocks:
                skipped_records += 1
                continue
            chrom_len = sizes[chrom]
            for b_start, b_end in blocks:
                add_interval_bases(bases_by_chr[chrom], chrom_len, b_start, b_end, window, step)
            add_interval_hit(hits_by_chr[chrom], chrom_len, span_start, span_end, window, step)
            used_alignments += 1

    with open(out_path, "w") as out:
        out.write("chr\tstart\tend\tdepth\tread_hits\n")
        for chrom in order:
            chrom_len = sizes[chrom]
            nwin = nwin_by_chr[chrom]
            bases = bases_by_chr[chrom]
            hits = hits_by_chr[chrom]
            for idx in range(nwin):
                start = idx * step
                end = start + window
                if end > chrom_len:
                    end = chrom_len
                denom = end - start
                depth = (bases[idx] / denom) if denom > 0 else 0.0
                out.write(f"{chrom}\t{start}\t{end}\t{depth:.10g}\t{hits[idx]}\n")
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
            total, used, skipped, windows = run_sam_to_windows(args.genome, args.sam, args.window, args.step, args.out)
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

KMER_OUT="${PREFIX}.kmer.windows.tsv"
READS_OUT="${PREFIX}.reads.windows.tsv"

info "SEXmer-map starting"
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

    MAPPER="bbmap.sh"
    if [[ "$SEQ_TYPE" == "ONT" || "$SEQ_TYPE" == "PacBio" ]]; then
        if command -v mapPacBio.sh &>/dev/null; then
            MAPPER="mapPacBio.sh"
            info "Long-read seq-type selected; using mapPacBio.sh for BBTools long-read mapping."
        else
            warn "mapPacBio.sh not found; falling back to bbmap.sh for long-read validation."
        fi
    fi

    run_bbmap_to_sam() {
        local r1="$1"
        local r2="${2:-}"
        if [[ -n "$r2" ]]; then
            info "Running BBMap paired-end validation: $(basename "$r1") , $(basename "$r2")"
            "$MAPPER" ref="$GENOME" path="$BBMAP_INDEX_DIR" in="$r1" in2="$r2" out="$SAM_TMP" overwrite=t threads="$THREADS" 2>&1 \
                | while IFS= read -r line; do info "  [bbmap] $line"; done
        else
            info "Running BBMap single/long-read validation: $(basename "$r1")"
            "$MAPPER" ref="$GENOME" path="$BBMAP_INDEX_DIR" in="$r1" out="$SAM_TMP" overwrite=t threads="$THREADS" 2>&1 \
                | while IFS= read -r line; do info "  [bbmap] $line"; done
        fi
    }

    if [[ ${#READ_FILES[@]} -eq 2 ]] && filename_pair_score "${READ_FILES[0]}" "${READ_FILES[1]}"; then
        info "Read mode detected: paired-end based on filename pattern."
        run_bbmap_to_sam "${READ_FILES[0]}" "${READ_FILES[1]}"
    else
        if [[ ${#READ_FILES[@]} -eq 2 ]]; then
            info "Read mode detected: two single/long-read files; filename pattern is not paired-end."
        else
            info "Read mode detected: single/long-read file set (${#READ_FILES[@]} file(s))."
        fi
        FIRST=1
        for rf in "${READ_FILES[@]}"; do
            TMP_ONE="${MAP_TMPDIR}/bbmap_one_${FIRST}.sam"
            info "Running BBMap for read file ${FIRST}/${#READ_FILES[@]}: $(basename "$rf")"
            "$MAPPER" ref="$GENOME" path="$BBMAP_INDEX_DIR" in="$rf" out="$TMP_ONE" overwrite=t threads="$THREADS" 2>&1 \
                | while IFS= read -r line; do info "  [bbmap] $line"; done
            if [[ "$FIRST" -eq 1 ]]; then
                cat "$TMP_ONE" > "$SAM_TMP"
            else
                awk 'BEGIN{OFS="\t"} /^@/ {next} {print}' "$TMP_ONE" >> "$SAM_TMP"
            fi
            FIRST=$(( FIRST + 1 ))
        done
    fi

    [[ -s "$SAM_TMP" ]] || { error "BBMap did not produce SAM output for read validation."; exit 1; }

    info "Converting BBMap SAM output to sliding-window read depth table..."
    python3 "$SCANNER" sam-windows \
        --genome "$GENOME" \
        --sam "$SAM_TMP" \
        --window "$WINDOW" \
        --step "$STEP" \
        --out "$READS_OUT"

    [[ -f "$READS_OUT" ]] || { error "Read-window parser did not produce output file: $READS_OUT"; exit 1; }
    output "Read validation window table written to: ${READS_OUT}"
fi

info "SEXmer-map complete."
