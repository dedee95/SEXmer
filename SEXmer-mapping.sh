#!/usr/bin/env bash
# SEXmer-mapping.sh - Map SEXmer marker k-mers to reference genome windows.
# Author: Dede Kurniawan

set -euo pipefail
export LC_ALL=C

# defaults
KMER_SIZE="auto"
WINDOW=10000
STEP=2500
COUNT_MODE="hit"
MEM="16G"
THREADS=8
PREFIX=""
TMPDIR_BASE="$(pwd)"
GENOME=""
MARKERS=""

# log helpers
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*" >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"  >&2; }

usage() {
    cat >&2 <<EOF

SEXmer-mapping.sh - Map SEXmer marker k-mers to reference genome windows.

Usage: SEXmer-mapping.sh <genome.fa> <markers.fa> --prefix <prefix> [OPTIONS]

Mandatory:
  <genome.fa>          Reference genome FASTA (.fa, .fasta, optionally .gz)
  <markers.fa>         Marker k-mer FASTA, e.g. MSK.fa or FSK.fa
  --prefix             Output filename prefix (output: <prefix>.windows.tsv)

Optional:
  -k, --kmer-size      K-mer size; auto-detected from markers if unset [default: ${KMER_SIZE}]
  -w, --window         Window size in bp                               [default: ${WINDOW}]
  -s, --step           Sliding step size in bp                         [default: ${STEP}]
  --count-mode         Count mode; currently only 'hit'                [default: ${COUNT_MODE}]
  --mem                RAM budget for logging/compatibility            [default: ${MEM}]
  -t, --threads        CPU threads for future/compatibility use        [default: ${THREADS}]
  --tmpdir             Parent directory for temporary work folder      [default: current dir]
  -h, --help           Show this help and exit

EOF
    exit 1
}

[[ $# -eq 0 ]] && usage

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)          PREFIX="$2";      shift 2 ;;
        -k|--kmer-size)    KMER_SIZE="$2";  shift 2 ;;
        -w|--window)       WINDOW="$2";     shift 2 ;;
        -s|--step)         STEP="$2";       shift 2 ;;
        --count-mode)      COUNT_MODE="$2"; shift 2 ;;
        --mem)             MEM="$2";        shift 2 ;;
        -t|--threads)      THREADS="$2";    shift 2 ;;
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

WINDOW="$(normalize_int "$WINDOW")"
STEP="$(normalize_int "$STEP")"
THREADS="$(normalize_int "$THREADS")"

if [[ "$KMER_SIZE" != "auto" ]]; then
    KMER_SIZE="$(normalize_int "$KMER_SIZE")"
    [[ "$KMER_SIZE" =~ ^[1-9][0-9]*$ ]] && [[ "$KMER_SIZE" -le 63 ]] || {
        error "--kmer-size must be 'auto' or an integer between 1 and 63."; exit 1; }
fi

[[ "$WINDOW" =~ ^[1-9][0-9]*$ ]] || {
    error "--window must be a positive integer."; exit 1; }

[[ "$STEP" =~ ^[1-9][0-9]*$ ]] || {
    error "--step must be a positive integer."; exit 1; }

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    error "--threads must be a positive integer."; exit 1; }

[[ "$MEM" =~ ^[0-9_]+[GgMm]$ ]] || {
    error "--mem must be a number followed by G or M (e.g. 16G, 512M)."; exit 1; }
MEM="${MEM//_/}"

[[ "$COUNT_MODE" == "hit" ]] || {
    error "--count-mode currently supports only 'hit'."; exit 1; }

[[ -d "$TMPDIR_BASE" ]] || { error "Temporary parent directory does not exist: ${TMPDIR_BASE}"; exit 1; }
[[ -w "$TMPDIR_BASE" ]] || { error "Temporary parent directory is not writable: ${TMPDIR_BASE}"; exit 1; }

OUT_DIR="$(dirname "$PREFIX")"
[[ "$OUT_DIR" == "." ]] && OUT_DIR="$(pwd)"
[[ -d "$OUT_DIR" ]] || { error "Output directory does not exist: ${OUT_DIR}"; exit 1; }
[[ -w "$OUT_DIR" ]] || { error "Output directory is not writable: ${OUT_DIR}"; exit 1; }

command -v python3 &>/dev/null || { error "python3 not found on PATH."; exit 1; }

MAPPING_TMPDIR="${TMPDIR_BASE}/sexmer_mapping_tmp_$$"
mkdir -p "$MAPPING_TMPDIR"
cleanup() { rm -rf "$MAPPING_TMPDIR"; }
trap cleanup EXIT

SCANNER="${MAPPING_TMPDIR}/sexmer_mapping_scanner.py"
cat > "$SCANNER" <<'PY'
#!/usr/bin/env python3
import argparse
import gzip
import multiprocessing as mp
import os
import sys
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

# 2-bit encoding: A=0, C=1, G=2, T=3.
# Non-ACGT bases reset the rolling k-mer state, so k-mers crossing N/ambiguous
# bases are never counted.
BASE_CODE: Dict[str, int] = {
    "A": 0, "C": 1, "G": 2, "T": 3,
    "a": 0, "c": 1, "g": 2, "t": 3,
}
RC_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")

# Globals used by multiprocessing workers. They are initialized once per worker.
_G_MARKERS: Optional[set] = None
_G_K: Optional[int] = None
_G_WINDOW: Optional[int] = None
_G_STEP: Optional[int] = None


def log(kind: str, msg: str) -> None:
    print(f"[{kind}] {msg}", file=sys.stderr)


def open_text(path: str):
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


def load_markers(path: str, k_arg: str) -> Tuple[set, int, int, int]:
    raw_count = 0
    indexed = set()
    detected_k = None

    for name, seq in parse_fasta(path):
        if not seq:
            raise ValueError(f"Empty marker sequence in FASTA record: {name}")
        seq = seq.strip()
        if detected_k is None:
            detected_k = len(seq) if k_arg == "auto" else int(k_arg)
            if detected_k < 1 or detected_k > 63:
                raise ValueError("k-mer size must be between 1 and 63")
        if len(seq) != detected_k:
            raise ValueError(
                f"Marker length mismatch: expected {detected_k}, found {len(seq)} in record '{name}'"
            )
        code = encode_kmer(seq)
        if code is None:
            raise ValueError(f"Marker contains non-ACGT base in record '{name}'")
        rc_code = encode_kmer(revcomp(seq))
        if rc_code is None:
            raise ValueError(f"Marker contains non-ACGT base in record '{name}'")
        raw_count += 1
        indexed.add(code)
        indexed.add(rc_code)

    if raw_count == 0 or detected_k is None:
        raise ValueError(f"No marker sequence found in FASTA: {path}")

    return indexed, detected_k, raw_count, len(indexed)


def window_starts(seq_len: int, k: int, step: int) -> List[int]:
    if seq_len < k:
        return []
    # Report windows that can contain at least one k-mer start position.
    return list(range(0, seq_len - k + 1, step))


def count_valid_sites_for_window(
    q_start: int,
    q_end: int,
    intervals: Sequence[Tuple[int, int]],
    start_hint: int,
) -> Tuple[int, int]:
    """Count valid k-mer start sites in [q_start, q_end] inclusive.

    intervals are sorted inclusive [start, end] valid k-mer-start intervals
    derived from contiguous ACGT runs. start_hint lets consecutive windows avoid
    rescanning old intervals.
    """
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


def scan_chrom_int2bit(
    chrom: str,
    seq: str,
    markers: set,
    k: int,
    window: int,
    step: int,
) -> Tuple[str, int, int, int, int, str]:
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

        # Local variable binding reduces Python lookup overhead in the inner loop.
        base_code_get = BASE_CODE.get
        marker_set = markers
        counts_local = counts
        step_local = step
        nwin_last = nwin - 1
        k_local = k
        mask_local = mask
        span_local = span
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

            # Count this k-mer start position in every overlapping window where:
            # window_start <= kmer_start <= window_end - k.
            left = kmer_start - span_local
            if left <= 0:
                i_min = 0
            else:
                i_min = (left + step_local - 1) // step_local
            i_max = kmer_start // step_local
            if i_max > nwin_last:
                i_max = nwin_last

            for idx in range(i_min, i_max + 1):
                # For the last partial window, confirm the k-mer start is inside
                # its effective k-mer range.
                start = starts[idx]
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
        valid_sites, interval_hint = count_valid_sites_for_window(
            start, end - k, valid_intervals, interval_hint
        )
        hits = counts[idx]
        density = (hits / valid_sites) if valid_sites > 0 else 0.0
        hits_per_10kb = density * 10000.0
        rows.append(
            f"{chrom}\t{start}\t{end}\t{hits}\t{valid_sites}\t{density:.10g}\t{hits_per_10kb:.10g}\n"
        )

    return chrom, nwin, valid_kmers, skipped_non_acgt, total_hits, "".join(rows)

def init_worker(markers: set, k: int, window: int, step: int) -> None:
    global _G_MARKERS, _G_K, _G_WINDOW, _G_STEP
    _G_MARKERS = markers
    _G_K = k
    _G_WINDOW = window
    _G_STEP = step


def scan_chrom_worker(item: Tuple[str, str]) -> Tuple[str, int, int, int, int, str]:
    chrom, seq = item
    assert _G_MARKERS is not None
    assert _G_K is not None
    assert _G_WINDOW is not None
    assert _G_STEP is not None
    return scan_chrom_int2bit(chrom, seq, _G_MARKERS, _G_K, _G_WINDOW, _G_STEP)


def run_scan(
    genome: str,
    markers: set,
    k: int,
    window: int,
    step: int,
    threads: int,
    out_path: str,
) -> Tuple[int, int, int, int, int, int]:
    chrom_count = 0
    total_bases = 0
    total_windows = 0
    total_valid_kmers = 0
    total_skipped_n = 0
    total_hits = 0

    with open(out_path, "w") as out:
        out.write("chrom\tstart\tend\thits\tvalid_kmers\tdensity\thits_per_10kb\n")

        if threads <= 1:
            for chrom, seq in parse_fasta(genome):
                chrom_count += 1
                total_bases += len(seq)
                log("Info", f"Scanning {chrom} ({len(seq)} bp)...")
                _chrom, nwin, valid, skipped, hits, rows = scan_chrom_int2bit(
                    chrom, seq, markers, k, window, step
                )
                out.write(rows)
                total_windows += nwin
                total_valid_kmers += valid
                total_skipped_n += skipped
                total_hits += hits
                log("Info", f"  {chrom}: windows={nwin}, valid_kmers={valid}, skipped_N/non-ACGT={skipped}, marker_hits={hits}")
        else:
            log("Info", f"Using chromosome/contig-level multiprocessing with {threads} worker(s).")
            ctx = mp.get_context("fork") if hasattr(os, "fork") else mp.get_context("spawn")
            with ctx.Pool(
                processes=threads,
                initializer=init_worker,
                initargs=(markers, k, window, step),
            ) as pool:
                for chrom, nwin, valid, skipped, hits, rows in pool.imap(scan_chrom_worker, parse_fasta(genome), chunksize=1):
                    # imap preserves FASTA order, so output order remains stable.
                    chrom_count += 1
                    # bases cannot be recovered from result directly; count from rows not ideal.
                    # Return rows only after worker completes; log concise per-contig stats.
                    out.write(rows)
                    total_windows += nwin
                    total_valid_kmers += valid
                    total_skipped_n += skipped
                    total_hits += hits
                    # total_bases is accumulated in a second lightweight pass below for clean logging.
                    log("Info", f"  {chrom}: windows={nwin}, valid_kmers={valid}, skipped_N/non-ACGT={skipped}, marker_hits={hits}")

            # Multiprocessing consumes FASTA records in workers; compute total bases with a cheap pass
            # to keep the output statistics identical to single-worker mode.
            for _chrom, seq in parse_fasta(genome):
                total_bases += len(seq)

    return chrom_count, total_bases, total_windows, total_valid_kmers, total_skipped_n, total_hits


def main() -> int:
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--kmer-size", required=True)
    ap.add_argument("--window", type=int, required=True)
    ap.add_argument("--step", type=int, required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    try:
        markers, k, raw_count, indexed_count = load_markers(args.markers, args.kmer_size)
        if args.window < k:
            raise ValueError(f"Window size ({args.window}) must be >= k-mer size ({k})")

        log("Info", "2-bit integer marker index built successfully.")
        log("Info", f"  K-mer size      : {k}")
        log("Info", f"  Marker sequences: {raw_count}")
        log("Info", f"  Indexed k-mers   : {indexed_count} (forward + reverse-complement, deduplicated)")
        log("Info", "  Scanner backend : rolling 2-bit integer exact matching")

        chrom_count, total_bases, total_windows, total_valid_kmers, total_skipped_n, total_hits = run_scan(
            args.genome, markers, k, args.window, args.step, args.threads, args.out
        )

        log("Info", "SEXmer mapping complete.")
        log("Info", f"  Chromosomes/contigs : {chrom_count}")
        log("Info", f"  Genome bases        : {total_bases}")
        log("Info", f"  Windows written     : {total_windows}")
        log("Info", f"  Valid genome k-mers : {total_valid_kmers}")
        log("Info", f"  Skipped N/non-ACGT  : {total_skipped_n}")
        log("Info", f"  Marker hits         : {total_hits}")
        return 0

    except Exception as exc:
        log("Error", str(exc))
        return 1


if __name__ == "__main__":
    sys.exit(main())

PY
chmod +x "$SCANNER"

OUT_TSV="${PREFIX}.windows.tsv"

info "SEXmer-mapping starting"
info "Parameters : kmer-size=${KMER_SIZE}, window=${WINDOW}, step=${STEP}, count-mode=${COUNT_MODE}"
info "Settings   : mem=${MEM}, threads=${THREADS}"
info "Genome     : ${GENOME}"
info "Markers    : ${MARKERS}"
info "Output     : ${OUT_TSV}"
info "Temp dir   : ${MAPPING_TMPDIR}"
info "Coordinate : 0-based, half-open windows"
info "Strand mode: forward + reverse-complement marker lookup"
info "Normalize  : density=hits/valid_kmers; hits_per_10kb=density*10000"

info "Running coordinate-aware k-mer window scanner..."
python3 "$SCANNER" \
    --genome "$GENOME" \
    --markers "$MARKERS" \
    --kmer-size "$KMER_SIZE" \
    --window "$WINDOW" \
    --step "$STEP" \
    --threads "$THREADS" \
    --out "$OUT_TSV"

[[ -f "$OUT_TSV" ]] || { error "Scanner did not produce output file: $OUT_TSV"; exit 1; }

output "Window table written to: ${OUT_TSV}"
info "SEXmer-mapping complete."
