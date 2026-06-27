#!/usr/bin/env bash
# SEXmer-reads_v3.sh - Extract reads containing SEXmer marker k-mers using pyahocorasick.
# Author: Dede Kurniawan

set -euo pipefail
export LC_ALL=C

# defaults
MARKERS=""
PREFIX="output"
KMER_SIZE="auto"
MIN_HIT=1
TMPDIR_BASE="$(pwd)"
READS=()

# log helpers
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*"  >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"   >&2; }

usage() {
    cat >&2 <<EOF_USAGE

SEXmer-reads_v3.sh - Extract reads containing exact SEXmer marker k-mers using pyahocorasick.

Usage: SEXmer-reads_v3.sh --markers <markers.fa> [OPTIONS] <reads...>

Mandatory:
  --markers            FASTA file containing marker k-mers
  <reads>              One or more FASTQ files (.fq, .fastq, .fq.gz, .fastq.gz)

Optional:
  --prefix             Output filename prefix                          [default: ${PREFIX}]
  --hit                Minimum exact marker k-mer hits per read         [default: ${MIN_HIT}]
  -k, --kmer-size      K-mer size; auto-detected from markers if unset  [default: auto]
  --tmpdir             Parent directory for the temporary work folder   [default: current dir]
  -h, --help           Show this help and exit

Output:
  Single/long-read mode : <prefix>.sexmer.fq.gz
  Paired-end mode       : <prefix>.sexmer_1.fq.gz and <prefix>.sexmer_2.fq.gz

Notes:
  - Matching is exact k-mer matching only.
  - Reverse-complement marker k-mers are always included.
  - For paired-end reads, both mates are written if either mate reaches --hit.
  - If exactly two read files are provided, paired-end mode is detected from filename
    pattern and/or read ID agreement. Otherwise files are treated as independent
    single/long-read FASTQ files.

EOF_USAGE
    exit 1
}

[[ $# -eq 0 ]] && usage

# argument parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        --markers)          MARKERS="$2";      shift 2 ;;
        --prefix)           PREFIX="$2";       shift 2 ;;
        --hit)              MIN_HIT="$2";      shift 2 ;;
        -k|--kmer-size)     KMER_SIZE="$2";    shift 2 ;;
        --tmpdir)           TMPDIR_BASE="$2";  shift 2 ;;
        -h|--help)          usage ;;
        -*) error "Unknown option '$1'"; usage ;;
        *)  READS+=("$1"); shift ;;
    esac
done

[[ -z "$MARKERS" ]] && { error "--markers is required."; usage; }
[[ ${#READS[@]} -eq 0 ]] && { error "At least one input read file is required."; usage; }

[[ -r "$MARKERS" ]] || { error "Cannot read marker FASTA file: $MARKERS"; exit 1; }

[[ "$MIN_HIT" =~ ^[1-9][0-9]*$ ]] || {
    error "--hit must be a positive integer."; exit 1; }

if [[ "$KMER_SIZE" != "auto" ]]; then
    [[ "$KMER_SIZE" =~ ^[1-9][0-9]*$ ]] && [[ "$KMER_SIZE" -le 63 ]] || {
        error "--kmer-size must be 'auto' or an integer between 1 and 63."; exit 1; }
fi

[[ -d "$TMPDIR_BASE" ]] || { error "Temporary parent directory does not exist: ${TMPDIR_BASE}"; exit 1; }
[[ -w "$TMPDIR_BASE" ]] || { error "Temporary parent directory is not writable: ${TMPDIR_BASE}"; exit 1; }

OUT_DIR="$(dirname "$PREFIX")"
[[ "$OUT_DIR" == "." ]] && OUT_DIR="$(pwd)"
[[ -d "$OUT_DIR" ]] || { error "Output directory does not exist: ${OUT_DIR}"; exit 1; }
[[ -w "$OUT_DIR" ]] || { error "Output directory is not writable: ${OUT_DIR}"; exit 1; }

for f in "${READS[@]}"; do
    [[ -r "$f" ]] || { error "Cannot read input file: $f"; exit 1; }
done

command -v python3 &>/dev/null || { error "python3 not found on PATH."; exit 1; }
command -v gzip &>/dev/null || { error "gzip not found on PATH."; exit 1; }

if ! python3 - <<'PY_CHECK'
try:
    import ahocorasick  # pyahocorasick imports as ahocorasick
except Exception:
    raise SystemExit(1)
PY_CHECK
then
    error "Python package 'pyahocorasick' is required for SEXmer-reads_v3.sh."
    error "Install it with: pip install pyahocorasick  OR  conda install -c conda-forge pyahocorasick"
    exit 1
fi

READS_TMPDIR="${TMPDIR_BASE}/sexmer_reads_tmp_$$"
mkdir -p "$READS_TMPDIR"
cleanup() { rm -rf "$READS_TMPDIR"; }
trap cleanup EXIT

SCANNER="${READS_TMPDIR}/sexmer_reads_scanner.py"
cat > "$SCANNER" <<'PY'
#!/usr/bin/env python3
import argparse
import gzip
import ahocorasick
import os
import re
import sys
from typing import Dict, List, Sequence, Tuple

DNA = set("ACGT")
RC_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")


def log(kind: str, msg: str) -> None:
    print(f"[{kind}] {msg}", file=sys.stderr)


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def open_gzip_out(path: str):
    return gzip.open(path, "wt", compresslevel=1)


def revcomp(seq: str) -> str:
    return seq.translate(RC_TRANS)[::-1].upper()


def parse_fasta(path: str) -> List[str]:
    seqs: List[str] = []
    chunks: List[str] = []
    with open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chunks:
                    seqs.append("".join(chunks).upper())
                    chunks = []
            else:
                chunks.append(line)
        if chunks:
            seqs.append("".join(chunks).upper())
    return seqs


def build_marker_automaton(path: str, k_arg: str) -> Tuple[ahocorasick.Automaton, int, int, int]:
    seqs = parse_fasta(path)
    if not seqs:
        raise ValueError(f"No marker sequence found in FASTA: {path}")

    if k_arg == "auto":
        k = len(seqs[0])
    else:
        k = int(k_arg)

    if k < 1 or k > 63:
        raise ValueError("k-mer size must be between 1 and 63")

    # Deduplicate forward and reverse-complement marker k-mers before building
    # the Aho-Corasick automaton. This keeps the index smaller and avoids
    # double-counting identical marker sequences.
    marker_set = set()
    raw_count = 0
    for seq in seqs:
        if len(seq) != k:
            raise ValueError(
                f"Marker length mismatch: expected {k}, found {len(seq)} for sequence '{seq[:40]}...'"
            )
        if any(base not in DNA for base in seq):
            raise ValueError(f"Marker contains non-ACGT base: '{seq[:40]}...'")
        raw_count += 1
        marker_set.add(seq)
        marker_set.add(revcomp(seq))

    automaton = ahocorasick.Automaton()
    for marker in marker_set:
        # Store the marker itself as the value. The value is not needed for
        # extraction, but it keeps debugging simple if we later expose hits.
        automaton.add_word(marker, marker)
    automaton.make_automaton()

    return automaton, k, raw_count, len(marker_set)


def fastq_iter(path: str):
    with open_text(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline()
            p = fh.readline()
            q = fh.readline()
            if not q:
                raise ValueError(f"Truncated FASTQ record in {path}")
            yield h, s, p, q


def normalize_read_id(header: str) -> str:
    h = header.strip()
    if h.startswith("@"):
        h = h[1:]
    h = h.split()[0]
    h = re.sub(r"/[12]$", "", h)
    return h


def filename_pair_score(path1: str, path2: str) -> bool:
    b1 = os.path.basename(path1)
    b2 = os.path.basename(path2)

    suffixes = [".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"]
    stem1, stem2 = b1, b2
    for suf in suffixes:
        if stem1.endswith(suf):
            stem1 = stem1[: -len(suf)]
        if stem2.endswith(suf):
            stem2 = stem2[: -len(suf)]

    patterns = [
        (r"(.+)_1$", r"\1_2"),
        (r"(.+)_R1$", r"\1_R2"),
        (r"(.+)_R1_001$", r"\1_R2_001"),
        (r"(.+)\.1$", r"\1.2"),
        (r"(.+)-1$", r"\1-2"),
    ]

    for pat, repl in patterns:
        m = re.match(pat, stem1)
        if m and re.sub(pat, repl, stem1) == stem2:
            return True
    return False


def header_pair_score(path1: str, path2: str, max_records: int = 1000) -> Tuple[bool, int, int]:
    total = 0
    same = 0
    it1 = fastq_iter(path1)
    it2 = fastq_iter(path2)

    for _ in range(max_records):
        try:
            r1 = next(it1)
            r2 = next(it2)
        except StopIteration:
            break
        total += 1
        if normalize_read_id(r1[0]) == normalize_read_id(r2[0]):
            same += 1

    if total == 0:
        return False, same, total
    return (same / total) >= 0.90, same, total


def count_hits(seq: str, automaton: ahocorasick.Automaton, k: int, min_hit: int) -> int:
    seq = seq.strip().upper()
    if len(seq) < k:
        return 0
    hits = 0
    # pyahocorasick performs exact multi-pattern search in C and reports
    # matches while streaming over the sequence. We stop as soon as min_hit is
    # reached because extraction only needs a pass/fail decision.
    for _end, _marker in automaton.iter(seq):
        hits += 1
        if hits >= min_hit:
            return hits
    return hits


def write_record(out, rec) -> None:
    out.write(rec[0])
    out.write(rec[1])
    out.write(rec[2])
    out.write(rec[3])


def scan_single(automaton: ahocorasick.Automaton, k: int, min_hit: int, read_files: Sequence[str], out_path: str) -> Dict[str, int]:
    total_reads = 0
    kept_reads = 0
    with open_gzip_out(out_path) as out:
        for path in read_files:
            file_total = 0
            file_kept = 0
            for rec in fastq_iter(path):
                file_total += 1
                total_reads += 1
                if count_hits(rec[1], automaton, k, min_hit) >= min_hit:
                    write_record(out, rec)
                    file_kept += 1
                    kept_reads += 1
            log("Info", f"  {os.path.basename(path)}: {file_kept} of {file_total} reads extracted.")
    return {"total_reads": total_reads, "kept_reads": kept_reads}


def scan_single_one(markers_path: str, k_arg: str, min_hit: int, read_file: str, out_path: str) -> None:
    automaton, k, raw_count, indexed_count = build_marker_automaton(markers_path, k_arg)
    stats = scan_single(automaton, k, min_hit, [read_file], out_path)
    print(f"k\t{k}")
    print(f"raw_markers\t{raw_count}")
    print(f"indexed_markers\t{indexed_count}")
    print(f"total_reads\t{stats['total_reads']}")
    print(f"kept_reads\t{stats['kept_reads']}")


def scan_paired(automaton: ahocorasick.Automaton, k: int, min_hit: int, r1_path: str, r2_path: str, out1_path: str, out2_path: str) -> Dict[str, int]:
    total_pairs = 0
    kept_pairs = 0
    with open_gzip_out(out1_path) as out1, open_gzip_out(out2_path) as out2:
        it1 = fastq_iter(r1_path)
        it2 = fastq_iter(r2_path)
        while True:
            try:
                rec1 = next(it1)
            except StopIteration:
                try:
                    next(it2)
                    raise ValueError("R2 has more records than R1")
                except StopIteration:
                    break
            try:
                rec2 = next(it2)
            except StopIteration:
                raise ValueError("R1 has more records than R2")

            total_pairs += 1
            hit1 = count_hits(rec1[1], automaton, k, min_hit)
            hit2 = count_hits(rec2[1], automaton, k, min_hit)
            if hit1 >= min_hit or hit2 >= min_hit:
                write_record(out1, rec1)
                write_record(out2, rec2)
                kept_pairs += 1
    return {"total_pairs": total_pairs, "kept_pairs": kept_pairs}


def detect_mode(read_files: Sequence[str]) -> None:
    if len(read_files) == 1:
        print("single")
        print("filename_pair\tno")
        print("header_pair\tno")
        print("header_same\t0")
        print("header_total\t0")
        return
    if len(read_files) != 2:
        print("single")
        print("filename_pair\tno")
        print("header_pair\tno")
        print("header_same\t0")
        print("header_total\t0")
        return

    fname_pair = filename_pair_score(read_files[0], read_files[1])
    h_pair, same, total = header_pair_score(read_files[0], read_files[1])
    is_pair = fname_pair or h_pair
    print("paired" if is_pair else "single")
    print(f"filename_pair\t{'yes' if fname_pair else 'no'}")
    print(f"header_pair\t{'yes' if h_pair else 'no'}")
    print(f"header_same\t{same}")
    print(f"header_total\t{total}")


def main() -> int:
    ap = argparse.ArgumentParser(add_help=False)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_detect = sub.add_parser("detect", add_help=False)
    p_detect.add_argument("reads", nargs="+")

    p_single = sub.add_parser("single", add_help=False)
    p_single.add_argument("--markers", required=True)
    p_single.add_argument("--kmer-size", required=True)
    p_single.add_argument("--hit", type=int, required=True)
    p_single.add_argument("--out", required=True)
    p_single.add_argument("reads", nargs="+")

    p_single_one = sub.add_parser("single-one", add_help=False)
    p_single_one.add_argument("--markers", required=True)
    p_single_one.add_argument("--kmer-size", required=True)
    p_single_one.add_argument("--hit", type=int, required=True)
    p_single_one.add_argument("--out", required=True)
    p_single_one.add_argument("read")

    p_paired = sub.add_parser("paired", add_help=False)
    p_paired.add_argument("--markers", required=True)
    p_paired.add_argument("--kmer-size", required=True)
    p_paired.add_argument("--hit", type=int, required=True)
    p_paired.add_argument("--r1", required=True)
    p_paired.add_argument("--r2", required=True)
    p_paired.add_argument("--out1", required=True)
    p_paired.add_argument("--out2", required=True)

    p_marker = sub.add_parser("marker-info", add_help=False)
    p_marker.add_argument("--markers", required=True)
    p_marker.add_argument("--kmer-size", required=True)

    args = ap.parse_args()

    try:
        if args.cmd == "detect":
            detect_mode(args.reads)
            return 0

        if args.cmd == "marker-info":
            automaton, k, raw_count, indexed_count = build_marker_automaton(args.markers, args.kmer_size)
            print(f"k\t{k}")
            print(f"raw_markers\t{raw_count}")
            print(f"indexed_markers\t{indexed_count}")
            return 0

        if args.cmd == "single-one":
            scan_single_one(args.markers, args.kmer_size, args.hit, args.read, args.out)
            return 0

        automaton, k, raw_count, indexed_count = build_marker_automaton(args.markers, args.kmer_size)

        if args.cmd == "single":
            stats = scan_single(automaton, k, args.hit, args.reads, args.out)
            print(f"k\t{k}")
            print(f"raw_markers\t{raw_count}")
            print(f"indexed_markers\t{indexed_count}")
            print(f"total_reads\t{stats['total_reads']}")
            print(f"kept_reads\t{stats['kept_reads']}")
            return 0

        if args.cmd == "paired":
            stats = scan_paired(automaton, k, args.hit, args.r1, args.r2, args.out1, args.out2)
            print(f"k\t{k}")
            print(f"raw_markers\t{raw_count}")
            print(f"indexed_markers\t{indexed_count}")
            print(f"total_pairs\t{stats['total_pairs']}")
            print(f"kept_pairs\t{stats['kept_pairs']}")
            return 0

    except Exception as exc:
        log("Error", str(exc))
        return 1

    log("Error", "Unknown scanner command")
    return 1


if __name__ == "__main__":
    sys.exit(main())
PY
chmod +x "$SCANNER"

info "SEXmer-reads v3 starting"
info "Parameters : hit=${MIN_HIT}, kmer-size=${KMER_SIZE}"
info "Marker file: ${MARKERS}"
info "Input reads: ${READS[*]}"
info "Temp dir   : ${READS_TMPDIR}"

# Load marker information once for validation and logging.
MARKER_INFO="${READS_TMPDIR}/marker_info.txt"
python3 "$SCANNER" marker-info \
    --markers "$MARKERS" \
    --kmer-size "$KMER_SIZE" \
    > "$MARKER_INFO"

DETECTED_K=$(awk -F'\t' '$1 == "k" { print $2 }' "$MARKER_INFO")
RAW_MARKERS=$(awk -F'\t' '$1 == "raw_markers" { print $2 }' "$MARKER_INFO")
INDEXED_MARKERS=$(awk -F'\t' '$1 == "indexed_markers" { print $2 }' "$MARKER_INFO")

info "pyahocorasick marker index built successfully."
info "  K-mer size      : ${DETECTED_K}"
info "  Marker sequences: ${RAW_MARKERS}"
info "  Indexed kmers   : ${INDEXED_MARKERS} (forward + reverse-complement, deduplicated)"
info "  Scanner backend : pyahocorasick"

# Detect read mode.
DETECT_INFO="${READS_TMPDIR}/detect_mode.txt"
python3 "$SCANNER" detect "${READS[@]}" > "$DETECT_INFO"
MODE=$(head -n 1 "$DETECT_INFO")
FILENAME_PAIR=$(awk -F'\t' '$1 == "filename_pair" { print $2 }' "$DETECT_INFO")
HEADER_PAIR=$(awk -F'\t' '$1 == "header_pair" { print $2 }' "$DETECT_INFO")
HEADER_SAME=$(awk -F'\t' '$1 == "header_same" { print $2 }' "$DETECT_INFO")
HEADER_TOTAL=$(awk -F'\t' '$1 == "header_total" { print $2 }' "$DETECT_INFO")

if [[ "$MODE" == "paired" ]]; then
    info "Read mode detected: paired-end"
    info "  Filename pair pattern: ${FILENAME_PAIR}"
    info "  Header ID agreement  : ${HEADER_SAME}/${HEADER_TOTAL} sampled read pairs"

    OUT1="${PREFIX}.sexmer_1.fq.gz"
    OUT2="${PREFIX}.sexmer_2.fq.gz"

    info "Extracting paired reads..."
    RUN_INFO="${READS_TMPDIR}/paired_run_info.txt"
    python3 "$SCANNER" paired \
        --markers "$MARKERS" \
        --kmer-size "$KMER_SIZE" \
        --hit "$MIN_HIT" \
        --r1 "${READS[0]}" \
        --r2 "${READS[1]}" \
        --out1 "$OUT1" \
        --out2 "$OUT2" \
        > "$RUN_INFO"

    TOTAL_PAIRS=$(awk -F'\t' '$1 == "total_pairs" { print $2 }' "$RUN_INFO")
    KEPT_PAIRS=$(awk -F'\t' '$1 == "kept_pairs" { print $2 }' "$RUN_INFO")

    output "R1 reads written to: ${OUT1}"
    output "R2 reads written to: ${OUT2}"
    info "SEXmer reads complete."
    info "  Mode           : paired-end"
    info "  Pair rule      : write both mates if either mate has >= ${MIN_HIT} exact marker hit(s)"
    info "  Input pairs    : ${TOTAL_PAIRS}"
    info "  Extracted pairs: ${KEPT_PAIRS}"
else
    if [[ ${#READS[@]} -eq 2 ]]; then
        info "Read mode detected: single/long-read"
        info "  Filename pair pattern: ${FILENAME_PAIR}"
        info "  Header ID agreement  : ${HEADER_SAME}/${HEADER_TOTAL} sampled records"
    else
        info "Read mode detected: single/long-read"
        info "  Input files: ${#READS[@]}"
    fi

    OUT="${PREFIX}.sexmer.fq.gz"
    info "Extracting single/long reads..."

    if [[ ${#READS[@]} -gt 1 ]]; then
        info "Running one streaming worker per input file..."
        PIDS=()
        TMP_OUTS=()
        TMP_STATS=()
        for i in "${!READS[@]}"; do
            tmp_out="${READS_TMPDIR}/single_${i}.fq.gz"
            tmp_stat="${READS_TMPDIR}/single_${i}.stats"
            TMP_OUTS+=("$tmp_out")
            TMP_STATS+=("$tmp_stat")
            python3 "$SCANNER" single-one \
                --markers "$MARKERS" \
                --kmer-size "$KMER_SIZE" \
                --hit "$MIN_HIT" \
                --out "$tmp_out" \
                "${READS[$i]}" \
                > "$tmp_stat" &
            PIDS+=("$!")
        done

        for pid in "${PIDS[@]}"; do
            wait "$pid"
        done

        cat "${TMP_OUTS[@]}" > "$OUT"

        TOTAL_READS=0
        KEPT_READS=0
        for stat in "${TMP_STATS[@]}"; do
            t=$(awk -F'\t' '$1 == "total_reads" { print $2 }' "$stat")
            k=$(awk -F'\t' '$1 == "kept_reads" { print $2 }' "$stat")
            TOTAL_READS=$(( TOTAL_READS + t ))
            KEPT_READS=$(( KEPT_READS + k ))
        done
    else
        RUN_INFO="${READS_TMPDIR}/single_run_info.txt"
        python3 "$SCANNER" single \
            --markers "$MARKERS" \
            --kmer-size "$KMER_SIZE" \
            --hit "$MIN_HIT" \
            --out "$OUT" \
            "${READS[@]}" \
            > "$RUN_INFO"
        TOTAL_READS=$(awk -F'\t' '$1 == "total_reads" { print $2 }' "$RUN_INFO")
        KEPT_READS=$(awk -F'\t' '$1 == "kept_reads" { print $2 }' "$RUN_INFO")
    fi

    output "Reads written to: ${OUT}"
    info "SEXmer reads complete."
    info "  Mode           : single/long-read"
    info "  Read rule      : write read if it has >= ${MIN_HIT} exact marker hit(s)"
    info "  Input reads    : ${TOTAL_READS}"
    info "  Extracted reads: ${KEPT_READS}"
fi
