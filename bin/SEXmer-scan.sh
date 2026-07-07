#!/usr/bin/env bash
# SEXmer scan - Scan all k-mer sequences and classify them as MSK, FSK, MBK, FBK, or neutral.
# Author: Dede Kurniawan 

set -euo pipefail
export LC_ALL=C

MIN_COUNT=3
MAX_COUNT=1000
FOLD_THRESHOLD=5
PREFIX="output"
OUTDIR="."
MEM="8G"
THREADS=8
MALE_INPUT=""
FEMALE_INPUT=""
NEUTRAL_MAX=100000
TMPDIR_BASE="$(pwd)"
SEED=42
PLOT_ENABLED=1
PLOT_FORMAT="png"

# log helpers
info()    { echo "[Info] $*"    >&2; }
output()  { echo "[Output] $*" >&2; }
warn()    { echo "[Warning] $*" >&2; }
error()   { echo "[Error] $*"  >&2; }

usage() {
    cat <<EOF

SEXmer scan - Scan all k-mer sequences and classify them as MSK, FSK, MBK, FBK, or neutral.

Usage: SEXmer scan -m <male_files> -f <female_files> [OPTIONS]

Mandatory:
  -m, --male           Specify male dump files (separated by commas)
  -f, --female         Specify female dump files (separated by commas)

Optional:
  --prefix             The prefix used on generated files              [default: output]
  -o, --outdir         Specify output directory name                   [default: current dir]
  --mem                Specify max MEM for this task (e.g. 8G)         [default: ${MEM}]
  -t, --threads        Specify CPU threads for this task               [default: ${THREADS}]
  --neutral-max        Maximum neutral k-mers to retain, 0=keep all    [default: ${NEUTRAL_MAX}]
  --tmpdir             Parent directory for the temporary work folder  [default: current dir]
  --min-count          Minimum k-mer count to retain                   [default: ${MIN_COUNT}]
  --max-count          Maximum pooled k-mer count within one sex       [default: ${MAX_COUNT}]
  --fold-threshold     Specify fold-change cutoff for MBK/FBK          [default: ${FOLD_THRESHOLD}]
  --seed               Specify random seed for neutral k-mer sampling  [default: ${SEED}]
  --no-plot            Do not generate any visualization
  --plot-format        Specify plot format: svg, png, or pdf           [default: ${PLOT_FORMAT}]
  -h, --help           Show this help message and exit

Categories:
  MSK      Male-specific k-mer
  FSK      Female-specific k-mer
  MBK      Male-biased k-mer
  FBK      Female-biased k-mer
  neutral  K-mer without sex-specific or sex-biased signal

EOF
}

if [[ $# -eq 0 ]]; then
    usage >&2
    exit 1
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--male)           MALE_INPUT="$2";     shift 2 ;;
        -f|--female)         FEMALE_INPUT="$2";   shift 2 ;;
        --prefix)            PREFIX="$2";         shift 2 ;;
        -o|--outdir)         OUTDIR="$2";         shift 2 ;;
        --mem)               MEM="$2";            shift 2 ;;
        -t|--threads)        THREADS="$2";        shift 2 ;;
        --neutral-max)       NEUTRAL_MAX="$2";    shift 2 ;;
        --tmpdir)            TMPDIR_BASE="$2";    shift 2 ;;
        --min-count)         MIN_COUNT="$2";      shift 2 ;;
        --max-count)         MAX_COUNT="$2";      shift 2 ;;
        --fold-threshold)    FOLD_THRESHOLD="$2"; shift 2 ;;
        --seed)              SEED="$2";           shift 2 ;;
        --no-plot)           PLOT_ENABLED=0;      shift ;;
        --plot-format)       PLOT_FORMAT="$2";    shift 2 ;;
        -h|--help)           usage; exit 0 ;;
        *) error "Unknown option '$1'"; usage >&2; exit 1 ;;
    esac
done

[[ -z "$MALE_INPUT" ]]   && { error "-m/--male files not specified.";   usage >&2; exit 1; }
[[ -z "$FEMALE_INPUT" ]] && { error "-f/--female files not specified."; usage >&2; exit 1; }

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

case "$PLOT_FORMAT" in
    svg|png|pdf) ;;
    *) error "--plot-format must be one of: svg, png, pdf."; exit 1 ;;
esac

[[ -d "$TMPDIR_BASE" ]] || { error "Temporary parent directory does not exist: $TMPDIR_BASE"; exit 1; }
[[ -w "$TMPDIR_BASE" ]] || { error "Temporary parent directory is not writable: $TMPDIR_BASE"; exit 1; }

mkdir -p "$OUTDIR" || { error "Failed to create output directory: $OUTDIR"; exit 1; }
[[ -d "$OUTDIR" ]] || { error "Output path is not a directory: $OUTDIR"; exit 1; }
[[ -w "$OUTDIR" ]] || { error "Output directory is not writable: $OUTDIR"; exit 1; }

if [[ "$OUTDIR" == "." ]]; then
    OUTPUT_PREFIX="$PREFIX"
else
    OUTPUT_PREFIX="${OUTDIR%/}/${PREFIX}"
fi

IFS=',' read -ra MALE_FILES   <<< "$MALE_INPUT"
IFS=',' read -ra FEMALE_FILES <<< "$FEMALE_INPUT"

for f in "${MALE_FILES[@]}" "${FEMALE_FILES[@]}"; do
    [[ -r "$f" ]] || { error "Cannot read input file: $f"; exit 1; }
done

if ! sort --version 2>&1 | grep -q 'GNU'; then
    error "GNU sort is required but not found. On macOS: brew install coreutils"
    exit 1
fi

if [[ "$PLOT_ENABLED" -eq 1 ]]; then
    command -v python3 &>/dev/null || { error "python3 not found on PATH. It is required for default SVG plotting; use --no-plot to skip plotting."; exit 1; }
    python3 - <<'PYCHECK' || { error "Python plotting dependencies missing. Required modules: pandas, numpy, matplotlib, scipy. Use --no-plot to skip plotting."; exit 1; }
import importlib.util
import sys
missing = [m for m in ("pandas", "numpy", "matplotlib", "scipy") if importlib.util.find_spec(m) is None]
if missing:
    print("Missing Python module(s): " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
PYCHECK
fi

# set up temp directory
SCAN_TMPDIR="${TMPDIR_BASE}/sexmer_scan_tmp_$$"
mkdir -p "$SCAN_TMPDIR"
cleanup() { rm -rf "$SCAN_TMPDIR"; }
trap cleanup EXIT

if [[ "$PLOT_ENABLED" -eq 1 ]]; then
    PLOTTER="${SCAN_TMPDIR}/sexmer_scan_plot.py"
    cat > "$PLOTTER" <<'PYPLOT'
#!/usr/bin/env python3

import argparse
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

# defaults
COLORS = {
    'MSK':     '#0e616d',
    'FSK':     '#f76e68',
    'MBK':     '#fec44a',
    'FBK':     '#1497a5',
    'neutral': '#dddddd',
}

CATEGORY_ORDER = ['FSK', 'MSK', 'FBK', 'MBK', 'neutral']

def info(msg):   print(f"[Info] {msg}",    file=sys.stderr)
def output(msg): print(f"[Output] {msg}",  file=sys.stderr)
def warn(msg):   print(f"[Warning] {msg}", file=sys.stderr)
def error(msg):  print(f"[Error] {msg}",   file=sys.stderr)

def usage():
    print("""
SEXmer-plot.py - Generate publication-ready figures from SEXmer scan output.

Usage: SEXmer-plot.py -i <input.kmers.tsv> [OPTIONS]

Mandatory:
  -i, --input          Input TSV file from SEXmer scan (*.kmers.tsv)

Optional:
  --prefix             Output filename prefix            [default: input basename]
  --format             Output format: svg, png, or pdf    [default: svg]
  -h, --help           Show this help and exit

""", file=sys.stderr)
    sys.exit(1)

def fmt_num(n):
    if n >= 1_000_000:
        return f"{n/1_000_000:.2f}M"
    if n >= 1_000:
        return f"{n/1_000:.1f}K"
    return str(int(n))

def clean_axes(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

def load_data(path):
    info("Reading input data...")
    try:
        df = pd.read_csv(
            path, sep='\t',
            usecols=['male_count', 'female_count', 'category'],
            dtype={'male_count': np.int32, 'female_count': np.int32, 'category': 'category'},
        )
    except Exception as e:
        error(f"Failed to read input file: {e}")
        sys.exit(1)

    if df.empty:
        error("Input file is empty or has no valid rows.")
        sys.exit(1)

    info(f"Loaded {len(df):,} k-mers.")
    return df

# scatter_panel
def scatter_panel(ax, mc, fc, cats):
    STYLE = {
        'MSK':     dict(s=35, alpha=0.80, zorder=5),
        'FSK':     dict(s=35, alpha=0.80, zorder=5),
        'MBK':     dict(s=35, alpha=0.80, zorder=4),
        'FBK':     dict(s=35, alpha=0.80, zorder=4),
        'neutral': dict(s=35, alpha=0.35, zorder=3),
    }

    for cat in CATEGORY_ORDER:
        mask = cats == cat
        if not mask.any():
            continue
        label = cat if cat != 'neutral' else None
        ax.scatter(mc[mask], fc[mask], c=COLORS[cat], edgecolors='none',
                   label=label, **STYLE[cat])

    ax.set_xlabel('k-mer counts in males',   fontsize=13)
    ax.set_ylabel('k-mer counts in females', fontsize=13)

    xmax  = float(np.nanmax(mc)) if len(mc) else 10
    ymax  = float(np.nanmax(fc)) if len(fc) else 10
    pad_x = xmax * 0.03
    pad_y = ymax * 0.03
    ax.set_xlim(-pad_x, xmax * 1.05)
    ax.set_ylim(-pad_y, ymax * 1.05)
    ax.set_aspect('equal')

    ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        order = ['FSK', 'MSK', 'FBK', 'MBK']
        pairs = [(h, l) for cat in order for h, l in zip(handles, labels) if l == cat]
        if pairs:
            oh, ol = zip(*pairs)
            ax.legend(oh, ol, loc='upper right', framealpha=0.95,
                      fontsize=10, title='Category', title_fontsize=11,
                      edgecolor='gray', fancybox=False)

    clean_axes(ax)
    ax.tick_params(axis='both', labelsize=10)

# bar_panel
def bar_panel(ax, df):
    counts = [len(df[df['category'] == c]) for c in ('FSK', 'MSK')]
    bars   = ax.bar(['FSK', 'MSK'], counts, color=[COLORS['FSK'], COLORS['MSK']],
                    width=0.6, edgecolor='black', linewidth=1.2)

    for bar, cnt in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2., bar.get_height(),
                fmt_num(cnt), ha='center', va='bottom',
                fontsize=12, fontweight='bold')

    ax.set_ylabel('k-mer counts', fontsize=13)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: fmt_num(v)))
    ax.yaxis.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_ylim(0, (max(counts) if counts else 1) * 1.25)

    clean_axes(ax)
    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=10)

# create_sexplot
def create_sexplot(df, out_path):
    info("Generating sex-specific k-mer scatter + bar plot...")

    mc   = df['male_count'].values.astype(np.float32)
    fc   = df['female_count'].values.astype(np.float32)
    cats = df['category'].values

    fig = plt.figure(figsize=(11, 7))
    gs  = GridSpec(1, 2, figure=fig, width_ratios=[2, 1], wspace=0.2)

    scatter_panel(fig.add_subplot(gs[0, 0]), mc, fc, cats)
    bar_panel(fig.add_subplot(gs[0, 1]), df)

    fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    output(f"Sex plot written to: {out_path}")

# add_reflines
def add_reflines(ax, lim):
    xv = np.array([0.0, lim])
    ax.plot(xv, xv * 1.0, color='red',   lw=1.5, ls='--', zorder=5, label='1:1')
    ax.plot(xv, xv * 2.0, color='blue',  lw=1.5, ls='--', zorder=5, label='2:1 (female-biased)')
    ax.plot(xv, xv * 0.5, color='green', lw=1.5, ls='--', zorder=5, label='1:2 (male-biased)')

def format_abundance_ax(ax, lim):
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel('k-mer counts in males',   fontsize=11)
    ax.set_ylabel('k-mer counts in females', fontsize=11)
    ax.set_aspect('equal')
    clean_axes(ax)

def smooth_density(ax, x, y, lim, bins=400, sigma=3):
    h, _, _   = np.histogram2d(x, y, bins=[bins, bins], range=[[0, lim], [0, lim]])
    h_smooth  = gaussian_filter(h, sigma=sigma)
    h_display = np.sqrt(h_smooth)

    blues     = matplotlib.colormaps.get_cmap('Blues')
    cmap_list = [(1, 1, 1)] + [blues(v) for v in np.linspace(0.15, 1.0, 255)]
    cmap      = mcolors.LinearSegmentedColormap.from_list('smoothScatter', cmap_list)

    ax.imshow(h_display.T, origin='lower', aspect='equal',
              extent=[0, lim, 0, lim], cmap=cmap, interpolation='bilinear')

# create_abundance_plot
def create_abundance_plot(df, out_path):
    info("Generating abundance scatter + density plot...")

    x = df['male_count'].values.astype(np.float64)
    y = df['female_count'].values.astype(np.float64)

    lim = float(np.ceil(np.percentile(np.concatenate([x, y]), 99)))
    if lim <= 0:
        lim = 10.0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), constrained_layout=True)
    fig.patch.set_facecolor('white')

    ax1.scatter(x, y, s=1.5, c='black', alpha=0.35, linewidths=0, rasterized=True)
    add_reflines(ax1, lim)
    format_abundance_ax(ax1, lim)

    smooth_density(ax2, x, y, lim)
    add_reflines(ax2, lim)
    format_abundance_ax(ax2, lim)

    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    output(f"Abundance plot written to: {out_path}")

# argument parsing
def parse_args():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-i', '--input',  default=None)
    parser.add_argument('--prefix',       default=None)
    parser.add_argument('--format',       choices=['svg', 'png', 'pdf'], default='svg')
    parser.add_argument('-h', '--help',   action='store_true')
    return parser.parse_args()

def main():
    args = parse_args()

    if args.help or args.input is None:
        usage()

    if not os.path.isfile(args.input):
        error(f"Input file not found: {args.input}")
        sys.exit(1)

    prefix = args.prefix or os.path.splitext(os.path.basename(args.input))[0]
    fmt    = args.format

    info("SEXmer scan plotter starting")
    info(f"Input  : {args.input}")
    info(f"Format : {fmt.upper()}")

    df = load_data(args.input)

    create_sexplot(df,        f"{prefix}_sexplot.{fmt}")
    create_abundance_plot(df, f"{prefix}_abundance.{fmt}")

    info("SEXmer scan plotter complete.")

if __name__ == "__main__":
    main()

PYPLOT
    chmod +x "$PLOTTER"
fi

# log run parameters
info "SEXmer scan starting"
info "Parameters: min-count=${MIN_COUNT}, max-count=${MAX_COUNT}, fold-threshold=${FOLD_THRESHOLD}, neutral-max=${NEUTRAL_MAX}, seed=${SEED}"
info "Settings  : mem=${MEM}, threads=${THREADS}"
info "Plotting  : enabled=${PLOT_ENABLED}, format=${PLOT_FORMAT}"
info "Temp dir  : ${SCAN_TMPDIR}"
info "Male files  (${#MALE_FILES[@]}): ${MALE_FILES[*]}"
info "Female files (${#FEMALE_FILES[@]}): ${FEMALE_FILES[*]}"

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
info "Identifying MSK candidates (male-consistent, pooled count ${MIN_COUNT}-${MAX_COUNT})..."

join "$SCAN_TMPDIR/male_consistent.txt" "$SCAN_TMPDIR/male_agg.txt" \
| awk -v mn="$MIN_COUNT" -v mx="$MAX_COUNT" '$2 >= mn && $2 <= mx' \
> "$SCAN_TMPDIR/msk_candidates.txt"

info "  $(wc -l < "$SCAN_TMPDIR/msk_candidates.txt") MSK candidates before cross-sex filter."

# STEP 4: Identify FSK candidates
info "Identifying FSK candidates (female-consistent, pooled count ${MIN_COUNT}-${MAX_COUNT})..."

join "$SCAN_TMPDIR/female_consistent.txt" "$SCAN_TMPDIR/female_agg.txt" \
| awk -v mn="$MIN_COUNT" -v mx="$MAX_COUNT" '$2 >= mn && $2 <= mx' \
> "$SCAN_TMPDIR/fsk_candidates.txt"

info "  $(wc -l < "$SCAN_TMPDIR/fsk_candidates.txt") FSK candidates before cross-sex filter."

# STEP 5: Cross-sex exclusion
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
info "Annotating k-mers with categories..."

awk -v ft="$FOLD_THRESHOLD" \
    -v nm="$N_MALE" \
    -v nf="$N_FEMALE" \
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
        else if (mc > 0 && fc > 0 && mc * nf < ft * fc * nm && fc * nm < ft * mc * nf) { cat = "neutral" }
        else                                                         { next            }

        print kmer "\t" mc "\t" fc "\t" cat
    }
' "$SCAN_TMPDIR/msk_kmers.txt" "$SCAN_TMPDIR/fsk_kmers.txt" \
  "$SCAN_TMPDIR/mbk_kmers.txt" "$SCAN_TMPDIR/fbk_kmers.txt" \
  "$SCAN_TMPDIR/union_counts.txt" \
> "$SCAN_TMPDIR/annotated_all.txt"

info "  $(wc -l < "$SCAN_TMPDIR/annotated_all.txt") total annotated k-mers."

# STEP 9: Split categories and sub-sample neutral k-mers
info "Splitting categories and sub-sampling neutral k-mers if needed..."

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
    awk -v max="$NEUTRAL_MAX" -v total="$NEUTRAL_TOTAL" -v seed="$SEED" '
    BEGIN { _a=1664525; _c=1013904223; _m=4294967296; _rng=seed+0 }
    function _rand() { _rng=(_a*_rng+_c)%_m; return _rng/_m }
    {
        remaining = total - NR + 1
        if (max > 0 && (remaining <= max || _rand() < max / remaining)) {
            print; max--
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
TSV_OUT="${OUTPUT_PREFIX}.kmers.tsv"
info "Writing main TSV output..."

{
    printf "kmer\tmale_count\tfemale_count\tcategory\n"
    cat "$SCAN_TMPDIR/non_neutral.txt" "$NEUTRAL_FILE"
} > "$TSV_OUT"

output "Main k-mer table written to: ${TSV_OUT}"

# STEP 11: Write FASTA outputs
MSK_FA="${OUTPUT_PREFIX}.MSK.fa"
FSK_FA="${OUTPUT_PREFIX}.FSK.fa"

info "Writing MSK FASTA output..."
awk -F'\t' 'BEGIN { n=0 }
     NR > 1 && $4 == "MSK" { n++; printf ">MSK_%d count=%s\n%s\n", n, $2, $1 }
' "$TSV_OUT" > "$MSK_FA"

info "Writing FSK FASTA output..."
awk -F'\t' 'BEGIN { n=0 }
     NR > 1 && $4 == "FSK" { n++; printf ">FSK_%d count=%s\n%s\n", n, $3, $1 }
' "$TSV_OUT" > "$FSK_FA"

output "MSK sequences written to: ${MSK_FA}"
output "FSK sequences written to: ${FSK_FA}"

# STEP 12: Generate visualization plots
if [[ "$PLOT_ENABLED" -eq 1 ]]; then
    info "Generating visualization plots from main TSV output..."
    python3 "$PLOTTER" -i "$TSV_OUT" --prefix "$OUTPUT_PREFIX" --format "$PLOT_FORMAT"
    SEX_PLOT="${OUTPUT_PREFIX}_sexplot.${PLOT_FORMAT}"
    ABUNDANCE_PLOT="${OUTPUT_PREFIX}_abundance.${PLOT_FORMAT}"
    [[ -f "$SEX_PLOT" ]] || { error "Plotter did not produce output file: $SEX_PLOT"; exit 1; }
    [[ -f "$ABUNDANCE_PLOT" ]] || { error "Plotter did not produce output file: $ABUNDANCE_PLOT"; exit 1; }
    output "Sex plot written to: ${SEX_PLOT}"
    output "Abundance plot written to: ${ABUNDANCE_PLOT}"
else
    info "Plot generation skipped (--no-plot)."
fi

# STEP 13: Tally category counts
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