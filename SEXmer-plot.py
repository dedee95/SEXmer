#!/usr/bin/env python3
# SEXmer-plot.py - Generate publication-ready figures from SEXmer scan output.
# Author: Dede Kurniawan

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
  --format             Output format: png or pdf         [default: png]
  -h, --help           Show this help and exit

Output files:
  {prefix}_sexplot.{format}    scatter + bar plot
  {prefix}_abundance.{format}  raw scatter + smooth density plot
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
    parser.add_argument('--format',       choices=['png', 'pdf'], default='png')
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

    info("SEXmer-plot starting")
    info(f"Input  : {args.input}")
    info(f"Format : {fmt.upper()}")

    df = load_data(args.input)

    create_sexplot(df,        f"{prefix}_sexplot.{fmt}")
    create_abundance_plot(df, f"{prefix}_abundance.{fmt}")

    info("SEXmer-plot complete.")

if __name__ == "__main__":
    main()
