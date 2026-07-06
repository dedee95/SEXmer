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
DEFAULT_FORMAT = "svg"

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

def choose_figsize(chr_count: int, figsize: Optional[Tuple[float, float]]) -> Tuple[float, float]:
    if figsize is not None:
        return figsize
    auto_width = min(max(14.0, chr_count * 0.55), 28.0)
    return (auto_width, 5.3)

def select_region(
    df: pd.DataFrame,
    chr_info: Sequence[dict],
    target_chr: Optional[str],
    start: Optional[int],
    end: Optional[int],
) -> Tuple[pd.DataFrame, List[dict]]:
    if target_chr is None:
        if start is not None or end is not None:
            raise ValueError("--start/--end require --chr")
        return df, list(chr_info)

    selected = df[df["chr"] == target_chr].copy()
    if selected.empty:
        available = ", ".join(str(item["chr"]) for item in chr_info[:20])
        suffix = " ..." if len(chr_info) > 20 else ""
        raise ValueError(f"Chromosome not found: {target_chr}. Available: {available}{suffix}")

    zoom_start = start if start is not None else float(selected["start"].min())
    zoom_end = end if end is not None else float(selected["end"].max())
    if zoom_end <= zoom_start:
        raise ValueError("--end must be larger than --start")

    selected = selected[(selected["pos"] >= zoom_start) & (selected["pos"] <= zoom_end)].copy()
    if selected.empty:
        raise ValueError(f"No data in region {int(zoom_start):,}-{int(zoom_end):,} on {target_chr}")

    base_info = next(item for item in chr_info if item["chr"] == target_chr)
    selected["cum_pos"] = selected["pos"]
    selected_info = [
        {
            "chr": target_chr,
            "start": float(zoom_start),
            "end": float(zoom_end),
            "center": (float(zoom_start) + float(zoom_end)) / 2.0,
            "color": base_info["color"],
        }
    ]
    return selected, selected_info

def create_plot(
    df: pd.DataFrame,
    chr_info: Sequence[dict],
    output_file: Path,
    y_label: str,
    plot_type: str,
    figsize: Optional[Tuple[float, float]],
    target_chr: Optional[str],
    start: Optional[int],
    end: Optional[int],
    ymax: Optional[float],
) -> None:
    plot_df, plot_chr_info = select_region(df, chr_info, target_chr, start, end)
    fig, ax = plt.subplots(figsize=choose_figsize(len(plot_chr_info), figsize), dpi=DPI)

    for item in plot_chr_info:
        chrom = item["chr"]
        chr_df = plot_df[plot_df["chr"] == chrom].sort_values("cum_pos")
        if chr_df.empty:
            continue
        x = chr_df["cum_pos"] / 1e6
        y = chr_df["plot_value"]
        if plot_type == "line":
            ax.fill_between(x, 0, y, color=item["color"], alpha=0.25, linewidth=0, rasterized=True)
            ax.plot(x, y, color=item["color"], linewidth=0.9, alpha=0.9, rasterized=True)
        elif plot_type == "bar":
            step = chr_df["cum_pos"].diff().median()
            width_mb = (step if pd.notna(step) and step > 0 else 100_000) / 1e6
            ax.bar(x, y, color=item["color"], width=width_mb, alpha=0.9, edgecolor="none", rasterized=True)
        else:
            ax.scatter(x, y, c=item["color"], s=10, alpha=0.75, edgecolors="none", rasterized=True)

    if target_chr is None:
        tick_positions = [0.0] + [item["end"] / 1e6 for item in plot_chr_info]
        tick_labels = [""] + [strip_chr_prefix(item["chr"]) for item in plot_chr_info]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontweight="bold", rotation=0)
        x_left = 0.0
        x_right = max(tick_positions) if tick_positions else float(plot_df["cum_pos"].max() / 1e6)
        ax.set_xlabel("Chromosome", fontsize=18, labelpad=10)
    else:
        x_left = (start if start is not None else float(plot_df["pos"].min())) / 1e6
        x_right = (end if end is not None else float(plot_df["pos"].max())) / 1e6
        ax.set_xlabel(f"{target_chr} position (Mb)", fontsize=18, labelpad=10)

    x_range = x_right - x_left
    x_pad = x_range * 0.01 if x_range > 0 else 0.5
    ax.set_xlim(left=x_left - x_pad, right=x_right + x_pad)

    ax.set_ylabel(y_label, fontsize=18, labelpad=10)

    y_upper = float(ymax) if ymax is not None else nice_y_upper(plot_df["plot_value"])
    if not math.isfinite(y_upper) or y_upper <= 0:
        raise ValueError("--ymax must be a positive number")

    y_lower = -0.04 * y_upper
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
    ax.spines["bottom"].set_position(("outward", 0))
    ax.spines["left"].set_bounds(0, y_upper)
    ax.spines["bottom"].set_bounds(x_left, x_right)

    fig.tight_layout()
    fig.savefig(output_file, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def default_output_stem(
    input_path: Path,
    input_type: str,
    target_chr: Optional[str],
    start: Optional[int],
    end: Optional[int],
) -> str:
    stem = input_path.name
    for suffix in (".kmer.windows.tsv", ".reads.windows.tsv", ".windows.tsv", ".tsv"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break

    parts = [stem, input_type, "manhattan"]
    if target_chr is not None:
        parts.append(str(target_chr))
    if start is not None or end is not None:
        s = start if start is not None else 0
        e = end if end is not None else 0
        parts.append(f"{s}-{e}")
    return ".".join(parts)

class SexmerArgumentParser(argparse.ArgumentParser):
    def format_help(self) -> str:
        return self.description or ""

def parse_args() -> argparse.Namespace:
    help_text = f"""
SEXmer manhattan v3 - Visualize SEXmer map window TSV outputs.

Usage: manhattan_v3.py <window.tsv> [window.tsv ...] [OPTIONS]

Mandatory:
  <window.tsv>         Specify SEXmer map window TSV file(s)
                       Accepts kmer.windows.tsv and reads.windows.tsv outputs

Optional:
  -c, --chr            Plot one chromosome only, e.g. chr5
  --start              Start position for zoom region in bp             [requires: --chr]
  --end                End position for zoom region in bp               [requires: --chr]
  --ymax               Specify maximum value for the Y-axis
  --figsize            Specify figure size in inches, e.g. --figsize 16 6
  -p, --plot-type      Specify plot type: scatter, line, or bar          [default: scatter]
                       Note: line plot is drawn with filled area under the line
  -f, --format         Specify output format: svg, png, or pdf           [default: {DEFAULT_FORMAT}]
  -o, --output         Specify output file name                         [single input only]
  -h, --help           Show this help message and exit
"""
    parser = SexmerArgumentParser(
        add_help=False,
        usage=argparse.SUPPRESS,
        description=help_text,
    )
    parser.add_argument("inputs", nargs="+", metavar="<window.tsv>")
    parser.add_argument("--chr", "-c")
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--ymax", type=float)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("WIDTH", "HEIGHT"))
    parser.add_argument("--plot-type", "-p", choices=["scatter", "line", "bar"], default="scatter")
    parser.add_argument("--format", "-f", choices=["svg", "png", "pdf"], default=DEFAULT_FORMAT)
    parser.add_argument("--output", "-o")
    parser.add_argument("--help", "-h", action="help")
    return parser.parse_args()

def main() -> int:
    args = parse_args()
    logger = setup_logger()

    try:
        if args.output and len(args.inputs) != 1:
            raise ValueError("--output can only be used with one input file")

        figsize = tuple(args.figsize) if args.figsize is not None else None

        for input_name in args.inputs:
            input_path = Path(input_name)
            if not input_path.exists():
                raise FileNotFoundError(f"Input file not found: {input_path}")

            logger.info(f"Reading input data: {input_path}")
            raw_df, input_type, y_label = load_sexmer_table(input_path)
            logger.info(f"Detected input type: {input_type}; Y-axis: {y_label}")

            logger.info("Preparing Manhattan plot data")
            plot_df, chr_info = prepare_manhattan_data(raw_df)

            if args.output:
                output_file = Path(args.output)
            else:
                output_file = input_path.parent / f"{default_output_stem(input_path, input_type, args.chr, args.start, args.end)}.{args.format}"

            logger.info(f"Generating {args.format.upper()} Manhattan plot")
            create_plot(
                df=plot_df,
                chr_info=chr_info,
                output_file=output_file,
                y_label=y_label,
                plot_type=args.plot_type,
                figsize=figsize,
                target_chr=args.chr,
                start=args.start,
                end=args.end,
                ymax=args.ymax,
            )
            log_output(logger, f"Plot file: {output_file}")

        logger.info("Done")
        return 0
    except Exception as exc:
        logger.error(str(exc))
        return 1

if __name__ == "__main__":
    sys.exit(main())