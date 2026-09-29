import argparse
import os

import matplotlib
matplotlib.use("Agg")  # headless-safe: this script only saves PNGs, never pops up a window
import matplotlib.pyplot as plt

from parser import parse_vcf
from statistics import compute_statistics

COLORS = {
    "quality": "#4C72B0",
    "chrom": "#55A868",
    "snp_indel": "#C44E52",
    "ts_tv": "#8172B2",
    "af": "#CCB974",
    "depth": "#64B5CD",
}


def _save_histogram(values, title, xlabel, path, color):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(values.dropna(), bins=20, color=color, edgecolor="black")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _save_bar(labels, values, title, ylabel, path, color):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(labels, values, color=color, edgecolor="black")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_quality_distribution(df, path):
    _save_histogram(df["qual"], "Variant Quality Distribution", "QUAL", path, COLORS["quality"])


def plot_variants_per_chromosome(df, path):
    counts = df["chrom"].value_counts().sort_index()
    _save_bar(counts.index.astype(str), counts.values, "Variants per Chromosome", "Count", path, COLORS["chrom"])


def plot_snp_indel_distribution(stats, path):
    counts = {k: v for k, v in stats["variant_counts"].items() if v > 0}
    _save_bar(list(counts.keys()), list(counts.values()), "SNP vs. Indel Distribution", "Count", path, COLORS["snp_indel"])


def plot_ts_tv_ratio(stats, path):
    ratio = stats["ts_tv_ratio"]
    ratio_label = f"{ratio:.2f}" if ratio == ratio else "N/A"
    _save_bar(
        ["Transitions", "Transversions"],
        [stats["ts_count"], stats["tv_count"]],
        f"Transition/Transversion Ratio (Ts/Tv = {ratio_label})",
        "Count",
        path,
        COLORS["ts_tv"],
    )


def plot_allele_frequency_distribution(df, path):
    _save_histogram(df["info_AF"], "Variant Allele Frequency Distribution", "Allele Frequency (AF)", path, COLORS["af"])


def plot_depth_distribution(df, path):
    _save_histogram(df["info_DP"], "Variant Depth Distribution", "Depth (DP)", path, COLORS["depth"])


def generate_all_plots(df, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    stats = compute_statistics(df)

    plots = {
        "quality_distribution.png": (lambda path: plot_quality_distribution(df, path), "qual"),
        "variants_per_chromosome.png": (lambda path: plot_variants_per_chromosome(df, path), "chrom"),
        "snp_indel_distribution.png": (lambda path: plot_snp_indel_distribution(stats, path), None),
        "ts_tv_ratio.png": (lambda path: plot_ts_tv_ratio(stats, path), None),
        "allele_frequency_distribution.png": (lambda path: plot_allele_frequency_distribution(df, path), "info_AF"),
        "depth_distribution.png": (lambda path: plot_depth_distribution(df, path), "info_DP"),
    }

    saved_paths = []
    for filename, (draw, required_column) in plots.items():
        if required_column is not None and required_column not in df.columns:
            print(f"Skipping {filename}: column '{required_column}' not present in this VCF")
            continue
        path = os.path.join(out_dir, filename)
        draw(path)
        saved_paths.append(path)

    return saved_paths


def main():
    arg_parser = argparse.ArgumentParser(description="Generate summary plots for a VCF file.")
    arg_parser.add_argument("vcf_path", help="Path to the input VCF file")
    arg_parser.add_argument("--out-dir", default="reports/figures", help="Directory to save the plots in (default: reports/figures)")
    args = arg_parser.parse_args()

    df = parse_vcf(args.vcf_path)
    saved_paths = generate_all_plots(df, args.out_dir)

    print(f"Saved {len(saved_paths)} plots to {args.out_dir}/:")
    for path in saved_paths:
        print(f"  {path}")


if __name__ == "__main__":
    main()
