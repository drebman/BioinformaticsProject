#!/usr/bin/env python3
"""
plot_gene_lengths.py

Generates:
1) Overlaid histograms of gene length distributions for each comparison pair
2) A scatter plot of genome size (bp) vs. median gene length (bp) across all genomes

Usage (from your repo root):
    python plot_gene_lengths.py \
        --results_dir results \
        --stats genome_stats.tsv \
        --outdir figures \
        --bins 60

Requirements:
    - pandas
    - matplotlib

Assumptions:
    - You have already run compare_gene_lengths.py to produce:
        results/buchnera_vs_ecoli.Buchnera.lengths.tsv
        results/buchnera_vs_ecoli.Ecoli_K12.lengths.tsv
        results/wigglesworthia_vs_salmonella.Wigglesworthia.lengths.tsv
        results/wigglesworthia_vs_salmonella.Salmonella_no75.lengths.tsv
        results/cpneumoniae_vs_isosphaera.C_pneumoniae_TW183.lengths.tsv
        results/cpneumoniae_vs_isosphaera.Isosphaera_pallida.lengths.tsv
    - You have genome_stats.tsv created by compute_genome_stats.py with columns:
        label, path, genome_size_bp, gc_content_pct
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def load_lengths(tsv_path):
    """Read a .lengths.tsv file and return a pandas Series of lengths (bp)."""
    df = pd.read_csv(tsv_path, sep="\t")
    # Expect columns: index, length_bp
    col = "length_bp"
    if col not in df.columns:
        # Support potential alternative column naming
        candidates = [c for c in df.columns if "length" in c.lower()]
        if not candidates:
            raise ValueError(f"Could not find 'length_bp' column in {tsv_path}")
        col = candidates[0]
    return df[col].dropna()

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def plot_pair_hist(lengths1, lengths2, label1, label2, out_path, bins=60):
    """Create an overlaid histogram for two genomes' gene lengths."""
    plt.figure(figsize=(8,6))
    plt.hist(lengths1, bins=bins, label=label1, alpha=0.5)
    plt.hist(lengths2, bins=bins, label=label2, alpha=0.5)
    plt.xlabel("Gene Length (bp)")
    plt.ylabel("Count")
    plt.title(f"Gene Length Distributions: {label1} vs {label2}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

def plot_genome_size_vs_median(stats_df, label_to_median, out_path):
    """Scatter of genome size vs median gene length for all available labels."""
    # Merge stats_df['label', 'genome_size_bp'] with medians
    med_df = pd.DataFrame(
        [(lab, med) for lab, med in label_to_median.items()],
        columns=["label", "median_gene_length"]
    )
    merged = stats_df.merge(med_df, on="label", how="inner")

    if merged.empty:
        raise ValueError("No overlap between genome_stats labels and computed medians.")

    plt.figure(figsize=(8,6))
    plt.scatter(merged["genome_size_bp"], merged["median_gene_length"], s=70)
    # annotate points
    for _, row in merged.iterrows():
        plt.text(row["genome_size_bp"], row["median_gene_length"], row["label"], fontsize=8, ha="left", va="bottom")
    plt.xlabel("Genome Size (bp)")
    plt.ylabel("Median Gene Length (bp)")
    plt.title("Genome Size vs. Median Gene Length")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

def main():
    ap = argparse.ArgumentParser(description="Plot gene length distributions and genome size trends.")
    ap.add_argument("--results_dir", default="results", help="Directory containing *.lengths.tsv files (default: results)")
    ap.add_argument("--stats", default="genome_stats.tsv", help="Path to genome_stats.tsv (default: genome_stats.tsv)")
    ap.add_argument("--outdir", default="figures", help="Directory to write figures (default: figures)")
    ap.add_argument("--bins", type=int, default=60, help="Number of histogram bins (default: 60)")
    args = ap.parse_args()

    ensure_dir(args.outdir)

    # Define your three comparison pairs and their labels exactly as produced by compare_gene_lengths.py
    pairs = [
        # (prefix, label1, label2)
        ("buchnera_vs_ecoli", "Buchnera", "Ecoli_K12"),
        ("wigglesworthia_vs_salmonella", "Wigglesworthia", "Salmonella_no75"),
        ("cpneumoniae_vs_isosphaera", "C_pneumoniae_TW183", "Isosphaera_pallida"),
    ]

    # Plot histograms per pair and collect medians
    label_to_median = {}
    for prefix, lab1, lab2 in pairs:
        f1 = os.path.join(args.results_dir, f"{prefix}.{lab1}.lengths.tsv")
        f2 = os.path.join(args.results_dir, f"{prefix}.{lab2}.lengths.tsv")
        if not os.path.exists(f1) or not os.path.exists(f2):
            print(f"[WARN] Missing lengths file(s) for pair {prefix}. Skipping histogram.", flush=True)
        else:
            len1 = load_lengths(f1)
            len2 = load_lengths(f2)
            out_png = os.path.join(args.outdir, f"{prefix}.hist.png")
            plot_pair_hist(len1, len2, lab1, lab2, out_png, bins=args.bins)
            print(f"[OK] Wrote {out_png}")

            # Save medians
            label_to_median[lab1] = float(len1.median())
            label_to_median[lab2] = float(len2.median())

    # Load genome stats and make scatter
    if os.path.exists(args.stats):
        stats_df = pd.read_csv(args.stats, sep="\t")
        # Expect columns: label, path, genome_size_bp, gc_content_pct
        if not {"label", "genome_size_bp"}.issubset(stats_df.columns):
            raise ValueError("genome_stats.tsv must contain columns: label, genome_size_bp")
        scatter_path = os.path.join(args.outdir, "genome_size_vs_median_gene_length.png")
        plot_genome_size_vs_median(stats_df, label_to_median, scatter_path)
        print(f"[OK] Wrote {scatter_path}")
    else:
        print(f"[WARN] Stats file not found: {args.stats}. Skipping genome-size scatter.", flush=True)

if __name__ == "__main__":
    main()
