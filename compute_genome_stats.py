#!/usr/bin/env python3
import argparse
import os
import sys

def parse_fasta(path):
    seqs = []
    cur = []
    try:
        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                if not line:
                    continue
                if line.startswith('>'):
                    if cur:
                        seqs.append(''.join(cur))
                        cur = []
                else:
                    cur.append(line.strip())
            if cur:
                seqs.append(''.join(cur))
    except FileNotFoundError:
        print(f"ERROR: File not found: {path}", file=sys.stderr)
        return []
    return seqs

def gc_content(seq):
    if not seq:
        return 0.0
    s = seq.upper()
    gc = s.count('G') + s.count('C')
    return gc / len(s)

def format_markdown_table(rows):
    # rows: list of dicts with keys: label, path, genome_size, gc_pct, genes
    header = "| Bacterium | FASTA Path | Genome Size (bp) | GC Content (%) | Genes compared |"
    sep    = "|---|---|---:|---:|---|"
    lines = [header, sep]
    for r in rows:
        lines.append(f"| {r['label']} | `{r['path']}` | {r['genome_size']} | {r['gc_pct']:.2f} | {r.get('genes','')} |")
    return "\n".join(lines)

def main():
    p = argparse.ArgumentParser(description="Compute genome size and GC% for FASTA files and print Markdown/TSV summaries (with 'Genes compared').")
    p.add_argument("files", nargs="+", help="Paths to FASTA files")
    p.add_argument("--labels", nargs="*", default=None, help="Optional labels for each file (same order as files). Default = basename.")
    p.add_argument("--genes", nargs="*", default=None, help="Optional genes-compared entry for each file (same order as files).")
    p.add_argument("--tsv", default="genome_stats.tsv", help="Path to write TSV summary (default: genome_stats.tsv)")
    p.add_argument("--md", default="genome_stats.md", help="Path to write Markdown table (default: genome_stats.md)")
    args = p.parse_args()

    if args.labels and len(args.labels) != len(args.files):
        print("ERROR: If provided, --labels must have same length as files.", file=sys.stderr)
        sys.exit(2)
    if args.genes and len(args.genes) != len(args.files):
        print("ERROR: If provided, --genes must have same length as files.", file=sys.stderr)
        sys.exit(2)

    rows = []
    for i, path in enumerate(args.files):
        label = args.labels[i] if args.labels else os.path.basename(path)
        genes = args.genes[i] if args.genes else ""
        seqs = parse_fasta(path)
        if not seqs:
            # Already printed error; continue to next file
            continue
        total_len = sum(len(s) for s in seqs)
        gc = gc_content(''.join(seqs))
        rows.append({
            "label": label,
            "path": path,
            "genome_size": total_len,
            "gc_pct": gc * 100.0,
            "genes": genes
        })

    if not rows:
        print("No valid FASTA inputs parsed; nothing to write.", file=sys.stderr)
        sys.exit(1)

    # Write TSV
    with open(args.tsv, "w", encoding="utf-8") as out:
        out.write("label\tpath\tgenome_size_bp\tgc_content_pct\tgenes_compared\n")
        for r in rows:
            out.write(f"{r['label']}\t{r['path']}\t{r['genome_size']}\t{r['gc_pct']:.4f}\t{r.get('genes','')}\n")

    # Write Markdown
    md = format_markdown_table(rows)
    with open(args.md, "w", encoding="utf-8") as out:
        out.write(md + "\n")

    # Also print Markdown table to stdout for convenience
    print(md)

if __name__ == "__main__":
    main()
