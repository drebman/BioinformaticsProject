import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------- FASTA I/O ----------
def read_fasta(path: str) -> str:
    seq = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith(">"):
                continue
            seq.append(s)
    return "".join(seq).upper()


# ---------- Needleman–Wunsch (simple global alignment) ----------
def needleman_wunsch(s1: str, s2: str, match=2, mismatch=-1, gap=-2) -> tuple[str, str]:
    n, m = len(s1), len(s2)
    score = np.zeros((n + 1, m + 1), dtype=int)
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)  # 0 diag, 1 up, 2 left

    for i in range(1, n + 1):
        score[i, 0] = score[i - 1, 0] + gap
        trace[i, 0] = 1
    for j in range(1, m + 1):
        score[0, j] = score[0, j - 1] + gap
        trace[0, j] = 2

    for i in range(1, n + 1):
        a = s1[i - 1]
        for j in range(1, m + 1):
            b = s2[j - 1]
            diag = score[i - 1, j - 1] + (match if a == b else mismatch)
            up = score[i - 1, j] + gap
            left = score[i, j - 1] + gap
            best = diag if diag >= up and diag >= left else (up if up >= left else left)
            score[i, j] = best
            trace[i, j] = 0 if best == diag else (1 if best == up else 2)

    i, j = n, m
    a1, a2 = [], []
    while i > 0 or j > 0:
        t = trace[i, j]
        if t == 0:
            a1.append(s1[i - 1]); a2.append(s2[j - 1]); i -= 1; j -= 1
        elif t == 1:
            a1.append(s1[i - 1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(s2[j - 1]); j -= 1
    return "".join(reversed(a1)), "".join(reversed(a2))


# ---------- Find deletions in symbiont along the reference (free-living) coordinate ----------
def find_losses(aln_sym: str, aln_free: str) -> list[tuple[int, int]]:
    free_pos = -1
    losses, in_del, start = [], False, None
    for a, b in zip(aln_sym, aln_free):
        if b != "-":                       # advance along reference coordinate
            free_pos += 1
            if a == "-":                   # deletion in symbiont
                if not in_del:
                    in_del = True; start = free_pos
            else:
                if in_del:                 # close the deletion run
                    losses.append((start, free_pos - 1))
                    in_del = False
    if in_del:
        losses.append((start, free_pos))
    return losses


# ---------- Plot & CSV (UPDATED: adds 'gene' and 'comparison' columns) ----------
def write_outputs(losses, aln_sym, aln_free, free_label, sym_label, gene, outbase):
    free_len = sum(1 for c in aln_free if c != "-")
    sym_len  = sum(1 for c in aln_sym  if c != "-")
    matches  = sum(1 for a, b in zip(aln_sym, aln_free) if a == b and a != "-" and b != "-")
    compared = sum(1 for a, b in zip(aln_sym, aln_free) if a != "-" and b != "-")
    identity = matches / compared if compared else 0.0
    total_loss = sum(e - s + 1 for s, e in losses)
    loss_frac  = total_loss / free_len if free_len else 0.0

    # CSV of loss segments (now includes 'gene' and 'comparison')
    rows = [{
        "gene": gene,
        "comparison": f"{sym_label} vs {free_label}",
        f"start_{free_label}": s,
        f"end_{free_label}": e,
        "length": e - s + 1
    } for (s, e) in losses]
    df = pd.DataFrame(rows, columns=[
        "gene", "comparison", f"start_{free_label}", f"end_{free_label}", "length"
    ])
    df.to_csv(outbase + "_loss_segments.csv", index=False)

    # Box diagram
    fig, ax = plt.subplots(figsize=(10, 1.8))
    ax.set_title(f"{gene} — Regions missing in {sym_label} relative to {free_label}")
    ax.add_patch(plt.Rectangle((0, 0.4), max(1, free_len), 0.2, fill=False, linewidth=2))
    for s, e in losses:
        ax.add_patch(plt.Rectangle((s, 0.3), e - s + 1, 0.4, alpha=0.85))
    ax.set_xlim(0, max(1, free_len))
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel(f"{free_label} coordinate (aa)")
    ax.text(max(1, free_len)*0.01, 0.9, f"{free_label} length: {free_len} aa", fontsize=9, va="center")
    ax.text(max(1, free_len)*0.35, 0.9, f"{sym_label} length: {sym_len} aa", fontsize=9, va="center")
    ax.text(max(1, free_len)*0.65, 0.9, f"Loss: {total_loss} aa ({loss_frac:.1%})", fontsize=9, va="center")
    fig.tight_layout()
    fig.savefig(outbase + "_loss_box_diagram.png", dpi=300, bbox_inches="tight")
    fig.savefig(outbase + "_loss_box_diagram.pdf", dpi=300, bbox_inches="tight")

    print(f"[saved] {outbase}_loss_box_diagram.png/.pdf and {outbase}_loss_segments.csv")
    print(f"       identity≈{identity:.4f}, loss={total_loss} aa ({loss_frac:.2%})")


def compare_pair(sym_fa: str, free_fa: str, gene: str,
                 sym_label: str, free_label: str, outbase: str):
    sym = read_fasta(sym_fa)
    ref = read_fasta(free_fa)
    aln_sym, aln_free = needleman_wunsch(sym, ref)        # global alignment
    losses = find_losses(aln_sym, aln_free)               # deletions in sym vs ref
    write_outputs(losses, aln_sym, aln_free, free_label, sym_label, gene, outbase)


# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Highlight losses in symbiont vs reference (free-living).")
    ap.add_argument("--sym", required=True, help="Symbiont FASTA (same gene)")
    ap.add_argument("--free", required=True, help="Reference (free-living) FASTA")
    ap.add_argument("--gene", required=True, help="Gene name (for titles and CSV)")
    ap.add_argument("--sym-label", default="Symbiont", help="Label for symbiont")
    ap.add_argument("--free-label", default="Free-living", help="Label for reference")
    ap.add_argument("--out", required=True, help="Output basename (no extension)")
    args = ap.parse_args()
    compare_pair(args.sym, args.free, args.gene, args.sym_label, args.free_label, args.out)


if __name__ == "__main__":
    main()
