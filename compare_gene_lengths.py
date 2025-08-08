#!/usr/bin/env python3
import argparse
import math
import statistics as stats
from collections import defaultdict

def parse_gff_attributes(attr_field):
    out = {}
    if not attr_field or attr_field == ".":
        return out
    parts = attr_field.strip().split(";")
    for p in parts:
        if not p:
            continue
        if "=" in p:
            k, v = p.split("=", 1)
        elif " " in p:
            k, v = p.split(" ", 1)
        else:
            k, v = p, ""
        out[k.strip()] = v.strip()
    return out

def extract_lengths_from_gff(path, feature_type="CDS", id_keys=("ID","locus_tag","gene","Name","protein_id")):
    lengths = []
    id_to_len = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = cols
            if ftype != feature_type:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            length = abs(end_i - start_i) + 1
            lengths.append(length)
            attrd = parse_gff_attributes(attrs)
            chosen_id = None
            for k in id_keys:
                if k in attrd and attrd[k] and attrd[k] != ".":
                    chosen_id = attrd[k]
                    break
            if chosen_id is not None and chosen_id not in id_to_len:
                id_to_len[chosen_id] = length
    return lengths, id_to_len

def mann_whitney_u(x, y):
    combined = [(val, 0) for val in x] + [(val, 1) for val in y]
    combined.sort(key=lambda t: t[0])
    ranks = {}
    i = 0
    N = len(combined)
    while i < N:
        j = i + 1
        while j < N and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j
    R1 = sum(ranks[i] for i in range(N) if combined[i][1] == 0)
    R2 = sum(ranks[i] for i in range(N) if combined[i][1] == 1)
    n1 = len(x); n2 = len(y)
    U1 = R1 - n1*(n1+1)/2.0
    U2 = R2 - n2*(n2+1)/2.0
    U = min(U1, U2)
    mu_U = n1*n2/2.0
    tie_counts = defaultdict(int)
    for val, _ in combined:
        tie_counts[val] += 1
    T = sum(c**3 - c for c in tie_counts.values())
    sigma_U = math.sqrt(n1*n2*(n1+n2+1 - T/((n1+n2)*(n1+n2-1))) / 12.0)
    if sigma_U == 0:
        return U, float('nan'), 1.0
    z = (U - mu_U + 0.5) / sigma_U if U < mu_U else (U - mu_U - 0.5) / sigma_U
    p = 2.0 * 0.5 * (1 - math.erf(abs(z)/math.sqrt(2)))
    return U, z, p

def wilcoxon_signed_rank(x, y):
    diffs = [a - b for a, b in zip(x, y)]
    pairs = [(abs(d), 1 if d > 0 else -1) for d in diffs if d != 0]
    if not pairs:
        return 0.0, float('nan'), 1.0
    pairs.sort(key=lambda t: t[0])
    N = len(pairs)
    ranks = [0.0]*N
    i = 0
    while i < N:
        j = i + 1
        while j < N and pairs[j][0] == pairs[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j
    Wpos = sum(r for r, (_, sgn) in zip(ranks, pairs) if sgn > 0)
    Wneg = sum(r for r, (_, sgn) in zip(ranks, pairs) if sgn < 0)
    W = min(Wpos, Wneg)
    mu = N*(N+1)/4.0
    tie_counts = defaultdict(int)
    for a, _ in pairs:
        tie_counts[a] += 1
    T = sum(c**3 - c for c in tie_counts.values())
    sigma = math.sqrt(N*(N+1)*(2*N+1)/24.0 - T/48.0)
    if sigma == 0:
        return W, float('nan'), 1.0
    z = (W - mu - 0.5) / sigma if W < mu else (W - mu + 0.5) / sigma
    p = 2.0 * 0.5 * (1 - math.erf(abs(z)/math.sqrt(2)))
    return W, z, p

def summarize(vec):
    if not vec:
        return {"n":0, "mean":float('nan'), "median":float('nan'), "iqr":(float('nan'), float('nan'))}
    q1 = stats.quantiles(vec, n=4, method="inclusive")[0]
    q3 = stats.quantiles(vec, n=4, method="inclusive")[2]
    return {"n": len(vec), "mean": stats.fmean(vec), "median": stats.median(vec), "iqr": (q1, q3)}

def write_tsv(path, rows, header):
    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")

def main():
    ap = argparse.ArgumentParser(description="Compare gene/CDS length distributions between two genomes using GFF3.")
    ap.add_argument("--gff1", required=True, help="Path to first genome GFF3 (e.g., symbiont)")
    ap.add_argument("--gff2", required=True, help="Path to second genome GFF3 (e.g., free-living)")
    ap.add_argument("--label1", default="genome1", help="Label for genome 1")
    ap.add_argument("--label2", default="genome2", help="Label for genome 2")
    ap.add_argument("--feature", choices=["CDS","gene"], default="CDS", help="Feature type to measure (default: CDS)")
    ap.add_argument("--pair_key", default="locus_tag", help="GFF attribute key to attempt pairing on (e.g., locus_tag, gene, Name, protein_id)")
    ap.add_argument("--out_prefix", default="comparison", help="Prefix for output files")
    args = ap.parse_args()

    lengths1, idmap1 = extract_lengths_from_gff(args.gff1, feature_type=args.feature)
    lengths2, idmap2 = extract_lengths_from_gff(args.gff2, feature_type=args.feature)

    write_tsv(f"{args.out_prefix}.{args.label1}.lengths.tsv",
              [(i+1, L) for i, L in enumerate(lengths1)], ["index","length_bp"])
    write_tsv(f"{args.out_prefix}.{args.label2}.lengths.tsv",
              [(i+1, L) for i, L in enumerate(lengths2)], ["index","length_bp"])

    U, z_u, p_u = mann_whitney_u(lengths1, lengths2)

    paired1 = []
    paired2 = []
    intersect_ids = set(idmap1.keys()).intersection(set(idmap2.keys()))
    if intersect_ids:
        for k in intersect_ids:
            paired1.append(idmap1[k])
            paired2.append(idmap2[k])
        W, z_w, p_w = wilcoxon_signed_rank(paired1, paired2)
    else:
        W, z_w, p_w = float('nan'), float('nan'), float('nan')

    s1 = summarize(lengths1)
    s2 = summarize(lengths2)

    with open(f"{args.out_prefix}.stats.txt", "w", encoding="utf-8") as out:
        out.write(f"Feature: {args.feature}\n")
        out.write(f"Genome 1 ({args.label1}): n={s1['n']}, mean={s1['mean']:.2f}, median={s1['median']:.2f}, IQR=({s1['iqr'][0]:.2f},{s1['iqr'][1]:.2f})\n")
        out.write(f"Genome 2 ({args.label2}): n={s2['n']}, mean={s2['mean']:.2f}, median={s2['median']:.2f}, IQR=({s2['iqr'][0]:.2f},{s2['iqr'][1]:.2f})\n\n")
        out.write("Unpaired comparison (Mann-Whitney U, two-sided normal approx):\n")
        out.write(f"U={U:.2f}, z={z_u:.3f}, p={p_u:.3g}\n\n")
        if not math.isnan(p_w):
            out.write(f"Paired comparison over {len(paired1)} shared IDs (Wilcoxon signed-rank, two-sided normal approx):\n")
            out.write(f"W={W:.2f}, z={z_w:.3f}, p={p_w:.3g}\n")
        else:
            out.write("Paired comparison: no shared IDs found across chosen attributes; skipped.\n")

    print(f"Wrote: {args.out_prefix}.{args.label1}.lengths.tsv")
    print(f"Wrote: {args.out_prefix}.{args.label2}.lengths.tsv")
    print(f"Wrote: {args.out_prefix}.stats.txt")
    if not math.isnan(p_w):
        print(f"Paired on {len(paired1)} shared IDs. See stats file for details.")
    else:
        print("No shared IDs for pairing; only unpaired test reported.")

if __name__ == "__main__":
    main()
