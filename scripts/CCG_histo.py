#!/usr/bin/env python3
"""
Plot CCG count histograms for each CCG7/CCG10 FASTA file.
Counts CCG/CGG directly from sequences using CAACAGCCGCCA (fwd) and AGGAGGCGG (rev) anchors.
"""

import os
import glob
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FASTA_DIR = "/scratch/SCWF00044/Ruban/Lymphoblastoids_proline_phasing/260510/PB045/ccg_7_10_fasta/"
PLOT_DIR = os.path.join(FASTA_DIR, "plots")
MAX_CCG = 15

os.makedirs(PLOT_DIR, exist_ok=True)


def count_ccg_in_read(seq):
    seq = seq.upper()
    best = 0

    # forward
    start = 0
    while True:
        pos = seq.find("CAACAGCCGCCA", start)
        if pos == -1:
            break
        ccg_start = pos + 12
        cctcct = seq.find("CCTCCT", ccg_start)
        if cctcct != -1:
            region = seq[ccg_start:cctcct]
        else:
            tcct = seq.find("TCCT", ccg_start)
            if tcct != -1:
                region = seq[ccg_start:tcct]
            else:
                start = pos + 1
                continue
        c = 0
        for i in range(0, len(region) - 2, 3):
            if region[i:i+3] == "CCG":
                c += 1
            else:
                break
        if c > best:
            best = c
        start = pos + 1

    # reverse
    start = 0
    while True:
        pos = seq.find("AGGAGGCGG", start)
        if pos == -1:
            break
        cgg_start = pos + 6
        tggcgg = seq.find("TGGCGG", cgg_start)
        if tggcgg != -1:
            region = seq[cgg_start:tggcgg]
        else:
            end = -1
            for b in ["CTGTTG", "GAAG"]:
                p = seq.find(b, cgg_start)
                if p != -1 and (end == -1 or p < end):
                    end = p
            if end != -1:
                region = seq[cgg_start:end]
            else:
                start = pos + 1
                continue
        c = 0
        for i in range(0, len(region) - 2, 3):
            if region[i:i+3] == "CGG":
                c += 1
            else:
                break
        if c > best:
            best = c
        start = pos + 1

    return best


def parse_fasta(path):
    rid = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if rid:
                    yield rid, ''.join(seq)
                rid = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if rid:
        yield rid, ''.join(seq)


def plot_histogram(counts, title, out_path):
    if not counts:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    bins = [x - 0.5 for x in range(1, MAX_CCG + 2)]
    vals, _, bars = ax.hist(counts, bins=bins, color="#2196F3",
                            edgecolor="white", alpha=0.85, rwidth=0.8)

    for bar, v in zip(bars, vals):
        if v > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f"{int(v)}", ha="center", va="bottom", fontsize=9)

    ax.set_xlabel("CCG count", fontsize=12)
    ax.set_ylabel("Number of reads", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xticks(range(1, MAX_CCG + 1))
    ax.set_xlim(0.5, MAX_CCG + 0.5)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    fastas = sorted(glob.glob(os.path.join(FASTA_DIR, "*.fasta")))
    if not fastas:
        print(f"No FASTA files found in {FASTA_DIR}")
        return

    for fasta_path in fastas:
        fname = os.path.basename(fasta_path)
        m = re.match(r'(bc\d+)_(CCG\d+)\.fasta', fname)
        if not m:
            continue

        bc = m.group(1)
        ccg_label = m.group(2)

        counts = []
        for rid, seq in parse_fasta(fasta_path):
            c = count_ccg_in_read(seq)
            if 0 < c <= MAX_CCG:
                counts.append(c)

        title = f"CCG distribution — {bc} {ccg_label} ({len(counts)} reads)"
        out_path = os.path.join(PLOT_DIR, f"{bc}_{ccg_label}_ccg_histogram.png")
        plot_histogram(counts, title, out_path)
        print(f"{fname} -> {out_path}  ({len(counts)} reads)")

    print("Done.")


if __name__ == "__main__":
    main()
