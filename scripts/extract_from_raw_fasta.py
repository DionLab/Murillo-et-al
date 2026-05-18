#!/usr/bin/env python3
"""
Count CCG repeats from raw demultiplexed HiFi FASTA files and extract
reads with exactly 7 or 10 CCG repeats.

Counts only exact CCG/CGG codons between anchors:
  FWD: CAACAGCCGCCA .... CCTCCT  (count CCG after CCA)
  REV: AGGAGG .... TGGCGG        (count CGG, exclude last before TGGCGG)
       AGGAGG .... CTGTTG        (no CCA present, count all CGG)

Scans all anchor hits per read and keeps the one with the highest CCG count,
so random genomic matches (which give 0) don't mask the real HTT hit.
"""

import os
import csv
import re
import glob

RAW_DIR = "/scratch/SCWF00044/Ruban/Lymphoblastoids_proline_phasing/260510/PB045/"
OUTPUT_DIR = "/scratch/SCWF00044/Ruban/Lymphoblastoids_proline_phasing/260510/PB045/ccg_7_10_fasta"

os.makedirs(OUTPUT_DIR, exist_ok=True)


def count_ccg(seq):
    count = 0
    for i in range(0, len(seq) - 2, 3):
        if seq[i:i+3] == "CCG":
            count += 1
        else:
            break
    return count


def count_cgg(seq):
    count = 0
    for i in range(0, len(seq) - 2, 3):
        if seq[i:i+3] == "CGG":
            count += 1
        else:
            break
    return count


def find_all(seq, pattern):
    positions = []
    start = 0
    while True:
        pos = seq.find(pattern, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    return positions


def analyse_read(seq):
    seq = seq.upper()
    best_strand = "none"
    best_count = 0

    # forward: find all CAACAGCCGCCA hits, count CCG after each, keep best
    for pos in find_all(seq, "CAACAGCCGCCA"):
        ccg_start = pos + 12
        cctcct_pos = seq.find("CCTCCT", ccg_start)
        if cctcct_pos != -1:
            region = seq[ccg_start:cctcct_pos]
        else:
            tcct_pos = seq.find("TCCT", ccg_start)
            if tcct_pos != -1:
                region = seq[ccg_start:tcct_pos]
            else:
                continue
        c = count_ccg(region)
        if c > best_count:
            best_count = c
            best_strand = "fwd"

    # reverse: find all AGGAGGCGG hits, count CGG after each, keep best
    for pos in find_all(seq, "AGGAGGCGG"):
        cgg_start = pos + 6
        tggcgg_pos = seq.find("TGGCGG", cgg_start)
        if tggcgg_pos != -1:
            # TGGCGG = TGG(CCA) + CGG(CCG before CCA), both inside the anchor
            region = seq[cgg_start:tggcgg_pos]
        else:
            # no CCA, use CTGTTG or GAAG as boundary, all CGGs are real
            end_pos = -1
            for boundary in ["CTGTTG", "GAAG"]:
                p = seq.find(boundary, cgg_start)
                if p != -1 and (end_pos == -1 or p < end_pos):
                    end_pos = p
            if end_pos != -1:
                region = seq[cgg_start:end_pos]
            else:
                continue
        if len(region) < 3:
            continue
        c = count_cgg(region)
        if c > best_count:
            best_count = c
            best_strand = "rev"

    return best_strand, best_count


def parse_fasta(fasta_path):
    rid = None
    seq = []
    with open(fasta_path) as f:
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


def write_fasta(seqs, path):
    with open(path, 'w') as f:
        for rid, s in sorted(seqs, key=lambda x: x[0]):
            f.write(f">{rid}\n")
            for i in range(0, len(s), 80):
                f.write(s[i:i+80] + '\n')


def get_barcode(filename):
    m = re.search(r'(bc\d+)', filename)
    return m.group(1) if m else filename.rsplit('.', 1)[0]


def main():
    fastas = sorted(glob.glob(os.path.join(RAW_DIR, "*.fasta")))
    if not fastas:
        print(f"No FASTA files found in {RAW_DIR}")
        return

    print(f"Found {len(fastas)} FASTA files\n")
    print(f"{'Barcode':<12} {'Total':>7} {'HTT':>6} {'CCG7':>6} {'CCG10':>6}")
    print('-' * 42)

    summary = []

    for fasta_path in fastas:
        fname = os.path.basename(fasta_path)
        bc = get_barcode(fname)

        ccg7_reads = []
        ccg10_reads = []
        csv_rows = []
        total = 0
        htt_count = 0

        for rid, seq in parse_fasta(fasta_path):
            total += 1
            strand, ccg_count = analyse_read(seq)
            csv_rows.append([rid, len(seq), strand, ccg_count])

            if ccg_count > 0:
                htt_count += 1
            if ccg_count == 7:
                ccg7_reads.append((rid, seq))
            elif ccg_count == 10:
                ccg10_reads.append((rid, seq))

        csv_path = os.path.join(OUTPUT_DIR, f"{bc}_repeat_counts.csv")
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['read_id', 'seq_length', 'strand', 'ccg_count'])
            w.writerows(csv_rows)

        if ccg7_reads:
            write_fasta(ccg7_reads, os.path.join(OUTPUT_DIR, f"{bc}_CCG7.fasta"))
        if ccg10_reads:
            write_fasta(ccg10_reads, os.path.join(OUTPUT_DIR, f"{bc}_CCG10.fasta"))

        print(f"{bc:<12} {total:>7} {htt_count:>6} {len(ccg7_reads):>6} {len(ccg10_reads):>6}")
        summary.append([bc, total, htt_count, len(ccg7_reads), len(ccg10_reads)])

    with open(os.path.join(OUTPUT_DIR, 'extraction_summary.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['barcode', 'total_reads', 'htt_reads', 'ccg7', 'ccg10'])
        w.writerows(summary)

    print(f"\nDone. Output in {OUTPUT_DIR}")


if __name__ == '__main__':
    main()