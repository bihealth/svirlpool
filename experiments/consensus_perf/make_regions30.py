"""Region set for the larger phasing ablation: random 5 Mb tiles covering
~30% of every autosome's non-N sequence, DISJOINT from the svp_improvements
5% blocks (+-1 Mb), so rules fitted on the 5% can be tested on it.

usage: make_regions30.py <reference.fa> <sweep_regions.v2.bed> <out.bed> [--frac 0.3] [--seed 30]
"""

import argparse
import random

import pysam

TILE = 5_000_000
PAD = 1_000_000


def main():
    p = argparse.ArgumentParser()
    p.add_argument("ref")
    p.add_argument("blocks")
    p.add_argument("out")
    p.add_argument("--frac", type=float, default=0.3)
    p.add_argument("--seed", type=int, default=30)
    a = p.parse_args()
    blocks = {}
    for line in open(a.blocks):
        c, s, e = line.split()[:3]
        blocks.setdefault(c, []).append((int(s) - PAD, int(e) + PAD))
    rng = random.Random(a.seed)
    fa = pysam.FastaFile(a.ref)
    out, tot_nonN, tot_sel = [], 0, 0
    for c in [f"chr{i}" for i in range(1, 23)]:
        L = fa.get_reference_length(c)
        tiles = []
        nonN = 0
        for s in range(0, L, TILE):
            e = min(L, s + TILE)
            seq = fa.fetch(c, s, e)
            n_ok = e - s - seq.count("N") - seq.count("n")
            nonN += n_ok
            if n_ok < 0.5 * (e - s):
                continue
            if any(s < be and bs < e for bs, be in blocks.get(c, [])):
                continue
            tiles.append((s, e, n_ok))
        rng.shuffle(tiles)
        sel, got = [], 0
        for s, e, n_ok in tiles:
            if got >= a.frac * nonN:
                break
            sel.append((s, e))
            got += n_ok
        tot_nonN += nonN
        tot_sel += got
        merged = []
        for s, e in sorted(sel):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = e
            else:
                merged.append([s, e])
        out += [(c, s, e) for s, e in merged]
    with open(a.out, "w") as f:
        for i, (c, s, e) in enumerate(out):
            f.write(f"{c}\t{s}\t{e}\tabl30_{i:03d}\n")
    print(f"{len(out)} blocks, {tot_sel / 1e6:.0f} Mb non-N selected = {tot_sel / tot_nonN:.3f} of autosomal non-N")


if __name__ == "__main__":
    main()
