"""Representative 15% region set for the end-to-end benchmark of the tiered
phasing: 5 Mb tiles, 15% of every autosome's non-N sequence, DISJOINT from the
svp_improvements 5% blocks (+-1 Mb, where the tiered rule was fitted).

Among --draws random selections, the one whose aggregate statistics deviate
least from the whole autosomal genome (all tiles, the 5% blocks included) is
kept. Statistics, per non-N Mb or as fractions of the non-N sequence: tandem
repeat (TRF) share, GC, T2TQ100 and V5 SVs (>= 50 bp) per Mb, and the share
covered by each truth set's benchmark regions. Score: the largest relative
deviation over the statistics.

usage: make_regions15.py <reference.fa> <sweep_regions.v2.bed> <trf.bed>
       <t2tq100.vcf.gz> <t2tq100.bed> <v5.vcf.gz> <v5.bed> <out.bed>
       [--frac 0.15] [--draws 5000] [--seed 15]
"""

import argparse
import random
from collections import defaultdict

import numpy as np
import pysam

TILE = 5_000_000
PAD = 1_000_000
CHROMS = [f"chr{i}" for i in range(1, 23)]


def bed_cover(path):
    """chrom -> sorted, merged intervals."""
    iv = defaultdict(list)
    for line in open(path):
        if line.startswith(("#", "track")):
            continue
        c, s, e = line.split()[:3]
        iv[c].append((int(s), int(e)))
    out = {}
    for c, v in iv.items():
        m = []
        for s, e in sorted(v):
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[c] = np.array(m, dtype=np.int64)
    return out


def covered(iv, s, e):
    if iv is None or len(iv) == 0:
        return 0
    a = np.clip(iv[:, 0], s, e)
    b = np.clip(iv[:, 1], s, e)
    return int((b - a).clip(min=0).sum())


def sv_positions(path):
    pos = defaultdict(list)
    for rec in pysam.VariantFile(path):
        alt = rec.alts[0] if rec.alts else ""
        svlen = rec.info.get("SVLEN")
        if isinstance(svlen, tuple):
            svlen = svlen[0]
        if svlen is None:
            svlen = len(alt) - len(rec.ref) if not alt.startswith("<") else 0
        if abs(int(svlen)) >= 50:
            pos[rec.chrom].append(rec.pos)
    return {c: np.sort(np.array(v)) for c, v in pos.items()}


def main():
    p = argparse.ArgumentParser()
    for k in ("ref", "blocks", "trf", "t2t_vcf", "t2t_bed", "v5_vcf", "v5_bed", "out"):
        p.add_argument(k)
    p.add_argument("--frac", type=float, default=0.15)
    p.add_argument("--draws", type=int, default=5000)
    p.add_argument("--seed", type=int, default=15)
    a = p.parse_args()
    blocks = defaultdict(list)
    for line in open(a.blocks):
        c, s, e = line.split()[:3]
        blocks[c].append((int(s) - PAD, int(e) + PAD))
    trf, t2t_bed, v5_bed = bed_cover(a.trf), bed_cover(a.t2t_bed), bed_cover(a.v5_bed)
    t2t, v5 = sv_positions(a.t2t_vcf), sv_positions(a.v5_vcf)
    fa = pysam.FastaFile(a.ref)
    tiles = []  # chrom, start, end, eligible, stats
    for c in CHROMS:
        L = fa.get_reference_length(c)
        for s in range(0, L, TILE):
            e = min(L, s + TILE)
            seq = fa.fetch(c, s, e).upper()
            n_ok = e - s - seq.count("N")
            if n_ok < 0.5 * (e - s):
                continue
            st = np.array(
                [
                    n_ok,
                    covered(trf.get(c), s, e),
                    seq.count("G") + seq.count("C"),
                    np.searchsorted(t2t.get(c, []), e)
                    - np.searchsorted(t2t.get(c, []), s),
                    np.searchsorted(v5.get(c, []), e)
                    - np.searchsorted(v5.get(c, []), s),
                    covered(t2t_bed.get(c), s, e),
                    covered(v5_bed.get(c), s, e),
                ],
                dtype=np.float64,
            )
            ok = not any(s < be and bs < e for bs, be in blocks.get(c, []))
            tiles.append((c, s, e, ok, st))
    names = ["trf", "gc", "t2t_sv_per_mb", "v5_sv_per_mb", "t2t_bed", "v5_bed"]

    def summary(stats):
        n = stats[:, 0].sum()
        return np.array(
            [
                stats[:, 1].sum() / n,
                stats[:, 2].sum() / n,
                stats[:, 3].sum() / n * 1e6,
                stats[:, 4].sum() / n * 1e6,
                stats[:, 5].sum() / n,
                stats[:, 6].sum() / n,
            ]
        )

    genome = summary(np.array([t[4] for t in tiles]))
    by_c = defaultdict(list)
    for i, t in enumerate(tiles):
        if t[3]:
            by_c[t[0]].append(i)
    need = {c: a.frac * sum(t[4][0] for t in tiles if t[0] == c) for c in CHROMS}
    rng = random.Random(a.seed)
    best = None
    for _ in range(a.draws):
        sel = []
        for c in CHROMS:
            idx = by_c[c][:]
            rng.shuffle(idx)
            got = 0.0
            for i in idx:
                n_i = tiles[i][4][0]
                # stop at the tile count closest to the target share
                if abs(got + n_i - need[c]) >= abs(got - need[c]):
                    break
                sel.append(i)
                got += n_i
        s = summary(np.array([tiles[i][4] for i in sel]))
        dev = np.abs(s / genome - 1)
        score = dev.max()
        if best is None or score < best[0]:
            best = (score, sel, s)
    score, sel, s = best
    sel = sorted(sel, key=lambda i: (CHROMS.index(tiles[i][0]), tiles[i][1]))
    merged = []
    for i in sel:
        c, st, e = tiles[i][:3]
        if merged and merged[-1][0] == c and st <= merged[-1][2]:
            merged[-1][2] = e
        else:
            merged.append([c, st, e])
    with open(a.out, "w") as f:
        for k, (c, st, e) in enumerate(merged):
            f.write(f"{c}\t{st}\t{e}\ttier15_{k:03d}\n")
    tot = sum(t[4][0] for t in tiles)
    got = sum(tiles[i][4][0] for i in sel)
    print(
        f"{len(sel)} tiles in {len(merged)} blocks, {got / 1e6:.0f} Mb non-N = {got / tot:.3f} of autosomal non-N"
    )
    print(f"best of {a.draws} draws, max relative deviation {score:.4f}")
    print(f"{'stat':16s} {'genome':>10s} {'subset':>10s} {'rel.dev':>8s}")
    for n, g, v in zip(names, genome, s, strict=True):
        print(f"{n:16s} {g:10.4f} {v:10.4f} {v / g - 1:+8.4f}")


if __name__ == "__main__":
    main()
