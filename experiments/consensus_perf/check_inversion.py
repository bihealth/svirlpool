"""Compare inverted alignments (parse_pair_both) with minimap2's own alignment
of the other direction.

usage: check_inversion.py <reads.fa> [...]
"""

import sys
import tempfile
from pathlib import Path

import numpy as np
from Bio import SeqIO

from svirlpool.localassembly import read_phasing as rp


def jacc(a, b):
    return len(a & b) / len(a | b) if a | b else 1.0


for path in sys.argv[1:]:
    reads = {r.id: r for r in SeqIO.parse(path, "fasta")}
    seqs = {n: str(r.seq) for n, r in reads.items()}
    params = rp.PhasingParams(align_pairs_once=False)
    with tempfile.TemporaryDirectory() as tmp:
        alns = rp.run_ava(reads, Path(tmp), params, threads=1, timeout=600)
    direct = {(a.query_name, a.reference_name): a for a in alns}
    stats = []
    for (q, t), a in direct.items():
        if (t, q) not in direct:
            continue
        fwd, inv = rp.parse_pair_both(a, seqs)
        ref = rp.parse_pair(direct[(t, q)], seqs)  # t aligned on q
        assert (inv.q, inv.t) == (ref.q, ref.t)
        shared = inv.mism.keys() & ref.mism.keys()
        stats.append(
            (
                jacc(set(inv.mism), set(ref.mism)),
                np.mean([inv.mism[k] == ref.mism[k] for k in shared])
                if shared
                else 1.0,
                jacc(inv.dels, ref.dels),
                abs(inv.t_start - ref.t_start) + abs(inv.t_end - ref.t_end),
                abs(inv.clip_left - ref.clip_left)
                + abs(inv.clip_right - ref.clip_right),
                abs(sum(s for _, s in inv.indels) - sum(s for _, s in ref.indels)),
                a.is_reverse,
            )
        )
        # the forward view must be exactly parse_pair
        f2 = rp.parse_pair(a, seqs)
        assert (fwd.mism, fwd.dels, fwd.indels) == (f2.mism, f2.dels, f2.indels)
    for strand in (0, 1):
        s = np.array([x for x in stats if x[6] == strand], dtype=float)
        if len(s) == 0:
            continue
        print(
            f"{Path(path).name} {'reverse' if strand else 'forward'}: pairs {len(s)}  "
            f"mism jaccard median {np.median(s[:, 0]):.3f} mean {s[:, 0].mean():.3f}  "
            f"base agreement {s[:, 1].mean():.4f}  "
            f"del jaccard median {np.median(s[:, 2]):.3f}  span diff median {np.median(s[:, 3]):.0f} "
            f"p90 {np.quantile(s[:, 3], 0.9):.0f}  clip diff median {np.median(s[:, 4]):.0f}  "
            f"net indel diff median {np.median(s[:, 5]):.0f} p90 {np.quantile(s[:, 5], 0.9):.0f}"
        )
