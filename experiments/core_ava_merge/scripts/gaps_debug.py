"""Print the alignment gaps inside the judged window for HG002 vs others at a locus.

usage: gaps_debug.py chr pos [chr pos ...]
"""

import pickle
import sys
from pathlib import Path

from svirlpool.localassembly import SVpatterns
from svirlpool.svcalling import core_alignment_merge as cam
from svirlpool.svcalling import svcomposite_merging as M

P = Path("/home/mayv_c/development/svp_mendel/proto")
comps = pickle.load(open(P / "all_svComposites.pkl", "rb"))
cores = pickle.load(open(P / "cores.pkl", "rb"))
params = cam.CoreAlnParams()
args = sys.argv[1:]
for chrom, pos in zip(args[::2], map(int, args[1::2]), strict=True):
    near = [s for s in comps if s.ref_start[0] == chrom and abs(s.ref_start[1] - pos) < 5000
            and issubclass(s.sv_type, (SVpatterns.SVpatternInsertion, SVpatterns.SVpatternDeletion))]
    for cc in M._candidate_groups_for_chromosome(chrom, near, 150):
        g = [near[i] for i in cc]
        if not any(abs(s.ref_start[1] - pos) < 300 for s in g):
            continue
        placed, oriented = cam.place_composites(g, cores)
        spans, drift = cam.core_reference_spans(g, placed)
        keys = {p.key for p in placed if p is not None}
        alns = cam.align_group_cores(oriented, keys, params, cam.GroupStats(), {k: v[0] for k, v in spans.items()}, drift)
        child = [(s, p) for s, p in zip(g, placed, strict=True) if p is not None and p.key[0] == "HG002" and abs(s.ref_start[1] - pos) < 300]
        others = [(s, p) for s, p in zip(g, placed, strict=True) if p is not None and abs(s.ref_start[1] - pos) < 300]
        for xs, x in child:
            for ys, y in others:
                if y.key == x.key or ys.sv_type is not xs.sv_type:
                    continue
                aln = alns.get((x.key, y.key))
                if aln is None:
                    continue
                d = cam.decide_pair(x, y, aln, len(oriented[x.key]), len(oriented[y.key]), params,
                                    spans.get(x.key), spans.get(y.key), alns.get((y.key, x.key)))
                w = d.detail.get("window")
                print(f"{chrom}:{pos} {x.key} {xs.get_size()} vs {y.key} {ys.get_size()}: indel={d.detail.get('indel')} window={w}")
                if w:
                    gaps = [(gq, ql, tl) for gq, gt, ql, tl in aln.gaps if (gq < w[1] and gq + max(ql, 1) > w[0])]
                    print("    gaps (q, qlen, tlen):", gaps)
