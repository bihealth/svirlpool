"""Inspect records with two different-size HG002 consensuses (child_collapse.json)."""

import json
import pickle
import random
import sys
from pathlib import Path

from svirlpool.localassembly import SVpatterns
from svirlpool.svcalling import core_alignment_merge as cam
from svirlpool.svcalling import svcomposite_merging as M

P = Path("/home/mayv_c/development/svp_mendel/proto")
d = Path(sys.argv[1])
n = int(sys.argv[2])
recs = json.load(open(d / "child_collapse.json"))
comps = pickle.load(open(P / "all_svComposites.pkl", "rb"))
cores = pickle.load(open(P / "cores.pkl", "rb"))
idx = {}
for s in comps:
    for p in s.svPatterns:
        idx[(p.samplename, p.consensusID, p.read_start, p.read_end)] = s
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
params = cam.CoreAlnParams()
for rec in random.sample(recs, n):
    pats = rec[3]
    members = [idx[(p[0], p[1], p[3], p[4])] for p in pats]
    chrom, pos = members[0].ref_start
    near = [s for s in comps if s.ref_start[0] == chrom and abs(s.ref_start[1] - pos) < 20000
            and not issubclass(s.sv_type, SVpatterns.SVpatternAdjacency)]
    group = None
    for cc in M._candidate_groups_for_chromosome(chrom, near, 150):
        g = [near[i] for i in cc]
        if any(m in g for m in members):
            group = g
            break
    log = []
    st = cam.GroupStats()
    new, rest = cam.merge_group_with_core_alignments(
        group, cores,
        lambda a, b: M.can_merge_svComposites_insertions(a, b, 0.1, 150, 0.7),
        lambda a, b: M.can_merge_svComposites_deletions(a, b, 0.1, 150, 0.7),
        "align", params, st, log,
    )
    print(f"\n### record {rec[0]} size {rec[1]} at {chrom}:{pos}")
    for m in members:
        print("   member", cam._pair_log_id(m), m.get_size())
    ids = {cam._pair_log_id(m) for m in members}
    for r in log:
        if r["a"] in ids and r["b"] in ids and r["final"]:
            print(f"   LINK {r['a']} {r['size_a']} -- {r['b']} {r['size_b']}: old={r['old']} aln={r['link_aln']} twin={r.get('twin')} "
                  f"corr={r.get('corr')} indel={r.get('indel')} reason={r.get('reason')}")
