"""Per-group minimap2 time with the module's settings, for one chromosome.

usage: mm2_profile.py CHROM [timeout]
"""

import pickle
import sys
import time
from pathlib import Path

from svirlpool.localassembly import SVpatterns
from svirlpool.svcalling import core_alignment_merge as cam
from svirlpool.svcalling import svcomposite_merging as M

P = Path("/home/mayv_c/development/svp_mendel/proto")
chrom = sys.argv[1]
timeout = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
cs, cores = pickle.load(open(P / f"chrom_{chrom}.pkl", "rb"))
params = cam.CoreAlnParams(timeout=timeout)
rows = []
for cc in M._candidate_groups_for_chromosome(chrom, cs, 150):
    g = [cs[i] for i in cc if issubclass(cs[i].sv_type, (SVpatterns.SVpatternInsertion, SVpatterns.SVpatternDeletion))]
    placed, oriented = cam.place_composites(g, cores)
    keys = {p.key for p in placed if p is not None}
    if len(keys) < 2:
        continue
    st = cam.GroupStats()
    t = time.perf_counter()
    alns = cam.align_group_cores(oriented, keys, params, st)
    dt = time.perf_counter() - t
    lens = sorted(len(oriented[k]) for k in keys)
    if dt > 2:
        d = P / "slow_groups"
        d.mkdir(exist_ok=True)
        with open(d / f"{chrom}_{g[0].ref_start[1]}.fa", "w") as f:
            for k in keys:
                f.write(f">{cam._core_name(k)}\n{oriented[k]}\n")
    rows.append((dt, len(keys), sum(lens), lens[-1], st.minimap2_timeouts, len(alns) // 2, g[0].ref_start))
rows.sort()
tot = sum(r[0] for r in rows)
print(f"{len(rows)} groups, total {tot:.1f}s")
acc = 0.0
for frac in (0.5, 0.9, 0.99):
    print(frac, rows[int(frac * len(rows))])
print("slowest:")
for r in rows[-15:]:
    print(r)
# share of time in groups with max core > 10 kb / total bp > 50 kb
for lim in (5000, 10000, 20000, 50000):
    sel = [r for r in rows if r[3] > lim]
    print(f"max core > {lim}: {len(sel)} groups, {sum(r[0] for r in sel):.1f}s")
