"""Which part of the rule loses the 30% T2TQ100 SVs, and what are they?"""

import time

import cr_validate as V
import crtools as T
import numpy as np
import pandas as pd

from svirlpool.candidateregions import signalstrength_to_crs as S2C

pd.set_option("display.width", 250)
work = V.RUN30 + "HG002/work/"
truth = V.truth30()
cur = T.load_crs(work + "crs.db")
base = T.covered(truth, cur)
big = truth["size"].values >= 1000
rows = []
for b, s in ((600, 0), (1200, 1), (1200, 2), (600, 1), (600, 2)):
    crs = V.form("HG002_30", work, buffer=b, support=s)
    now = T.covered(truth, crs)
    gone = base & ~now
    rows.append(dict(**V.cost_row(f"b{b} support{s}", crs), lost=int(gone.sum()), lost_1kb=int((gone & big).sum())))
print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.1f}"))

# the lost SVs of support2 b600: their signals in the current CRs and those signals' support
crs = V.form("HG002_30", work, buffer=600, support=2)
gone = base & ~T.covered(truth, crs)
sig = T.load_signals(work + "signalstrength.tsv.gz")
t = time.time()
support = np.zeros(len(sig), dtype=np.int32)
for _c, g in sig.groupby("chr"):
    support[g.index.values] = S2C.signal_support(
        starts=g.start.values, sizes=g["size"].values, sv_types=g["type"].values,
        repeat_ids=g.rep.values, reads=pd.factorize(g.read.values)[0])
print(f"support of {len(sig)} signals in {time.time() - t:.0f} s")
sig["support"] = support
per_chr = sig.groupby("chr").size().sort_values(ascending=False)
print("signals per chr (top):", per_chr.head(5).to_dict())
rep_sizes = sig[sig.rep >= 0].groupby(["chr", "rep"]).size().sort_values(ascending=False)
print("largest repeat groups:", rep_sizes.head(5).to_dict())
trf = pd.read_csv(T.TRF, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
out = []
for r in truth[gone].itertuples():
    near = sig[(sig.chr == r.chr) & (sig.end >= r.start - 200) & (sig.start <= r.end + 200)]
    near = near[np.abs(near["size"]) >= 30]
    tr = trf[(trf.c == r.chr) & (trf.e > r.start) & (trf.s < r.end + 1)]
    out.append({"chr": r.chr, "pos": r.start, "size": r.size, "n_sig": len(near), "n_reads": near.read.nunique(),
                    "max_support": int(near.support.max()) if len(near) else -1,
                    "max_strength": float(near.strength.max()) if len(near) else 0.0,
                    "sizes": ",".join(str(abs(x)) for x in sorted(near["size"].values)[:8]),
                    "in_tr": len(tr) > 0, "cov": int(near["cov"].median()) if len(near) else 0})
out = pd.DataFrame(out)
print(out.to_string(index=False))
print(out.describe().round(1).to_string())
