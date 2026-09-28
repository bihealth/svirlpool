"""Where the consensus breaks with added read errors: phasing vs assembly.

usage: breakdown.py <run_dir>   (after eval_rates.py wrote <run_dir>/eval.tsv)
"""

import glob
import json
import sys
from pathlib import Path

import pandas as pd

run = Path(sys.argv[1])
df = pd.read_csv(run / "eval.tsv", sep="\t")
meta = []
for d in sorted(run.glob("rate_*")):
    rate = float(d.name.split("_")[1])
    for f in glob.glob(str(d / "chunk_*.jsonl")):
        for line in open(f):
            cd = json.loads(line)["consensus_dicts"]
            if not cd:
                continue
            c = next(iter(cd.values()))
            m = c.get("clustering_meta_data") or {}
            meta.append(
                {
                    "rate": rate,
                    "crID": int(c["ID"].split(".")[0]),
                    "n_snv": m.get("n_snv_sites", 0),
                    "n_sv": m.get("n_sv_sites", 0),
                    "n_lowq": m.get("n_low_quality", 0),
                    "n_unassigned": m.get("n_unassigned", 0),
                }
            )
df = df.merge(pd.DataFrame(meta), on=["rate", "crID"], how="left")
e = df[df.evaluable]
rows = []
for rate, g in df.groupby("rate"):
    ge = e[e.rate == rate]
    ph = g[g.method == "phased"]
    rows.append(
        {
            "rate": rate,
            "rec_het": round(ge[ge.cat == "het"].all_recovered.mean(), 3),
            "rec_cpx": round(ge[ge.cat == "cpx_het"].all_recovered.mean(), 3),
            "rec_hom": round(ge[ge.cat == "hom"].all_recovered.mean(), 3),
            "rec_ref": round(ge[ge.cat == "ref"].all_recovered.mean(), 3),
            "phased": len(ph),
            "acc_phased": round(ph.pair_acc.mean(), 4),
            "med_snv": g.n_snv.median(),
            "med_sv": g.n_sv.median(),
            "lowq/unassigned": f"{g.n_lowq.sum():.0f}/{g.n_unassigned.sum():.0f}",
            "status": dict(g.phasing_status.value_counts()),
        }
    )
pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 80)
print(pd.DataFrame(rows).to_string(index=False))
print("\ncontainers per category:", dict(e[e.rate == 0].cat.value_counts()))
