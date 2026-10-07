"""What a per-locus time limit on the consensus stage would cost.

A locus is one consensus container (one job of the batch loop). Its time is
the wall clock from its 'PROGRESS ... Processing container' line to the next
container or the end of the batch, escalation retries included (they run
inline). A limit of T s drops every container that needed longer, and with
it every call whose consensuses all come from dropped containers.

Per variant and limit:
  containers containers over the limit, summed over the trio
  crs        their candidate regions
  calls      family VCF calls lost (all of a call's CONSENSUSIDs dropped)
  HG002 TP / FP   truvari bench (T2TQ100 all, before refine) TPs and FPs lost

usage: time_limits.py [svp_tiered15 dir] [--variants cr_base,cr_s1_sat] [--limits 30,60,120,180,300]
"""

from __future__ import annotations

import argparse
import glob
import json
import re
import sqlite3
import subprocess
from datetime import datetime

import pandas as pd

TS = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}) - \S+ - \w+ - (.*)$")
START = re.compile(r"PROGRESS \[\d+/\d+\] Processing container \(representative crID (\d+)\)")
SAMPLES = ("HG002", "HG003", "HG004")


def container_seconds(workdir: str) -> pd.Series:
    """Wall seconds per container (representative crID) of one run."""
    secs: dict[int, float] = {}
    for path in glob.glob(workdir + "/consensus/*/consensus.batch_*.log"):
        if path.endswith(".diag.log"):
            continue
        cur, t0 = None, None
        for line in open(path):
            m = TS.match(line)
            if not m:
                continue
            t = datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S,%f")
            msg = m.group(2)
            start = START.search(msg)
            end = msg.startswith("Wrote ") and "container results" in msg
            if (start or end) and cur is not None:
                secs[cur] = secs.get(cur, 0.0) + (t - t0).total_seconds()
                cur = None
            if start:
                cur, t0 = int(start.group(1)), t
    return pd.Series(secs, dtype=float)


def crs_per_container(workdir: str) -> pd.Series:
    """Number of candidate regions per container (representative crID)."""
    db = sqlite3.connect(workdir + "/crs_containers.db")
    return pd.Series({cid: len(json.loads(data)["crs"])
                      for cid, data in db.execute("select crID, data from containers")})


def call_sources(vcf: str, sample: str | None = None) -> list[set[tuple[str, int]]]:
    """Per call: the (sample, container) pairs of its CONSENSUSIDs."""
    out = subprocess.run(["bash", "-c", f"zcat {vcf} | grep -v '^#' | cut -f8"],
                         capture_output=True, text=True, check=True).stdout
    calls = []
    for info in out.splitlines():
        m = re.search(r"CONSENSUSIDs=([^;]+)", info)
        src = set()
        for cid in m.group(1).split(",") if m else ():
            smp, _, cons = cid.rpartition(":")
            smp = smp or sample
            src.add((smp, int(cons.split(".")[0])))  # 664.0, 1168.rescue
        calls.append(src)
    return calls


def lost(calls: list[set], dropped: set) -> int:
    return sum(1 for src in calls if src and src <= dropped)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", nargs="?", default="/home/mayv_c/development/svp_tiered15")
    ap.add_argument("--variants", default="cr_base,cr_s1_sat")
    ap.add_argument("--limits", default="30,60,120,180,300")
    a = ap.parse_args()
    limits = [float(x) for x in a.limits.split(",")]
    rows, dist = [], []
    for v in a.variants.split(","):
        base = f"{a.root}/results/{v}/20x"
        secs = {s: container_seconds(f"{base}/{s}/work") for s in SAMPLES}
        ncrs = {s: crs_per_container(f"{base}/{s}/work") for s in SAMPLES}
        n_crs = sum(int(x.sum()) for x in ncrs.values())
        family = call_sources(f"{base}/family.vcf.gz")
        tv = f"{base}/truvari/T2TQ100/all"
        tp, fp = call_sources(f"{tv}/tp-comp.vcf.gz"), call_sources(f"{tv}/fp.vcf.gz")
        n_loci = sum(len(x) for x in secs.values())
        allsec = pd.concat(secs.values())
        dist.append({"variant": v, "containers": n_loci, "crs": n_crs, "median_s": allsec.median(),
                     "p90_s": allsec.quantile(.9), "p99_s": allsec.quantile(.99),
                     "max_s": allsec.max(), "total_h": allsec.sum() / 3600})
        for lim in [float("inf")] + limits:
            dropped = {(s, c) for s, x in secs.items() for c in x.index[x > lim]}
            kept_h = sum(x.clip(upper=lim).sum() for x in secs.values()) / 3600
            rows.append({
                "variant": v, "limit_s": "none" if lim == float("inf") else int(lim),
                "containers_lost": len(dropped), "containers_lost_pct": 100 * len(dropped) / n_loci,
                "crs_lost": sum(int(ncrs[s].get(c, 0)) for s, c in dropped),
                "calls_lost": lost(family, dropped), "calls": len(family),
                "HG002_TP_lost": lost(tp, dropped), "HG002_TP": len(tp),
                "HG002_FP_lost": lost(fp, dropped), "HG002_FP": len(fp),
                "consensus_h": kept_h,
            })
    pd.set_option("display.width", 250)
    print(pd.DataFrame(dist).to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print()
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.2f}"))


if __name__ == "__main__":
    main()
