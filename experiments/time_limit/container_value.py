"""Are the slow consensus containers the ones that resolve badly?

Per HG002 container (representative crID) of one svp_* run:
  seconds   wall clock from its PROGRESS line to the next container / the end
            of the batch (escalation retries included)
  last      a tool still timed out at the last escalation level (degraded)
  dropped   dropped at --container-time-limit
  repr      Q100 representation: per haplotype the best identity of any of
            its core consensuses (consensus_q100.py), mean of the two
            haplotypes; NaN when no consensus was located
  tp / fp   calls of the HG002 VCF that truvari (T2TQ100, all regions, before
            refine) counts as TP / FP and that come (partly) from it

Prints the containers binned by seconds, and what a limit of each --limits
value would drop: containers, hours, their representation, TPs, FPs.

usage: container_value.py <svp root> <variant> [--limits 120,180,300] [--truth T2TQ100]
"""

from __future__ import annotations

import argparse
import glob
import re
import subprocess
from datetime import datetime

import numpy as np
import pandas as pd

TS = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}) - \S+ - \w+ - (.*)$")
START = re.compile(r"PROGRESS \[\d+/\d+\] Processing container \(representative crID (\d+)\)")
LAST = re.compile(r"Container (\d+): .* timed out at the last level")
DROPPED = re.compile(r"Container (\d+): dropped at the container time limit")
HAPS = ("MATERNAL", "PATERNAL")
BINS = [0, 20, 60, 120, 180, 300, 600, np.inf]


def containers(work: str) -> pd.DataFrame:
    secs: dict[int, float] = {}
    last, dropped = set(), set()
    for path in glob.glob(work + "/consensus/*/consensus.batch_*.log"):
        if path.endswith(".diag.log"):
            continue
        cur, t0 = None, None
        for line in open(path):
            m = TS.match(line)
            if not m:
                continue
            msg = m.group(2)
            if x := LAST.search(msg):
                last.add(int(x.group(1)))
            if x := DROPPED.search(msg):
                dropped.add(int(x.group(1)))
            start = START.search(msg)
            end = msg.startswith("Wrote ") and "container results" in msg
            if (start or end) and cur is not None:
                t = datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S,%f")
                secs[cur] = secs.get(cur, 0.0) + (t - t0).total_seconds()
                cur = None
            if start:
                cur = int(start.group(1))
                t0 = datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S,%f")
    c = pd.DataFrame({"seconds": pd.Series(secs, dtype=float)})
    c["last"] = c.index.isin(last)
    c["dropped"] = c.index.isin(dropped)
    return c


def representation(tsv: str) -> pd.Series:
    d = pd.read_csv(tsv, sep="\t")
    ok = d[d.status == "ok"].copy()
    for h in HAPS:
        ok[f"id_{h}"] = 1 - ok[f"ed_{h}"] / ok.len.clip(lower=1)
    g = ok.groupby("crID")
    best = pd.DataFrame({h: g[f"id_{h}"].max() for h in HAPS})
    return best.mean(axis=1, skipna=True)


def call_containers(vcf: str) -> list[set[int]]:
    out = subprocess.run(["bash", "-c", f"zcat {vcf} | grep -v '^#' | cut -f8"],
                         capture_output=True, text=True, check=True).stdout
    calls = []
    for info in out.splitlines():
        m = re.search(r"CONSENSUSIDs=([^;]+)", info)
        calls.append({int(x.rpartition(":")[2].split(".")[0]) for x in m.group(1).split(",")}
                     if m else set())
    return calls


def per_container(calls: list[set[int]]) -> pd.Series:
    n: dict[int, int] = {}
    for src in calls:
        for c in src:
            n[c] = n.get(c, 0) + 1
    return pd.Series(n, dtype=int)


def summarize(g: pd.DataFrame) -> dict:
    return {
        "containers": len(g), "hours": g.seconds.sum() / 3600,
        "last_level": int(g["last"].sum()), "dropped": int(g.dropped.sum()),
        "located": int(g.repr.notna().sum()), "repr_mean": g.repr.mean(),
        "repr_ge_0.99": (g.repr >= 0.99).sum() / max(1, g.repr.notna().sum()),
        "tp": int(g.tp.sum()), "fp": int(g.fp.sum()),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("variant")
    ap.add_argument("--limits", default="120,180,300")
    ap.add_argument("--truth", default="T2TQ100")
    a = ap.parse_args()
    base = f"{a.root}/results/{a.variant}/20x"
    c = containers(f"{base}/HG002/work")
    c["repr"] = representation(f"{base}/HG002/consensus_q100.tsv").reindex(c.index)
    tv = f"{base}/truvari/{a.truth}/all"
    c["tp"] = per_container(call_containers(f"{tv}/tp-comp.vcf.gz")).reindex(c.index).fillna(0)
    c["fp"] = per_container(call_containers(f"{tv}/fp.vcf.gz")).reindex(c.index).fillna(0)
    pd.set_option("display.width", 250)
    fmt = lambda x: f"{x:.3f}"  # noqa: E731
    c["bin"] = pd.cut(c.seconds, BINS, right=False)
    rows = [dict(seconds=str(b), **summarize(g)) for b, g in c.groupby("bin", observed=True)]
    rows.append(dict(seconds="all", **summarize(c)))
    print(f"== {a.variant}: HG002 containers by wall seconds")
    print(pd.DataFrame(rows).to_string(index=False, float_format=fmt))
    rows = []
    for lim in (float(x) for x in a.limits.split(",")):
        rows.append(dict(limit_s=int(lim), **summarize(c[c.seconds > lim]),
                         hours_saved=(c.seconds - c.seconds.clip(upper=lim)).sum() / 3600))
    print(f"\n== {a.variant}: containers over each limit (what the limit would drop)")
    print(pd.DataFrame(rows).to_string(index=False, float_format=fmt))


if __name__ == "__main__":
    main()
