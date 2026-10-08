"""Would choosing the escalation level from the assembly size avoid retries?

Per lamassemble call of a run's consensus batch logs (README section 8):
input reads and bp, its level (threads, timeout), its wall seconds (to the next
log line) and whether it timed out. Per container: the level-1 time thrown
away when it escalated. Then a replay of the rule "start an assembly of more
than B bp at the second level": how many level-1 timeouts it predicts, how
many calls it moves needlessly, and the wall and thread-seconds of the
discarded level-1 attempts it would save.

usage: call_levels.py <work dir> [<work dir> ...] [--out calls.tsv]
"""

from __future__ import annotations

import argparse
import glob
import re
from datetime import datetime

import numpy as np
import pandas as pd

TS = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}) - (\S+) - (\w+) - (.*)$")
START = re.compile(r"PROGRESS \[\d+/\d+\] Processing container \(representative crID (\d+)\)")
LAM = re.compile(
    r"Running lamassemble on \S+ for (\S+) \(both strands: \w+\) with timeout of (\d+) seconds "
    r"\((\d+) reads, (\d+) bp\)"
)
ESC = re.compile(r"Container \d+: (.*) timed out \((\d+) thread\(s\), (\d+) s\); escalating")
LAM_LAST = re.compile(r"lamassemble timed out for")


def _t(s: str) -> float:
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S,%f").timestamp()


def parse(path: str, run: str) -> tuple[list[dict], list[dict]]:
    calls, escal = [], []
    crid, level, lam, level_start = None, 0, None, 0.0
    for line in open(path):
        m = TS.match(line)
        if not m:
            continue
        t, msg = _t(m.group(1)), m.group(4)
        if lam is not None:  # the call in progress ended with this line
            lam["seconds"] = t - lam.pop("t0")
            lam["timed_out"] = bool(ESC.search(msg) and "lamassemble" in msg) or bool(LAM_LAST.search(msg))
            calls.append(lam)
            lam = None
        if s := START.search(msg):
            crid, level, level_start = int(s.group(1)), 0, t
        elif lm := LAM.search(msg):
            lam = {"run": run, "crID": crid, "level": level, "timeout": int(lm.group(2)),
                   "reads": int(lm.group(3)), "bp": int(lm.group(4)), "t0": t}
        elif e := ESC.search(msg):
            escal.append({"run": run, "crID": crid, "level": level, "tool": e.group(1),
                          "threads": int(e.group(2)), "level_seconds": t - level_start})
            level, level_start = level + 1, t
    return calls, escal


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work", nargs="+")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    calls, escal = [], []
    for w in a.work:
        run = "/".join(w.rstrip("/").split("/")[-4:-1])
        for p in glob.glob(w + "/consensus/*/consensus.batch_*.log"):
            if not p.endswith(".diag.log"):
                c, e = parse(p, run)
                calls += c
                escal += e
    c, e = pd.DataFrame(calls), pd.DataFrame(escal)
    if a.out:
        c.to_csv(a.out, sep="\t", index=False)
        e.to_csv(a.out.replace(".tsv", ".escalations.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250)
    fmt = lambda x: f"{x:.3f}"  # noqa: E731
    l1 = c[c.level == 0]
    print(f"{len(c)} lamassemble calls, {len(l1)} at level 1; {int(l1.timed_out.sum())} of these timed out; "
          f"escalations: {e.tool.value_counts().to_dict()}")

    print("\n== level-1 lamassemble calls: timeout rate by input size")
    l1 = l1.assign(bin=pd.cut(l1.bp, [0, 2e4, 5e4, 1e5, 1.5e5, 2e5, 3e5, np.inf]))
    print(l1.groupby("bin", observed=True).agg(calls=("bp", "size"), timed_out=("timed_out", "sum"),
          rate=("timed_out", "mean"), median_s=("seconds", "median"), p99_s=("seconds", lambda x: x.quantile(.99)))
          .to_string(float_format=fmt))
    from sklearn.metrics import roc_auc_score

    print(f"AUC of the bp for a level-1 timeout: {roc_auc_score(l1.timed_out, l1.bp):.3f}; "
          f"of the reads: {roc_auc_score(l1.timed_out, l1.reads):.3f}")

    l2 = c[c.level == 1]
    if len(l2) and l2.timed_out.any():
        l2 = l2.assign(bin=pd.cut(l2.bp, [0, 5e4, 1e5, 2e5, 3e5, np.inf]))
        print("\n== level-2 lamassemble calls (4 threads, 60 s): timeout rate by input size")
        print(l2.groupby("bin", observed=True).agg(calls=("bp", "size"), timed_out=("timed_out", "sum"),
              rate=("timed_out", "mean"), median_s=("seconds", "median")).to_string(float_format=fmt))

    # replay: start an assembly of > B bp at level 2. The level-1 attempts of the
    # containers whose level-1 lamassemble timeout came from such a call are saved;
    # level-1 calls of > B bp that finished are moved to 4 threads needlessly.
    lam_esc = e[(e.tool == "lamassemble") & (e.level == 0)].set_index(["run", "crID"])
    trig = l1[l1.timed_out].set_index(["run", "crID"])
    print("\n== replay: start assemblies of more than B bp at level 2 (4 threads, 60 s)")
    rows = []
    for kb in (20, 30, 50, 75, 100, 150, 200):
        hit = trig[trig.bp > kb * 1e3]
        saved = lam_esc.loc[lam_esc.index.intersection(hit.index)]
        moved = l1[(l1.bp > kb * 1e3) & ~l1.timed_out]
        rows.append({"B_kb": kb, "timeouts_predicted": len(hit), "of": len(trig),
                     "level1_wall_saved_h": saved.level_seconds.sum() / 3600,
                     "needless_moves": len(moved), "their_level1_s": moved.seconds.sum() / 3600})
    print(pd.DataFrame(rows).to_string(index=False, float_format=fmt))
    print(f"(all discarded level-1 time of lamassemble escalations: {lam_esc.level_seconds.sum() / 3600:.3f} h; "
          f"of phasing escalations: {e[(e.tool != 'lamassemble') & (e.level == 0)].level_seconds.sum() / 3600:.3f} h)")


if __name__ == "__main__":
    main()
