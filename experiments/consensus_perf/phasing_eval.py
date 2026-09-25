"""Time ``read_phasing.phase_reads`` and score it against the trio read labels.

Reads the per-container read sets written by ``dump_phasing_reads.py``.

usage: phasing_eval.py <reads_dir> <out.tsv> [--procs 12] [--crIDs FILE|a,b,c]
                       [--compare BASE.tsv] [param=value ...]

``param=value`` sets ``read_phasing.PhasingParams`` fields. Prints a summary:
time, trio pair accuracy of the assigned reads, fraction of containers phased
perfectly, fraction of trio-labelled reads assigned, allele count histogram
and, with --compare, how many containers changed their partition.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from Bio import SeqIO

from svirlpool.localassembly import read_phasing

logging.disable(logging.WARNING)
TRIO = Path(__file__).parent.parent / "ava_phasing" / "data" / "trio_read_labels.tsv.gz"


def load_trio() -> dict[int, dict[str, str]]:
    out: dict[int, dict[str, str]] = {}
    with gzip.open(TRIO, "rt") as f:
        h = f.readline().rstrip("\n").split("\t")
        ic, ir, il = h.index("crID"), h.index("read_name"), h.index("label")
        for line in f:
            x = line.rstrip("\n").split("\t")
            if x[il] in ("pat", "mat"):
                out.setdefault(int(x[ic]), {})[x[ir]] = x[il]
    return out


def pair_accuracy(groups: dict[str, int], trio: dict[str, str]):
    lab = [(g, trio[r]) for r, g in groups.items() if r in trio]
    m = len(lab)
    if m < 2:
        return m, float("nan")
    ok = sum(
        (lab[i][0] == lab[j][0]) == (lab[i][1] == lab[j][1])
        for i in range(m)
        for j in range(i + 1, m)
    )
    return m, ok / (m * (m - 1) / 2)


def run_one(args):
    path, params = args
    reads = {r.id: r for r in SeqIO.parse(path, "fasta")}
    t0 = time.perf_counter()
    res = read_phasing.phase_reads(reads, params=params, timeout=20)
    return {
        "crID": int(path.stem),
        "n_reads": len(reads),
        "status": res.status,
        "n_alleles": res.n_alleles,
        "n_snv": res.n_snv_sites,
        "n_sv": res.n_sv_sites,
        "n_lowq": len(res.low_quality),
        "secs": round(time.perf_counter() - t0, 3),
        "groups": json.dumps(res.groups, sort_keys=True),
    }


def ari(a: dict[str, int], b: dict[str, int]) -> float:
    from sklearn.metrics import adjusted_rand_score

    common = sorted(a.keys() & b.keys())
    if len(common) < 2:
        return float("nan")
    return float(adjusted_rand_score([a[x] for x in common], [b[x] for x in common]))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("reads_dir", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--procs", type=int, default=12)
    p.add_argument("--crIDs")
    p.add_argument("--compare", type=Path)
    p.add_argument("params", nargs="*")
    a = p.parse_args()
    kw = {}
    for x in a.params:
        k, v = x.split("=", 1)
        kw[k] = eval(v)
    params = read_phasing.PhasingParams(**kw)
    files = sorted(a.reads_dir.glob("*.fa"), key=lambda f: int(f.stem))
    if a.crIDs:
        src = Path(a.crIDs)
        ids = {
            int(x)
            for x in (src.read_text().split() if src.exists() else a.crIDs.split(","))
        }
        files = [f for f in files if int(f.stem) in ids]
    # largest first: better load balance
    files.sort(key=lambda f: -f.stat().st_size)
    t0 = time.perf_counter()
    with Pool(a.procs) as pool:
        rows = list(pool.imap_unordered(run_one, [(f, params) for f in files]))
    wall = time.perf_counter() - t0
    rows.sort(key=lambda r: r["crID"])
    trio = load_trio()
    for r in rows:
        m, acc = pair_accuracy(json.loads(r["groups"]), trio.get(r["crID"], {}))
        r["trio_assigned"] = m
        r["trio_total"] = len(trio.get(r["crID"], {}))
        r["acc"] = round(acc, 4) if acc == acc else ""
    cols = list(rows[0].keys())
    cols.remove("groups")
    cols.append("groups")
    with open(a.out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    secs = np.array([r["secs"] for r in rows])
    accs = np.array([r["acc"] for r in rows if r["acc"] != ""])
    print(f"params: {kw}")
    print(
        f"containers {len(rows)}  sum secs {secs.sum():.1f}  wall {wall:.1f} ({a.procs} procs)  "
        f"median {np.median(secs):.2f}  p99 {np.quantile(secs, 0.99):.1f}  max {secs.max():.1f}"
    )
    print(
        f"pair acc (assigned) mean {accs.mean():.4f}  perfect {np.mean(accs == 1):.3f}  "
        f"trio reads assigned {sum(r['trio_assigned'] for r in rows) / max(1, sum(r['trio_total'] for r in rows)):.3f}"
    )
    print("status", dict(Counter(r["status"] for r in rows)))
    print("n_alleles", dict(sorted(Counter(r["n_alleles"] for r in rows).items())))
    if a.compare:
        base = {}
        with open(a.compare) as f:
            h = f.readline().rstrip("\n").split("\t")
            for line in f:
                x = dict(zip(h, line.rstrip("\n").split("\t"), strict=True))
                base[int(x["crID"])] = x
        same = diff_n = 0
        aris = []
        for r in rows:
            b = base.get(r["crID"])
            if b is None:
                continue
            if b["groups"] == r["groups"]:
                same += 1
            ga, gb = json.loads(r["groups"]), json.loads(b["groups"])
            v = ari(ga, gb)
            if v == v:
                aris.append(v)
            diff_n += int(b["n_alleles"]) != r["n_alleles"]
        bsecs = sum(float(b["secs"]) for b in base.values())
        print(
            f"vs {a.compare.name}: identical partitions {same}/{len(rows)}  "
            f"n_alleles differ {diff_n}  mean ARI {np.mean(aris):.4f}  "
            f"time {secs.sum():.1f} vs {bsecs:.1f} ({secs.sum() / bsecs:.2f}x)"
        )


if __name__ == "__main__":
    main()
