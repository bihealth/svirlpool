"""Run the AVA phasing prototype over many containers and tabulate.

usage: batch_eval.py <crIDs file|comma list> <out.tsv> [flank] [key=value ...]
"""

from __future__ import annotations

import gzip
import json
import logging
import sys
import time
import traceback
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

logging.disable(logging.INFO)
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import proto  # noqa: E402

TMP = Path("/home/mayv_c/.claude/jobs/bd7213cd/tmp/batch")
CUR = json.load(open(HERE / "data/flank200/current_clusters.json"))


def load_trio(path=HERE / "data/trio_read_labels.tsv.gz") -> dict[int, dict[str, str]]:
    if not path.exists():
        return {}
    out: dict[int, dict[str, str]] = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        ic, ir, il = header.index("crID"), header.index("read_name"), header.index("label")
        for line in f:
            x = line.rstrip("\n").split("\t")
            out.setdefault(int(x[ic]), {})[x[ir]] = x[il]
    return out


TRIO = load_trio()


def ari(a: list, b: list) -> float:
    from sklearn.metrics import adjusted_rand_score

    return float(adjusted_rand_score(a, b)) if len(a) > 1 else float("nan")


def partition_vs_trio(names, labels, trio: dict[str, str]):
    """Pair-level accuracy of a read partition against trio labels.

    Returns (n_labelled, pair_acc, same_hap_split, diff_hap_joined):
    over pairs of trio-labelled reads, fraction where 'same cluster' equals
    'same parent'; plus the two error kinds as fractions of their pair class.
    """
    lab = [(labels[i], trio[n]) for i, n in enumerate(names) if trio.get(n) in ("pat", "mat")]
    m = len(lab)
    if m < 2:
        return m, float("nan"), float("nan"), float("nan")
    ok = same_split = same_tot = diff_join = diff_tot = 0
    for i in range(m):
        for j in range(i + 1, m):
            sc = lab[i][0] == lab[j][0]
            sh = lab[i][1] == lab[j][1]
            ok += sc == sh
            if sh:
                same_tot += 1
                same_split += not sc
            else:
                diff_tot += 1
                diff_join += sc
    tot = m * (m - 1) // 2
    return (m, ok / tot, same_split / same_tot if same_tot else float("nan"),
            diff_join / diff_tot if diff_tot else float("nan"))


def run_one(args):
    crID, flank, pkw, kw = args
    t0 = time.time()
    try:
        r = proto.analyse(crID, flank, TMP / f"{crID}_{flank}", **dict(pkw))
    except Exception:
        return {"crID": crID, "error": traceback.format_exc(limit=1).splitlines()[-1]}
    names, A, D = r["names"], r["A"], r["D"]
    n = len(names)
    iu = np.triu_indices(n, 1)
    W = proto.phase_weights(A, D)
    labels = proto.correlation_cluster(W)
    if kw_get(kw, "refine", True):
        labels, meta = proto.refine_clusters(labels, names, r["sites"], W,
                                             min_disc=kw_get(kw, "min_disc", 2))
        labels = [l if l >= 0 else 10_000 + i for i, l in enumerate(labels)]
    sizes = sorted(Counter(labels).values(), reverse=True)
    min_group = kw_get(kw, "min_group", 3)
    big = [s for s in sizes if s >= min_group]
    row = {
        "crID": crID,
        "flank": flank,
        "n_reads": n,
        "raw_sites": len(r["raw_sites"]),
        "sites": len(r["sites"]),
        "sv_sites": len(r["sv_sites"]),
        "lowq": len(r["lowq"]),
        "pairs_decided": round(float(((A + D)[iu] > 0).mean()) if n > 1 else 0.0, 3),
        "n_groups": len(big),
        "grouped_frac": round(sum(big) / n, 3) if n else 0,
        "sizes": ",".join(map(str, sizes)),
        "secs": round(time.time() - t0, 1),
    }
    # reads in groups < min_group are unassigned (-1): they are re-added to a
    # consensus by alignment later, as the module already does
    cnt = Counter(labels)
    alabels = [l if cnt[l] >= min_group else -1 for l in labels]
    trio = TRIO.get(crID, {})
    asg = [i for i, l in enumerate(alabels) if l >= 0]
    m, acc, ss, dj = partition_vs_trio([names[i] for i in asg], [alabels[i] for i in asg], trio)
    lab_trio = [nm for nm in names if trio.get(nm) in ("pat", "mat")]
    row.update(assigned_trio_frac=round(m / len(lab_trio), 3) if lab_trio else float("nan"),
               acc_asg=round(acc, 3), split_err_asg=round(ss, 3), join_err_asg=round(dj, 3))
    # current method
    cur = CUR.get(str(crID), {})
    cur_lab = {rn: k for k, (cid, rns) in enumerate(sorted(cur.items())) for rn in rns}
    row["cur_n"] = len(cur)
    common = [i for i, nm in enumerate(names) if nm in cur_lab]
    row["ari_vs_cur"] = round(ari([labels[i] for i in common], [cur_lab[names[i]] for i in common]), 3) if common else float("nan")
    # trio
    trio = TRIO.get(crID, {})
    m, acc, ss, dj = partition_vs_trio(names, labels, trio)
    row.update(trio_n=m, trio_pat=sum(v == "pat" for v in trio.values()),
               trio_mat=sum(v == "mat" for v in trio.values()),
               acc_new=round(acc, 3), split_err_new=round(ss, 3), join_err_new=round(dj, 3))
    cn = [names[i] for i in common]
    m2, acc2, ss2, dj2 = partition_vs_trio(cn, [cur_lab[x] for x in cn], trio)
    row.update(acc_cur=round(acc2, 3), split_err_cur=round(ss2, 3), join_err_cur=round(dj2, 3))
    # per-read labels for later analysis
    row["_labels"] = json.dumps(dict(zip(names, map(int, alabels))))
    return row


def kw_get(kw, k, d):
    return dict(kw).get(k, d)


def main():
    src = sys.argv[1]
    crIDs = ([int(x) for x in open(src).read().split()] if Path(src).exists()
             else [int(x) for x in src.split(",")])
    out = Path(sys.argv[2])
    flank = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
    kw = []
    for a in sys.argv[4:]:
        k, v = a.split("=")
        kw.append((k, eval(v)))
    proto_kw = tuple((k, v) for k, v in kw if k not in ("min_group", "refine", "min_disc"))
    eval_kw = tuple((k, v) for k, v in kw if k in ("min_group", "refine", "min_disc"))
    rows = []
    with Pool(8) as pool:
        for row in pool.imap_unordered(run_one, [(c, flank, proto_kw, eval_kw) for c in crIDs]):
            rows.append(row)
    rows.sort(key=lambda r: r["crID"])
    cols = sorted({k for r in rows for k in r}, key=lambda k: (k.startswith("_"), list(rows[0].keys()).index(k) if k in rows[0] else 99))
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"wrote {len(rows)} rows to {out}")


if __name__ == "__main__":
    main()
