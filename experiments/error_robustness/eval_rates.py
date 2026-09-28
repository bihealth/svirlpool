"""Evaluate the consensus runs of ``run_rates.py`` against the truth.

Per rate and container:

* phasing: method / status of the consensuses, their number, and the trio
  pair accuracy of the reads' consensus membership (HG002 trio read labels);
* alleles: every (core) consensus is aligned to the reference around its
  locus (minimap2 map-ont, --eqx); its net indel over the locus span is
  matched to the T2TQ100 haplotype nets of ``container_truth.tsv`` like
  ``ava_phasing/container_truth.py`` does (both ~ref, or same sign and size
  ratio >= 0.7, one consensus per truth allele); a container is "recovered"
  when every truth allele is, "extra" when it has more distinct consensus
  alleles than truth alleles;
* consensus accuracy: mismatches + indels < 5 bp per aligned reference bp;
* cost: CPU seconds of the container (all attempts), escalations, unresolved.

usage: eval_rates.py <run_dir> <out.tsv> [--procs 12]
"""

from __future__ import annotations

import argparse
import glob
import gzip
import itertools
import json
import subprocess
import tempfile
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

HERE = Path(__file__).parent
DATA = HERE.parent / "ava_phasing" / "data"
REFERENCE = "/home/mayv_c/biodata/local/references/GRCh38/GRCh38.fa"
MIN_SV = 20
MIN_CONS_INDEL = 5
WINDOW = 10_000


def load_trio() -> dict[int, dict[str, str]]:
    out: dict[int, dict[str, str]] = {}
    with gzip.open(DATA / "trio_read_labels.tsv.gz", "rt") as f:
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
        return float("nan")
    ok = sum(
        (lab[i][0] == lab[j][0]) == (lab[i][1] == lab[j][1])
        for i in range(m)
        for j in range(i + 1, m)
    )
    return ok / (m * (m - 1) / 2)


def allele_match(a, b, r=0.7, ref_tol=MIN_SV):
    if abs(a) < ref_tol and abs(b) < ref_tol:
        return True
    if a == 0 or b == 0 or (a > 0) != (b > 0):
        return False
    return min(abs(a), abs(b)) / max(abs(a), abs(b)) >= r


def cons_stats(segs, lo, hi):
    """net indel over [lo, hi] and small-difference counts of one consensus."""
    net = 0
    diff = aligned = 0
    for a in segs:
        rpos = a.reference_start
        for op, n in a.cigartuples:
            if op in (7, 8):
                if op == 8:
                    diff += n
                aligned += n
                rpos += n
            elif op == 2:
                if n >= MIN_CONS_INDEL:
                    ov = min(hi, rpos + n) - max(lo, rpos)
                    if ov > 0:
                        net -= ov
                else:
                    diff += n
                rpos += n
            elif op == 1:
                if n >= MIN_CONS_INDEL:
                    if lo <= rpos <= hi:
                        net += n
                else:
                    diff += n

    def qint(a):
        ct = a.cigartuples
        lead = ct[0][1] if ct[0][0] in (4, 5) else 0
        trail = ct[-1][1] if ct[-1][0] in (4, 5) else 0
        ql = a.query_alignment_length
        return (trail, trail + ql) if a.is_reverse else (lead, lead + ql)

    segs = sorted(segs, key=lambda a: a.reference_start)
    for p, q in zip(segs, segs[1:], strict=False):
        if p.is_reverse != q.is_reverse or not (lo <= p.reference_end <= hi):
            continue
        rgap = q.reference_start - p.reference_end
        pq, qq = qint(p), qint(q)
        qgap = (pq[0] - qq[1]) if p.is_reverse else (qq[0] - pq[1])
        if abs(qgap - rgap) >= MIN_CONS_INDEL:
            net += qgap - rgap
    return net, diff, aligned


def align_container(args):
    crID, chrom, lo, hi, seqs = args
    if not seqs:
        return crID, {}
    ws = max(0, lo - WINDOW)
    with pysam.FastaFile(REFERENCE) as ref:
        wseq = ref.fetch(chrom, ws, hi + WINDOW)
    out = {}
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/w.fa", "w").write(f">w\n{wseq}\n")
        with open(f"{d}/c.fa", "w") as f:
            for cid, s in seqs.items():
                f.write(f">{cid}\n{s}\n")
        subprocess.run(
            f"minimap2 -a -x map-ont --eqx -t 1 {d}/w.fa {d}/c.fa > {d}/o.sam",
            shell=True,
            check=True,
            stderr=subprocess.DEVNULL,
        )
        segs: dict[str, list] = {}
        with pysam.AlignmentFile(f"{d}/o.sam") as sam:
            for a in sam:
                if a.is_unmapped or a.is_secondary:
                    continue
                if a.reference_start < hi - ws and a.reference_end > lo - ws:
                    segs.setdefault(a.query_name, []).append(a)
        for cid in seqs:
            s = segs.get(cid)
            out[cid] = cons_stats(s, lo - ws, hi - ws) if s else None
    return crID, out


def load_rate(d: Path):
    cons: dict[int, dict[str, dict]] = {}
    for f in glob.glob(str(d / "chunk_*.jsonl")):
        for line in open(f):
            cd = json.loads(line)["consensus_dicts"]
            for cid, c in cd.items():
                cons.setdefault(int(cid.split(".")[0]), {})[cid] = c
    runs = pd.read_csv(d / "containers.tsv", sep="\t")
    # lamassemble "using X out of Y sequences" of each container's last attempt
    lam: dict[int, list[tuple[int, int]]] = {}
    for f in glob.glob(str(d / "chunk_*.log")):
        cur = None
        for line in open(f, errors="replace"):
            if line.startswith("CONTAINER "):
                cur = int(line.split()[1])
                lam[cur] = []
            elif cur is not None and line.startswith("lamassemble: using "):
                x = line.split()
                lam[cur].append((int(x[2]), int(x[5])))
    return cons, runs, lam


def evaluate_rate(d: Path, truth: pd.DataFrame, trio, procs: int) -> pd.DataFrame:
    cons, runs, lam = load_rate(d)
    crIDs = sorted(runs.crID.unique())
    jobs = []
    for c in crIDs:
        t = truth.loc[c]
        seqs = {cid: x["consensus_sequence"] for cid, x in cons.get(c, {}).items()}
        jobs.append((c, t.chr, int(t.locus_start), int(t.locus_end), seqs))
    with Pool(procs) as pool:
        aln = dict(pool.imap_unordered(align_container, jobs))
    rows = []
    for c in crIDs:
        t = truth.loc[c]
        r = runs[runs.crID == c]
        last = r.iloc[-1]
        cs = cons.get(c, {})
        # "rep.k", or e.g. "rep.rescue" for the rescue assembly of the reads
        cids = sorted(cs, key=lambda x: (not x.split(".")[1].isdigit(), x))
        metas = [cs[x].get("clustering_meta_data") or {} for x in cids]
        groups = {}
        for k, x in enumerate(cids):
            for iv in cs[x]["intervals_cutread_alignments"]:
                groups.setdefault(iv[2], k)
        st = [aln[c].get(x) for x in cids]
        nets = [s[0] for s in st if s is not None]
        diff = sum(s[1] for s in st if s is not None)
        aligned = sum(s[2] for s in st if s is not None)
        cat = t.T2T_cat
        evaluable = cat in ("het", "hom", "cpx_het", "ref") and t.T2T_bench_frac >= 1
        tn = [int(x) for x in t.T2T_net.split(",")]
        tn = tn[:1] if cat in ("ref", "hom") else tn
        best = 0
        for perm in itertools.permutations(range(len(nets)), min(len(nets), len(tn))):
            best = max(
                best, sum(allele_match(nets[j], tn[i]) for i, j in enumerate(perm))
            )
        distinct = []
        for n in nets:
            if not any(allele_match(n, m) for m in distinct):
                distinct.append(n)
        rows.append(
            {
                "crID": c,
                "trf": bool(t.trf),
                "cat": cat,
                "evaluable": evaluable,
                "n_truth_alleles": len(tn),
                "method": metas[0].get("method", "-") if metas else "none",
                "phasing_status": metas[0].get("phasing_status", "-") if metas else "-",
                "n_consensus": len(cids),
                "n_unaligned": sum(s is None for s in st),
                "cons_net": ",".join("NA" if s is None else str(s[0]) for s in st)
                or "-",
                "recovered": best,
                "all_recovered": best == len(tn),
                "extra": len(distinct) > len(tn),
                "diff": diff,
                "aligned": aligned,
                "pair_acc": pair_accuracy(groups, trio.get(c, {})),
                # lamassemble reports only calls that leave sequences out
                "n_reads": len(groups),
                "lam_dropped": sum(n - u for u, n in lam.get(c, [])),
                "lam_one_read": sum(u == 1 and n > 1 for u, n in lam.get(c, [])),
                "secs": r.secs.sum(),
                "attempts": len(r),
                "unresolved": last.timed_out != "-",
                "timed_out": ";".join(x for x in r.timed_out if x != "-") or "-",
            }
        )
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for rate, g in df.groupby("rate"):
        e = g[g.evaluable]
        row = {
            "rate": rate,
            "containers": len(g),
            "no_cons": int((g.n_consensus == 0).sum()),
            "phased": int((g.method == "phased").sum()),
            "single": int((g.method == "phased_single").sum()),
            "cons/cont": round(g.n_consensus.mean(), 2),
            "pair_acc": round(g.pair_acc.mean(), 4),
            "alleles_rec": round(e.recovered.sum() / e.n_truth_alleles.sum(), 3),
            "cont_rec": round(e.all_recovered.mean(), 3),
            "rec_het": round(e[e.cat != "hom"].all_recovered.mean(), 3),
            "rec_nonTRF": round(e[~e.trf].all_recovered.mean(), 3),
            "rec_TRF": round(e[e.trf].all_recovered.mean(), 3),
            "extra": round(e.extra.mean(), 3),
            "cons_diff%": round(100 * g["diff"].sum() / max(1, g.aligned.sum()), 2),
            "lam_dropped": round(g.lam_dropped.sum() / max(1, g.n_reads.sum()), 3),
            "lam_1read": int((g.lam_one_read > 0).sum()),
            "cpu_s": round(g.secs.sum()),
            "escalated": int((g.attempts > 1).sum()),
            "unresolved": int(g.unresolved.sum()),
            "max_s": round(g.secs.max(), 1),
        }
        out.append(row)
    return pd.DataFrame(out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("run_dir", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--procs", type=int, default=12)
    a = p.parse_args()
    truth = pd.read_csv(DATA / "container_truth.tsv", sep="\t").set_index("crID")
    trio = load_trio()
    parts = []
    for d in sorted(a.run_dir.glob("rate_*")):
        if not (d / "containers.tsv").exists():
            continue
        df = evaluate_rate(d, truth, trio, a.procs)
        df.insert(0, "rate", float(d.name.split("_")[1]))
        parts.append(df)
    df = pd.concat(parts)
    df.to_csv(a.out, sep="\t", index=False)
    s = summarize(df)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 50)
    print(s.to_string(index=False))
    s.to_csv(a.out.with_suffix(".summary.tsv"), sep="\t", index=False)
    e = df[df.evaluable]
    print(
        f"\nevaluable containers (T2TQ100 benchmark, simple categories): {e.crID.nunique()} of {df.crID.nunique()}"
    )
    print("categories:", dict(Counter(e.drop_duplicates("crID").cat)))
    # stability vs no added errors
    base = df[df.rate == 0].set_index("crID")
    if len(base):
        print(
            "\nvs rate 0 (containers): same #consensus / same recovered / lost an allele / gained one"
        )
        for rate, g in df[df.rate > 0].groupby("rate"):
            g = g.set_index("crID")
            b = base.loc[g.index]
            ev = g.evaluable
            print(
                f"  {rate:.2f}: {np.mean(g.n_consensus == b.n_consensus):.3f} / "
                f"{np.mean((g.recovered == b.recovered)[ev]):.3f} / "
                f"{int(((g.recovered < b.recovered) & ev).sum())} / "
                f"{int(((g.recovered > b.recovered) & ev).sum())}"
            )


if __name__ == "__main__":
    main()
