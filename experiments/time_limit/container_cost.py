"""Where does a consensus container's time go, and what predicts it?

Per container (representative crID) of one svirlpool run, from its consensus
batch logs, the wall seconds split into stages (all escalation levels summed):

  prep       copy number, read fetch, cutting the reads
  ref_phase  reference-SNV phasing (--phasing-sites tiered/reference)
  ava        all-vs-all phasing: orienting the reads, minimap2, sites, clustering
  post_phase phasing done -> first lamassemble (writing the cluster reads)
  lam        lamassemble (LAST all-vs-all, layout, anchors, MAFFT, consensus)
  align      minimap2 of the reads to their consensus, re-adding unused reads

and the features the run logs before or while assembling (reads, cut bp,
alleles, AVA reads, ...). Joined with features of the container itself
(crs_containers.db: CR spans, signals, repeat IDs, depth), of the reference
(TRF overlap, k-mer self-repetitiveness of the CR sequence) and, for the
sample with truth, its value (container_value.py: Q100 representation,
truvari TP / FP calls).

usage: container_cost.py <work dir> <out.tsv> [--sample-results DIR] [--reference FA] [--trf BED]
  <work dir>          .../<sample>/work of a svirlpool run (consensus/*/consensus.batch_*.log)
  --sample-results    .../<variant>/20x of the tuning results, for the truth of HG002
"""

from __future__ import annotations

import argparse
import bisect
import glob
import json
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import container_value as CV  # noqa: E402

TS = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}) - (\S+) - (\w+) - (.*)$")
START = re.compile(r"PROGRESS \[\d+/\d+\] Processing container \(representative crID (\d+)\)")
FETCH = re.compile(r"fetch_for_cr\((\S+)\) found (\d+) alns")
PHASE = re.compile(
    r"read phasing: status=(\w+) alleles=(\d+) sizes=\[([\d, ]*)\] unassigned=(\d+) "
    r"low_quality=(\d+) snv_sites=(\d+) sv_sites=(\d+)"
)
AVA_N = re.compile(r"--secondary=yes -N (\d+)")
LAM = re.compile(r"Running lamassemble on \S+ for \S+ \(both strands: \w+\) with timeout of (\d+) seconds")
LAM_SIZE = re.compile(r"\((\d+) reads, (\d+) bp\)$")
TIME_DROP = re.compile(r"Container \d+: dropped at the container time limit")
SIZE_DROP = re.compile(r"Container \d+: dropped at --max-assembly-bp: (\d+) bp")
M_RETRY = re.compile(r"lamassemble: \d+ of \d+ reads linked at -m (\d+)")
ESC = re.compile(r"Container \d+: (.*) timed out \((\d+) thread")
LAST = re.compile(r"Container \d+: (.*) timed out at the last level")
LAM_TO = re.compile(r"lamassemble timed out for")


def _t(s: str) -> float:
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S,%f").timestamp()


def parse_batch(path: str) -> list[dict]:
    """One record per container of one batch log."""
    out: list[dict] = []
    cur: dict | None = None
    state, t_prev = "prep", 0.0
    lam_start: float | None = None
    lam_limit = 0
    phase_lines = 0  # read phasing result lines in the current attempt
    attempt = 0
    level_threads = [1]  # threads of each escalation level
    level_start = 0.0

    def close(t_end: float):
        nonlocal cur
        if cur is not None:
            cur["sec_" + state] += t_end - t_prev
            cur["seconds"] = t_end - cur["t0"]
            cur[f"sec_level{attempt}"] += t_end - level_start
            # thread-seconds: wall of each level x its threads, an upper bound
            # of the CPU (a level waiting on I/O uses less)
            cur["thread_seconds"] = sum(
                cur.get(f"sec_level{i}", 0.0) * level_threads[min(i, len(level_threads) - 1)]
                for i in range(attempt + 1)
            )
            out.append(cur)
            cur = None

    for line in open(path):
        m = TS.match(line)
        if not m:
            continue
        t, src, level, msg = _t(m.group(1)), m.group(2), m.group(3), m.group(4)
        if msg.startswith("escalation levels (threads, timeout s):"):
            level_threads = [int(x) for x in re.findall(r"\((\d+), \d+\)", msg)]
        if (s := START.search(msg)) or (msg.startswith("Wrote ") and "container results" in msg):
            close(t)
            if s:
                cur = defaultdict(float)
                cur.update(crID=int(s.group(1)), t0=t, batch=path, attempts=1)
                state, t_prev, attempt, phase_lines = "prep", t, 0, 0
                level_start = t
            continue
        if cur is None:
            continue
        cur["sec_" + state] += t - t_prev
        t_prev = t
        if lam_start is not None and state == "lam":
            # the lamassemble call that started at lam_start ended at t
            took = t - lam_start
            cur["lam_calls"] += 1
            cur["lam_max_s"] = max(cur["lam_max_s"], took)
            # below the last level the escalation line is the only trace
            if LAM_TO.search(msg) or ("lamassemble timed out" in msg and level == "WARNING"):
                cur["lam_timeouts"] += 1
                cur["lam_overshoot_s"] = max(cur["lam_overshoot_s"], took - lam_limit)
            lam_start = None
        if f := FETCH.search(msg):
            if attempt == 0:
                cur["n_crs_fetched"] += 1
                cur["alns"] += int(f.group(2))
            state = "prep"
        elif msg.startswith("number of reads:") and attempt == 0:
            cur["n_reads"] = int(msg.split(":")[1])
        elif msg.startswith("summed trimmed reads bp:"):
            if attempt == 0:
                cur["cut_bp"] = int(msg.split(":")[1])
            state = "ref_phase"
        elif p := PHASE.search(msg):
            phase_lines += 1
            sizes = [int(x) for x in p.group(3).split(",") if x.strip()]
            # with tiered sites the first result is the reference one
            tag = "ref" if (phase_lines == 1 and src.endswith("read_phasing") and state == "ref_phase") else "ava"
            if attempt == 0 or f"{tag}_status" not in cur:
                cur[f"{tag}_status"] = p.group(1)
                cur[f"{tag}_alleles"] = int(p.group(2))
                cur[f"{tag}_maxsize"] = max(sizes, default=0)
                cur[f"{tag}_unassigned"] = int(p.group(4))
                cur[f"{tag}_lowq"] = int(p.group(5))
                cur[f"{tag}_snv"] = int(p.group(6))
                cur[f"{tag}_sv"] = int(p.group(7))
            state = "ava" if tag == "ref" else "post_phase"
        elif src.endswith("read_phasing") and (a := AVA_N.search(msg)):
            cur["ava_reads"] = int(a.group(1))
            state = "ava"
        elif msg.startswith("read phasing: reused"):
            state = "post_phase"
        elif lm := LAM.search(msg):
            lam_start, lam_limit = t, int(lm.group(1))
            state = "lam"
            if sz := LAM_SIZE.search(msg):
                bp = int(sz.group(2))
                cur["max_call_bp"] = max(cur["max_call_bp"], bp)
                cur["max_call_reads"] = max(cur["max_call_reads"], int(sz.group(1)))
        elif M_RETRY.search(msg):
            cur["m_retries"] += 1
        elif TIME_DROP.search(msg):
            cur["time_dropped"] = 1
        elif sd := SIZE_DROP.search(msg):
            cur["size_dropped"] = 1
            cur["size_dropped_bp"] = int(sd.group(1))
        elif src.endswith("util.util") and msg.startswith("minimap2"):
            state = "align"
        elif level == "WARNING" and (e := ESC.search(msg)):
            cur["escalations"] += 1
            cur["esc_tools"] = ";".join(sorted(set(filter(None, [cur.get("esc_tools") or "", e.group(1)]))))
            cur["attempts"] += 1
            cur[f"sec_level{attempt}"] += t - level_start
            level_start = t
            attempt += 1
            phase_lines = 0
            state = "prep"
        elif level == "WARNING" and LAST.search(msg):
            cur["last_level_timeout"] = 1
        elif "read phasing: all-vs-all alignment failed" in msg:
            cur["ava_failed"] = 1
            state = "post_phase"
    return out


def log_table(work: str) -> pd.DataFrame:
    rows = []
    for path in sorted(glob.glob(work + "/consensus/*/consensus.batch_*.log")):
        if not path.endswith(".diag.log"):
            rows.extend(parse_batch(path))
    d = pd.DataFrame(rows).set_index("crID")
    d["batch"] = d.batch.str.extract(r"batch_(\d+)\.log").astype(int)
    for c in d.columns:
        if c.startswith("sec_") or c in ("lam_calls", "lam_timeouts", "escalations", "alns",
                                         "n_crs_fetched", "lam_max_s", "lam_overshoot_s"):
            d[c] = d[c].fillna(0)
    d["ava_ran"] = d.ava_reads.notna()
    return d.drop(columns="t0")


def db_table(db: str) -> pd.DataFrame:
    rows = {}
    for cid, data in sqlite3.connect(db).execute("select crID, data from containers"):
        c = json.loads(data)
        crs = c["crs"]
        sig = [s for cr in crs for s in cr["sv_signals"]]
        meds = [float(np.median([s["coverage"] for s in cr["sv_signals"]])) for cr in crs if cr["sv_signals"]]
        sizes = [abs(s["size"]) for s in sig]
        rows[cid] = {
            "n_crs": len(crs),
            "n_chr": len({cr["chr"] for cr in crs}),
            "chr": crs[0]["chr"],
            "start": min(cr["referenceStart"] for cr in crs),
            "span_sum": sum(cr["referenceEnd"] - cr["referenceStart"] for cr in crs),
            "span_max": max(cr["referenceEnd"] - cr["referenceStart"] for cr in crs),
            "n_signals": len(sig),
            "n_signal_reads": len({s["readname"] for s in sig}),
            "signals_per_read": len(sig) / max(1, len({s["readname"] for s in sig})),
            "repeat_share": np.mean([s.get("repeatID", -1) not in (-1, None) for s in sig]) if sig else 0.0,
            "n_repeat_ids": len({s["repeatID"] for s in sig if s.get("repeatID", -1) not in (-1, None)}),
            "depth": max(meds, default=0.0),
            "sig_size_sum": float(sum(sizes)),
            "sig_size_max": float(max(sizes, default=0)),
            "n_connecting": len(c.get("connecting_reads") or []),
            "crs": [(cr["chr"], cr["referenceStart"], cr["referenceEnd"]) for cr in crs],
        }
    return pd.DataFrame.from_dict(rows, orient="index")


class Trf:
    def __init__(self, bed: str):
        b = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "s", "e"])
        self.by = {c: (g.s.to_numpy(), g.e.to_numpy()) for c, g in b.sort_values(["chr", "s"]).groupby("chr")}

    def overlap(self, chrom: str, s: int, e: int, pad: int = 0) -> tuple[int, int]:
        """(bp of [s, e) covered by TRF, length of the longest TRF interval
        overlapping [s - pad, e + pad))."""
        if chrom not in self.by:
            return 0, 0
        ss, ee = self.by[chrom]
        i = max(0, bisect.bisect_left(ss, s - pad - 200_000))
        cov = np.zeros(max(1, e - s), dtype=bool)
        longest = 0
        while i < len(ss) and ss[i] < e + pad:
            if ee[i] > s - pad:
                longest = max(longest, int(ee[i] - ss[i]))
                a, b = max(s, ss[i]) - s, min(e, ee[i]) - s
                if b > a:
                    cov[a:b] = True
            i += 1
        return int(cov.sum()), longest


def kmer_repetitiveness(seq: str, k: int = 15) -> tuple[float, int]:
    """Share of the k-mer positions of `seq` whose k-mer occurs again in it, and
    the largest count of one k-mer: how many offsets a read of this sequence
    aligns to itself at (tandem copies)."""
    seq = seq.upper()
    n = len(seq) - k + 1
    if n <= 0:
        return 0.0, 0
    cnt = Counter(seq[i : i + k] for i in range(n))
    rep = sum(v for v in cnt.values() if v > 1)
    return rep / n, max(cnt.values())


def ref_table(db: pd.DataFrame, reference: str, trf_bed: str, flank: int = 1000) -> pd.DataFrame:
    import pysam

    fa = pysam.FastaFile(reference)
    trf = Trf(trf_bed)
    rows = {}
    for cid, crs in db.crs.items():
        span = cov = longest = 0
        rep_w = maxk = window = 0
        for chrom, s, e in crs:
            c, lg = trf.overlap(chrom, s, e, pad=flank)
            span += e - s
            cov += c
            longest = max(longest, lg)
            lo, hi = max(0, s - flank), min(fa.get_reference_length(chrom), e + flank)
            if hi - lo > 200_000:
                hi = lo + 200_000
            r, mk = kmer_repetitiveness(fa.fetch(chrom, lo, hi))
            rep_w += r * (hi - lo)
            window += hi - lo
            maxk = max(maxk, mk)
        rows[cid] = {
            "trf_frac": cov / max(1, span),
            "trf_longest": longest,
            "ref_kmer_rep": rep_w / max(1, window),
            "ref_kmer_max": maxk,
        }
    return pd.DataFrame.from_dict(rows, orient="index")


def truth_table(results: str, sample: str = "HG002", truth: str = "T2TQ100") -> pd.DataFrame:
    t = pd.DataFrame({"repr": CV.representation(f"{results}/{sample}/consensus_q100.tsv")})
    tv = f"{results}/truvari/{truth}/all"
    t = t.join(pd.DataFrame({
        "tp": CV.per_container(CV.call_containers(f"{tv}/tp-comp.vcf.gz")),
        "fp": CV.per_container(CV.call_containers(f"{tv}/fp.vcf.gz")),
    }), how="outer")
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("out")
    ap.add_argument("--sample-results", default="")
    ap.add_argument("--reference", default="")
    ap.add_argument("--trf", default="")
    a = ap.parse_args()
    d = log_table(a.work)
    db = db_table(a.work + "/crs_containers.db")
    d = d.join(db, how="left")
    if a.reference and a.trf:
        d = d.join(ref_table(db.loc[db.index.intersection(d.index)], a.reference, a.trf), how="left")
    if a.sample_results:
        t = truth_table(a.sample_results)
        d = d.join(t, how="left")
        d[["tp", "fp"]] = d[["tp", "fp"]].fillna(0)
    d = d.drop(columns="crs")
    d.index.name = "crID"
    d.to_csv(a.out, sep="\t")
    print(f"{len(d)} containers, {d.seconds.sum() / 3600:.1f} h; written to {a.out}")


if __name__ == "__main__":
    main()
