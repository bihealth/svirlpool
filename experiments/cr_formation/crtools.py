"""Candidate-region (CR) formation study: loaders and a scorer for any CR set.

Inputs from the svp_tiered15 ps_tiered HG002 run (15% regions). Truth (Q100)
is used only to SCORE CR sets, never to form them.
"""

from __future__ import annotations

import gzip
import json
import pickle
import sqlite3
import subprocess

import numpy as np
import pandas as pd

S = "/home/mayv_c/development/svp_tiered15/"
W = S + "results/ps_tiered/20x/HG002/work/"
TV = S + "results/ps_tiered/20x/truvari/T2TQ100/all/"
TRF = "/home/mayv_c/biodata/local/references/GRCh38/GRCh38.trf.bed"
CUT_FLANK = 200


def load_signals(path=W + "signalstrength.tsv.gz") -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as f:
        for line in f:
            d = json.loads(line.split("\t", 3)[3].strip().strip('"').replace('""', '"'))
            rows.append((d["chr"], d["ref_start"], d["ref_end"], d["size"], d["sv_type"],
                         d["readname"], d["coverage"], d["repeatID"], d["strength"]))
    return pd.DataFrame(rows, columns=["chr", "start", "end", "size", "type", "read",
                                       "cov", "rep", "strength"])


def load_crs(path=W + "crs.db") -> pd.DataFrame:
    rows = []
    for crID, blob in sqlite3.connect(path).execute("select crID, candidate_region from candidate_regions"):
        d = pickle.loads(blob)
        sig = d["sv_signals"]
        rows.append({"crID": crID, "chr": d["chr"], "start": d["referenceStart"], "end": d["referenceEnd"],
                         "n_sig": len(sig), "n_reads": len({s["readname"] for s in sig}),
                         "cov": float(np.median([s["coverage"] for s in sig])) if sig else 0.0,
                         "signals": sig})
    return pd.DataFrame(rows)


def load_truth() -> pd.DataFrame:
    """The truvari base set of T2TQ100 'all' (tp-base + fn), with TP status."""
    rows = []
    for fn, tp in (("tp-base.vcf.gz", True), ("fn.vcf.gz", False)):
        out = subprocess.run(["bash", "-c", f"zcat {TV}{fn} | grep -v '^#' | cut -f1,2,4,5"],
                             capture_output=True, text=True).stdout
        for line in out.splitlines():
            c, p, r, a = line.split("\t")
            p = int(p)
            rows.append((c, p, p + max(1, len(r) - 1) if len(r) > len(a) else p + 1,
                         abs(len(r) - len(a)), tp))
    return pd.DataFrame(rows, columns=["chr", "start", "end", "size", "tp"])


def covered(truth: pd.DataFrame, crs: pd.DataFrame, flank: int = CUT_FLANK) -> np.ndarray:
    """Whether each truth SV overlaps a CR (+- flank): a necessary condition to call it."""
    hit = np.zeros(len(truth), dtype=bool)
    for c, g in crs.groupby("chr"):
        s = np.sort(g.start.values - flank)
        order = np.argsort(g.start.values)
        e = (g.end.values + flank)[order]
        emax = np.maximum.accumulate(e)
        m = truth.chr.values == c
        ts, te = truth.start.values[m], truth.end.values[m]
        k = np.searchsorted(s, te, side="right")  # CRs starting at or before the SV end
        ok = (k > 0) & (emax[np.maximum(k - 1, 0)] >= ts)
        hit[np.where(m)[0]] = ok
    return hit


def cost(crs: pd.DataFrame) -> pd.Series:
    """Proxy of the bp of cut reads to assemble: depth x (span + 2 flanks)."""
    return crs["cov"] * (crs.end - crs.start + 2 * CUT_FLANK)


def evaluate(name: str, crs: pd.DataFrame, truth: pd.DataFrame, heavy: float = 100_000) -> dict:
    hit = covered(truth, crs)
    c = cost(crs)
    span = crs.end - crs.start
    return {
        "variant": name, "n_crs": len(crs), "span_med": float(span.median()),
        "span_p99": float(span.quantile(.99)), "n_span_gt5k": int((span > 5000).sum()),
        "cost_Mbp": float(c.sum() / 1e6), "n_heavy": int((c >= heavy).sum()),
        "heavy_cost_Mbp": float(c[c >= heavy].sum() / 1e6),
        "tp_covered": float(hit[truth.tp.values].mean()),
        "tp_lost": int((~hit[truth.tp.values]).sum()),
        "fn_covered": float(hit[~truth.tp.values].mean()),
    }
