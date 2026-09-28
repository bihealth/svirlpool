"""Phasing alone on dumped read sets with added errors, per parameter variant.

Reads: ``consensus_perf/dump_phasing_reads.py --orient`` output (<crID>.fa,
the reads as ``phase_reads`` gets them). Errors are added to every read
(``consensus.add_read_errors``). Scores the het containers (T2TQ100 het /
cpx_het): phased (>= 2 alleles) and trio pair accuracy; the hom / ref ones:
wrongly phased.

usage: phase_diag.py <reads_dir> <crIDs.txt> [--rates 0.04,0.06,0.08]
                     [--variants name:key=val;key=val ...] [--procs 12]
"""

from __future__ import annotations

import argparse
import functools
import logging
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from svirlpool.localassembly import read_phasing
from svirlpool.localassembly.consensus import add_read_errors

logging.disable(logging.WARNING)
HERE = Path(__file__).parent
DATA = HERE.parent / "ava_phasing" / "data"

VARIANTS = {
    "base": {},
    "conc0.8": {"min_conc": 0.8},
    "no_indel_mask": {"indel_window": 0},
    "min_alt5": {"min_alt": 5},
    "recur2": {"recurrence": 2},
}


def run_one(args):
    path, rate, variant, kw = args
    kw = dict(kw)
    min_conc = kw.pop("min_conc", None)
    if min_conc is not None:
        read_phasing.recurrent_sites = functools.partial(_RECURRENT, min_conc=min_conc)
    else:
        read_phasing.recurrent_sites = _RECURRENT
    reads = {r.id: r for r in SeqIO.parse(path, "fasta")}
    if rate > 0:
        reads = {n: add_read_errors(r, start=0, rate=rate) for n, r in reads.items()}
    res = read_phasing.phase_reads(
        reads, params=read_phasing.PhasingParams(**kw), timeout=60
    )
    return {
        "crID": int(path.stem),
        "rate": rate,
        "variant": variant,
        "status": res.status,
        "n_alleles": res.n_alleles,
        "n_snv": res.n_snv_sites,
        "groups": res.groups,
    }


_RECURRENT = read_phasing.recurrent_sites


def main():
    import eval_rates  # noqa: PLC0415  (same directory)

    p = argparse.ArgumentParser()
    p.add_argument("reads_dir", type=Path)
    p.add_argument("crIDs", type=Path)
    p.add_argument("--rates", default="0.04,0.06,0.08")
    p.add_argument("--variants", default=",".join(VARIANTS))
    p.add_argument("--procs", type=int, default=12)
    a = p.parse_args()
    ids = [int(x) for x in a.crIDs.read_text().split()]
    rates = [float(x) for x in a.rates.split(",")]
    jobs = [
        (a.reads_dir / f"{c}.fa", r, v, VARIANTS[v])
        for r in rates
        for v in a.variants.split(",")
        for c in ids
        if (a.reads_dir / f"{c}.fa").exists()
    ]
    with Pool(a.procs) as pool:
        rows = list(pool.imap_unordered(run_one, jobs, chunksize=4))
    truth = pd.read_csv(DATA / "container_truth.tsv", sep="\t").set_index("crID")
    trio = eval_rates.load_trio()
    for r in rows:
        r["cat"] = truth.loc[r["crID"], "T2T_cat"]
        r["acc"] = eval_rates.pair_accuracy(r.pop("groups"), trio.get(r["crID"], {}))
    df = pd.DataFrame(rows)
    out = []
    for (rate, v), g in df.groupby(["rate", "variant"]):
        het = g[g.cat.isin(["het", "cpx_het"])]
        hom = g[g.cat.isin(["hom", "ref"])]
        out.append(
            {
                "rate": rate,
                "variant": v,
                "het_phased": round((het.status == "phased").mean(), 3),
                "het_acc": round(np.nanmean(het.acc), 4),
                "hom_phased": round((hom.status == "phased").mean(), 3),
                "no_info": int((g.status == "no_information").sum()),
                "med_snv": g.n_snv.median(),
            }
        )
    pd.set_option("display.width", 200)
    print(pd.DataFrame(out).to_string(index=False))


if __name__ == "__main__":
    main()
