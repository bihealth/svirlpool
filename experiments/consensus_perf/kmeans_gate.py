"""Could the KMeans gate be a fast path in front of the read phasing?

In phased mode the KMeans clustering of the summed indels only runs when the
phasing path produces no consensus. Its acceptance gate
(``consensus.kmeans_partition``) is cheap; the phasing's all-vs-all is not. Per
container this measures, on the reads the consensus stage uses:

* the gate: accepted k (0 = rejected) and its time, plus the time of the
  summed-indel features it needs (computed in every mode anyway);
* the phasing: status, allele count, time (1 thread, as the first
  escalation level; --timeout generous so the time is the real cost);
* both partitions at the SV-allele level against trio read labels, with the
  T2TQ100 category deciding which haplotypes carry the same allele (hom / ref:
  one allele), as ``../ava_phasing/allele_eval.py``: purity, all alleles
  recovered, fraction of labelled reads assigned. The phasing partition is the
  one production uses: its clusters when phased, else one cluster of all
  reads but the low-quality ones (--phasing-fallback single).

usage: kmeans_gate.py <workdir> <out.tsv> [--procs 24] [--timeout 300] [--limit N]
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from phasing_eval import load_trio

from svirlpool.localassembly import consensus, read_phasing
from svirlpool.localassembly import read_cache as read_cache_mod
from svirlpool.signalprocessing import copynumber_tracks

logging.disable(logging.WARNING)
HERE = Path(__file__).parent
TRUTH = HERE.parent / "ava_phasing" / "data" / "container_truth.tsv"
BUFFER_CLIPPED = 500  # consensus --buffer-clipped-sequence default
PHASING_FLANK = 10000
MAX_CN = 4  # containers above are skipped by the consensus stage
# as in process_consensus_container
VARIANCE_THRESHOLD = 29.0
DISTANCE_THRESHOLD = 29.0


def one_batch(args):
    wd, crIDs, cfg, timeout = args
    containers = consensus.load_crs_containers_from_db(
        path_db=wd / "crs_containers.db", crIDs=crIDs
    )
    items = sorted(
        containers.items(),
        key=lambda it: (
            min(cr.chr for cr in it[1]["crs"]),
            min(cr.referenceStart for cr in it[1]["crs"]),
            it[0],
        ),
    )
    cache = read_cache_mod.ReadSequenceCache(path_alignments=Path(cfg["alignments"]))
    rows = []
    try:
        for rep, container in items:
            crs = {cr.crID: cr for cr in container["crs"]}
            row = {
                "crID": rep,
                "chr": min(cr.chr for cr in crs.values()),
                "start": min(cr.referenceStart for cr in crs.values()),
            }
            cn = max(
                2,
                copynumber_tracks.query_copynumber_from_regions(
                    bgzip_bed=wd / "copy_number_tracks.bed.gz",
                    regions=[(c.chr, c.referenceStart, c.referenceEnd) for c in crs.values()],
                ),
            )
            row["cn"] = cn
            if cn > MAX_CN:
                row["skipped"] = True
                rows.append(row)
                cache.advance(window_start=float("inf"))
                continue
            row["skipped"] = False
            alns, recs = {}, {}
            for cr in crs.values():
                a, s = cache.fetch_for_cr(cr)
                alns[cr.crID] = a
                recs.update(s)
            iv = consensus.get_read_alignment_intervals_in_cr(
                crs=list(crs.values()), dict_alignments=alns,
                buffer_clipped_length=BUFFER_CLIPPED,
            )
            cutreads = consensus.trim_reads(
                dict_alignments=alns,
                intervals=consensus.get_max_extents_of_read_alignments_on_cr(iv),
                read_records=recs,
            )
            row["n_reads"] = len(cutreads)

            signals = [s for cr in crs.values() for s in cr.sv_signals]
            row["sig_trf_frac"] = (
                sum(s.repeatID != -1 for s in signals) / len(signals) if signals else 0.0
            )
            t0 = time.perf_counter()
            indels = consensus.summed_indel_distribution(alns=alns, crs=crs)
            row["t_indels"] = time.perf_counter() - t0
            row["indels"] = json.dumps(
                {rn: v for rn, v in indels.items() if rn in cutreads}, sort_keys=True
            )
            t0 = time.perf_counter()
            part = consensus.kmeans_partition(
                dict_summed_indels=indels, pool=cutreads, max_k=cn,
                variance_threshold=VARIANCE_THRESHOLD,
                distance_threshold=DISTANCE_THRESHOLD,
            )
            row["t_gate"] = time.perf_counter() - t0
            if part is None:
                row["gate_k"] = 0
                row["kmeans_groups"] = "{}"
            else:
                names, labels, k = part
                row["gate_k"] = k
                row["kmeans_groups"] = json.dumps(
                    {rn: int(lb) for rn, lb in zip(names, labels, strict=True)
                     if rn in cutreads}, sort_keys=True)

            # the phasing exactly as consensus_while_phasing prepares it
            t0 = time.perf_counter()
            iv = consensus.get_read_alignment_intervals_in_cr(
                crs=list(crs.values()), dict_alignments=alns,
                buffer_clipped_length=BUFFER_CLIPPED, flank=PHASING_FLANK,
            )
            preads = consensus.trim_reads(
                dict_alignments=alns,
                intervals=consensus.get_max_extents_of_read_alignments_on_cr(iv),
                read_records=recs,
            )
            preads = consensus.orient_reads_to_reference(
                {rn: r for rn, r in preads.items() if rn in cutreads}, alns
            )
            res = read_phasing.phase_reads(reads=preads, threads=1, timeout=timeout)
            row["t_phase"] = time.perf_counter() - t0
            row["phase_status"] = res.status
            row["phase_n_alleles"] = res.n_alleles
            if res.status == "phased":
                groups = {rn: cid for cid, rns in res.clusters().items() for rn in rns}
            else:
                lowq = set(res.low_quality)
                groups = {rn: 0 for rn in cutreads if rn not in lowq}
            row["phase_groups"] = json.dumps(groups, sort_keys=True)
            rows.append(row)
            cache.advance(window_start=float("inf"))
    finally:
        cache.close()
    return rows


def allele_eval(labels: dict[str, int], trio: dict[str, str], cat: str, min_cl=3):
    """As ../ava_phasing/allele_eval.py: pat = A1, mat = A1 for hom/ref else A2."""
    amap = {"pat": "A1", "mat": "A1" if cat in ("hom", "ref") else "A2"}
    rows = [(labels[r], amap[trio[r]]) for r in labels if trio.get(r) in amap]
    n_lab = sum(1 for r in trio if trio[r] in amap)
    if len(rows) < 4:
        return None
    by_cl = defaultdict(Counter)
    for cl, al in rows:
        by_cl[cl][al] += 1
    purity = sum(c.most_common(1)[0][1] for c in by_cl.values()) / len(rows)
    sizes = Counter(cl for cl, _ in rows)
    maj = {cl: c.most_common(1)[0][0] for cl, c in by_cl.items() if sizes[cl] >= min_cl}
    recovered = {al for _, al in rows} <= set(maj.values())
    return {"purity": purity, "recovered": recovered, "assigned": len(rows) / max(n_lab, 1)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("workdir", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--procs", type=int, default=24)
    p.add_argument("--timeout", type=int, default=300)
    p.add_argument("--limit", type=int, default=0, help="first N batches only")
    a = p.parse_args()
    cfg = json.load(open(a.workdir / "config.json"))
    jobs = []
    with open(a.workdir / "consensus_batches.tsv") as f:
        ic = f.readline().rstrip("\n").split("\t").index("crIDs")
        for line in f:
            crIDs = [int(x) for x in line.rstrip("\n").split("\t")[ic].split(",")]
            jobs.append((a.workdir, crIDs, cfg, a.timeout))
    if a.limit:
        jobs = jobs[: a.limit]
    # finer jobs for load balance: one container list per job is a batch; split
    jobs = [
        (wd, crIDs[i : i + 10], c, t)
        for wd, crIDs, c, t in jobs
        for i in range(0, len(crIDs), 10)
    ]
    t0 = time.perf_counter()
    with Pool(a.procs) as pool:
        rows = [r for rs in pool.imap_unordered(one_batch, jobs) for r in rs]
    print(f"wall {time.perf_counter() - t0:.0f} s, {len(rows)} containers")

    d = pd.DataFrame(rows).sort_values("crID")
    truth = pd.read_csv(TRUTH, sep="\t").set_index("crID")
    trio = load_trio()
    ev = []
    for r in d.itertuples():
        out = {"crID": r.crID}
        if r.crID in truth.index:
            t = truth.loc[r.crID]
            out["truth_locus_ok"] = t.chr == r.chr and abs(int(t.start) - int(r.start)) < 5000
            out["cat"] = t.T2T_cat
            out["bench"] = t.T2T_bench_frac >= 1
            out["trf"] = bool(t.trf)
            if not r.skipped:
                for name, col in (("km", "kmeans_groups"), ("ph", "phase_groups")):
                    e = allele_eval(json.loads(getattr(r, col)), trio.get(r.crID, {}), t.T2T_cat)
                    for k, v in (e or {}).items():
                        out[f"{name}_{k}"] = v
        ev.append(out)
    d = d.merge(pd.DataFrame(ev), on="crID", how="left")
    d.to_csv(a.out, sep="\t", index=False)
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
