"""Timing and identity check of the reference-site phasing (arm B,
``ref_read_phasing.phase_reads_reference``) on frozen inputs.

  dump  <workdir> <eval.tsv> <out.pkl> [--top 150] [--random 350]
        the inputs of arm B (as ref_vs_ava_eval.py) for the slowest --top
        containers by B_secs and --random others, alignments as SAM text
  run   <in.pkl> <out.json> [--repeat 1]
        arm B on each container: seconds (best of --repeat), status, groups,
        low-quality reads, discordance and a digest of the SNV / SV sites
  diff  <a.json> <b.json>
        containers whose results differ, and the time of both

Run ``run`` once with the old and once with the new code on the PYTHONPATH.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import sys
import time
from pathlib import Path

import pandas as pd
import pysam

from svirlpool.localassembly import read_phasing, ref_read_phasing


def dump(a):
    from kmeans_gate import BUFFER_CLIPPED, MAX_CN, PHASING_FLANK

    from svirlpool.localassembly import consensus
    from svirlpool.localassembly import read_cache as read_cache_mod
    from svirlpool.signalprocessing import copynumber_tracks

    d = pd.read_csv(a.eval, sep="\t")
    d = d[~d.skipped.astype(bool)]
    top = d.nlargest(a.top, "B_secs").crID
    rest = (
        d[~d.crID.isin(top)]
        .sample(n=min(a.random, len(d) - len(top)), random_state=1)
        .crID
    )
    crIDs = sorted(set(top) | set(rest))
    cfg = json.load(open(a.workdir / "config.json"))
    containers = consensus.load_crs_containers_from_db(
        path_db=a.workdir / "crs_containers.db", crIDs=crIDs
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
    out = {"header": None, "items": []}
    try:
        for rep, container in items:
            crs = {cr.crID: cr for cr in container["crs"]}
            chrom = min(cr.chr for cr in crs.values())
            cn = max(
                2,
                copynumber_tracks.query_copynumber_from_regions(
                    bgzip_bed=a.workdir / "copy_number_tracks.bed.gz",
                    regions=[
                        (c.chr, c.referenceStart, c.referenceEnd) for c in crs.values()
                    ],
                ),
            )
            if cn > MAX_CN:
                cache.advance(window_start=float("inf"))
                continue
            alns, recs = {}, {}
            for cr in crs.values():
                al, s = cache.fetch_for_cr(cr)
                alns[cr.crID] = al
                recs.update(s)
            iv = consensus.get_read_alignment_intervals_in_cr(
                crs=list(crs.values()),
                dict_alignments=alns,
                buffer_clipped_length=BUFFER_CLIPPED,
            )
            cutreads = consensus.trim_reads(
                dict_alignments=alns,
                intervals=consensus.get_max_extents_of_read_alignments_on_cr(iv),
                read_records=recs,
            )
            iv = consensus.get_read_alignment_intervals_in_cr(
                crs=list(crs.values()),
                dict_alignments=alns,
                buffer_clipped_length=BUFFER_CLIPPED,
                flank=PHASING_FLANK,
            )
            preads = consensus.trim_reads(
                dict_alignments=alns,
                intervals=consensus.get_max_extents_of_read_alignments_on_cr(iv),
                read_records=recs,
            )
            preads = {rn: r for rn, r in preads.items() if rn in cutreads}
            flat = [x for al in alns.values() for x in al]
            if out["header"] is None and flat:
                out["header"] = flat[0].header.to_dict()
            on_chr = [c for c in crs.values() if c.chr == chrom]
            out["items"].append(
                {
                    "crID": rep,
                    "reads": preads,
                    "alns": [x.to_string() for x in flat],
                    "chrom": chrom,
                    "start": min(c.referenceStart for c in on_chr),
                    "end": max(c.referenceEnd for c in on_chr),
                    "flank": PHASING_FLANK,
                }
            )
            cache.advance(window_start=float("inf"))
    finally:
        cache.close()
    out["reference"] = cfg["reference"]
    pickle.dump(out, open(a.out, "wb"))
    print(f"{len(out['items'])} containers -> {a.out}")


def _digest(sites):
    h = hashlib.sha1()
    for s in sorted(sites, key=lambda s: (s.pos, s.kind, s.t)):
        h.update(
            repr(
                (s.t, s.pos, s.kind, sorted(s.agree), sorted(s.disagree), s.depth)
            ).encode()
        )
    return h.hexdigest()[:12]


def run(a):
    data = pickle.load(open(a.inp, "rb"))
    header = pysam.AlignmentHeader.from_dict(data["header"])
    ref = pysam.FastaFile(data["reference"])
    p = read_phasing.PhasingParams()
    rows = {}
    # the site lists, through the same steps as phase_reads_reference
    orig = read_phasing.phase_from_sites
    seen = {}

    def spy(all_reads, reads, s_snv, s_sv, lowq, params):
        seen["snv"], seen["sv"] = _digest(s_snv), _digest(s_sv)
        return orig(all_reads, reads, s_snv, s_sv, lowq, params)

    ref_read_phasing.phase_from_sites = spy
    for it in data["items"]:
        alns = [pysam.AlignedSegment.fromstring(s, header) for s in it["alns"]]
        best = float("inf")
        for _ in range(a.repeat):
            t0 = time.perf_counter()
            res = ref_read_phasing.phase_reads_reference(
                reads=it["reads"],
                alns=alns,
                ref=ref,
                chrom=it["chrom"],
                start=it["start"],
                end=it["end"],
                flank=it["flank"],
                params=p,
            )
            best = min(best, time.perf_counter() - t0)
        rows[str(it["crID"])] = {
            "secs": best,
            "status": res.status,
            "n_alleles": res.n_alleles,
            "groups": res.groups,
            "lowq": sorted(res.low_quality),
            "discordance": res.discordance,
            "n_snv": res.n_snv_sites,
            "n_sv": res.n_sv_sites,
            "snv": seen.get("snv"),
            "sv": seen.get("sv"),
        }
        seen.clear()
    json.dump(rows, open(a.out, "w"))
    s = sum(r["secs"] for r in rows.values())
    print(f"{len(rows)} containers, {s:.1f} s")


def diff(a):
    x, y = json.load(open(a.a)), json.load(open(a.b))
    keys = sorted(set(x) & set(y), key=int)
    nd = 0
    for k in keys:
        rx, ry = dict(x[k]), dict(y[k])
        tx, ty = rx.pop("secs"), ry.pop("secs")
        dx, dy = rx.pop("discordance"), ry.pop("discordance")
        same_d = (dx is None and dy is None) or (
            dx is not None and dy is not None and abs(dx - dy) < 1e-12
        )
        if rx != ry or not same_d:
            nd += 1
            fields = [f for f in rx if rx[f] != ry[f]] + (
                [] if same_d else ["discordance"]
            )
            print(f"crID {k}: differs in {fields}  ({tx:.2f} s vs {ty:.2f} s)")
    ta = sum(x[k]["secs"] for k in keys)
    tb = sum(y[k]["secs"] for k in keys)
    mx = max(x[k]["secs"] for k in keys), max(y[k]["secs"] for k in keys)
    print(
        f"{len(keys)} containers, {nd} differ; total {ta:.1f} s -> {tb:.1f} s ({ta / tb:.1f}x); "
        f"max {mx[0]:.2f} s -> {mx[1]:.2f} s"
    )
    return nd


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("dump")
    d.add_argument("workdir", type=Path)
    d.add_argument("eval", type=Path)
    d.add_argument("out", type=Path)
    d.add_argument("--top", type=int, default=150)
    d.add_argument("--random", type=int, default=350)
    r = sub.add_parser("run")
    r.add_argument("inp", type=Path)
    r.add_argument("out", type=Path)
    r.add_argument("--repeat", type=int, default=1)
    c = sub.add_parser("diff")
    c.add_argument("a", type=Path)
    c.add_argument("b", type=Path)
    a = p.parse_args()
    if a.cmd == "dump":
        dump(a)
    elif a.cmd == "run":
        run(a)
    else:
        sys.exit(1 if diff(a) else 0)


if __name__ == "__main__":
    main()
