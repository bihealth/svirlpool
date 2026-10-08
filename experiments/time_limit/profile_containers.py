"""Profile the consensus of single containers step by step.

Runs `consensus.crs_containers_to_consensus` on the given containers of a
svirlpool work dir (one uncapped level: 1 thread, no tool timeout, no
container limit) with timers around the steps of lamassemble and of the read
phasing, and work counters per lamassemble call:

  lam.lastdb / lam.lastal   the LAST subprocesses (wall)
  lam.parse                 _pairwise_alignments minus its subprocesses: MAF
                            parsing in Python
  lam.layout / lam.anchors  greedy layout and per-column anchors (Python)
  lam.mafft                 disttbfast
  lam.consensus             column-wise consensus (Python)
  ava.minimap2 / ava.parse  the phasing's all-vs-all and its pair parsing
  ref_phase                 reference-SNV phasing
  minimap_paf               reads to their consensus (final_consensus, re-adding)

Each step records wall and CPU seconds (this process plus its children), so
the share spent waiting (I/O, scheduling) is wall - CPU.

usage: profile_containers.py <work dir> <out prefix> --crids 2561,3464 \
    --bam BAM --reference FA --mat MAT [--tmp DIR] [--threads 1]
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import resource
import time
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path

from svirlpool.localassembly import (
    consensus,
    lamassemble,
    read_phasing,
    ref_read_phasing,
)
from svirlpool.util import util

STEPS: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0, 0])  # wall, cpu, calls
CALLS: list[dict] = []  # one per lamassemble call
_current: dict | None = None


def _cpu() -> float:
    s = resource.getrusage(resource.RUSAGE_SELF)
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return s.ru_utime + s.ru_stime + c.ru_utime + c.ru_stime


@contextmanager
def step(name: str):
    w0, c0 = time.perf_counter(), _cpu()
    try:
        yield
    finally:
        w, c = time.perf_counter() - w0, _cpu() - c0
        s = STEPS[name]
        s[0] += w
        s[1] += c
        s[2] += 1
        if _current is not None and name.startswith("lam."):
            _current[name] = _current.get(name, 0.0) + w


def timed(module, attr: str, name: str):
    f = getattr(module, attr)

    def wrapper(*a, **kw):
        with step(name):
            return f(*a, **kw)

    setattr(module, attr, wrapper)


def instrument() -> None:
    # LAST and MAFFT subprocesses run through _Deadline.run
    run = lamassemble._Deadline.run

    def deadline_run(self, cmd, **kw):
        tool = os.path.basename(str(cmd[0]))
        with step("lam." + ("mafft" if tool in ("disttbfast", "mafft") else tool)):
            return run(self, cmd, **kw)

    lamassemble._Deadline.run = deadline_run

    pw = lamassemble._pairwise_alignments

    def pairwise(params, scores, sequences, tmpdir, threads, deadline):
        before = {k: STEPS[k][0] for k in ("lam.lastdb", "lam.lastal")}
        cbefore = {k: STEPS[k][1] for k in ("lam.lastdb", "lam.lastal")}
        w0, c0 = time.perf_counter(), _cpu()
        out = pw(params, scores, sequences, tmpdir, threads, deadline)
        sub_w = sum(STEPS[k][0] - before[k] for k in before)
        sub_c = sum(STEPS[k][1] - cbefore[k] for k in cbefore)
        s = STEPS["lam.parse"]
        s[0] += time.perf_counter() - w0 - sub_w
        s[1] += _cpu() - c0 - sub_c
        s[2] += 1
        if _current is not None:
            _current["alignments"] = _current.get("alignments", 0) + len(out)
            _current["aligned_cols"] = _current.get("aligned_cols", 0) + sum(len(a[3]) for a in out)
            _current["lastal_passes"] = _current.get("lastal_passes", 0) + 1
        return out

    lamassemble._pairwise_alignments = pairwise

    lay = lamassemble._layout_of_seqs

    def layout(params, num_seqs, alns):
        with step("lam.layout"):
            r = lay(params, num_seqs, alns)
        if _current is not None:
            _current["kept"] = len(r[2])
        return r

    lamassemble._layout_of_seqs = layout

    anc = lamassemble._pairwise_anchors

    def anchors(*a, **kw):
        with step("lam.anchors"):
            r = list(anc(*a, **kw))
        if _current is not None:
            _current["anchors"] = _current.get("anchors", 0) + len(r)
        return iter(r)

    lamassemble._pairwise_anchors = anchors
    timed(lamassemble, "consensus_sequence", "lam.consensus")

    asm = lamassemble.assemble

    def assemble(sequences, *a, **kw):
        global _current
        _current = {"n_seqs": len(sequences), "bp": sum(len(s) for _, s in sequences),
                    "max_len": max((len(s) for _, s in sequences), default=0),
                    **kmer_mass([s for _, s in sequences])}
        w0, c0 = time.perf_counter(), _cpu()
        try:
            return asm(sequences, *a, **kw)
        finally:
            _current["wall"] = time.perf_counter() - w0
            _current["cpu"] = _cpu() - c0
            _current["crID"] = CONTAINER
            CALLS.append(_current)
            _current = None

    lamassemble.assemble = assemble

    timed(read_phasing, "run_ava", "ava.minimap2+best")
    timed(read_phasing, "parse_pair_both", "ava.parse")
    timed(read_phasing, "parse_pair", "ava.parse")
    timed(read_phasing, "snv_sites", "ava.snv_sites")
    timed(read_phasing, "sv_sites", "ava.sv_sites")
    timed(read_phasing, "phase_from_sites", "phase_from_sites")
    timed(ref_read_phasing, "phase_reads_reference", "ref_phase")
    timed(util, "align_reads_with_minimap_paf", "minimap_paf")
    timed(consensus, "orient_reads_to_reference", "ava.orient")
    timed(consensus, "process_consensus_container", "container")


CONTAINER = -1


def kmer_mass(seqs: list[str], k: int = 15) -> dict:
    """k-mer collision mass of the reads of one assembly: sum over k-mers of
    count^2, i.e. the seed hits of an all-vs-all alignment before LAST's -m
    cap; `repeat_x` = mass / (reads^2 x mean length), ~1 for a unique locus
    and ~ the tandem copies a read holds in a satellite."""
    from collections import Counter

    cnt: Counter = Counter()
    for s in seqs:
        s = s.upper()
        cnt.update(s[i : i + k] for i in range(len(s) - k + 1))
    mass = sum(v * v for v in cnt.values())
    n = len(seqs)
    mean_len = sum(map(len, seqs)) / max(1, n)
    return {"kmer_mass": mass, "repeat_x": mass / max(1.0, n * n * mean_len)}


def main():
    global CONTAINER
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("out")
    ap.add_argument("--crids", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--mat", required=True)
    ap.add_argument("--tmp", default=None)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--sites", default="tiered")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, filename=a.out + ".log",
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    instrument()
    work = Path(a.work)
    rows = []
    for crid in (int(x) for x in a.crids.split(",")):
        CONTAINER = crid
        STEPS.clear()
        w0, c0 = time.perf_counter(), _cpu()
        consensus.crs_containers_to_consensus(
            samplename=json.load(open(work / "config.json"))["samplename"],
            input=work / "crs_containers.db",
            copy_number_tracks=work / "copy_number_tracks.bed.gz",
            output=Path(a.out + f".{crid}.jsonl"),
            lamassemble_mat=a.mat,
            path_alignments=Path(a.bam),
            threads=a.threads,
            buffer_clipped_sequence=500,  # the CLI default (--buffer-clipped-sequence)
            consensus_method="lamassemble-onestrand",
            reference=Path(a.reference),
            escalation=[(a.threads, 10**6)],
            crIDs=[crid],
            tmp_dir_path=a.tmp,
            max_padding_size=100000,
            phasing_sites=a.sites,
            clustering_strategy="accurate",
        )
        row = {"crID": crid, "wall": time.perf_counter() - w0, "cpu": _cpu() - c0}
        for k, (w, c, n) in STEPS.items():
            row[f"{k}.wall"], row[f"{k}.cpu"], row[f"{k}.n"] = w, c, n
        rows.append(row)
        print(json.dumps({k: (round(v, 2) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
    import pandas as pd

    pd.DataFrame(rows).to_csv(a.out + ".steps.tsv", sep="\t", index=False)
    pd.DataFrame(CALLS).to_csv(a.out + ".lamcalls.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
