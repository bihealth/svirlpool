"""Consensus stage on a subsample of containers at added read error rates.

For every rate, runs ``consensus.crs_containers_to_consensus`` (the work
dir's config; ``--mode`` overrides its clustering mode) on the same random
subsample of containers, with ``added_error_rate`` = rate: the reads cut for
the phasing and for the assembly get that many additional random errors per
base
(``consensus.add_read_errors``). Writes, per rate, ``<out>/rate_<r>/``:
the consensus JSONL of every chunk and ``containers.tsv`` (per container:
wall seconds, attempts, timed-out tools, number of consensuses).

usage: run_rates.py <workdir> <out_dir> [--fraction 0.1] [--seed 0]
                    [--rates 0,0.01,...] [--procs 12] [--chunk 10]
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import random
import sqlite3
import time
from multiprocessing import Pool
from pathlib import Path

from svirlpool.localassembly import consensus, tool_timeouts

RATES = [0.0] + [round(0.01 * i, 2) for i in range(1, 11)]
_PROCESS_CONTAINER = consensus.process_consensus_container


def sample_crIDs(workdir: Path, fraction: float, seed: int) -> list[int]:
    con = sqlite3.connect(f"file:{workdir / 'crs_containers.db'}?mode=ro", uri=True)
    rows = []
    for k, v in con.execute("select crID, data from containers"):
        crs = json.loads(v)["crs"]
        rows.append((crs[0]["chr"], min(c["referenceStart"] for c in crs), k))
    ids = sorted(k for _, _, k in rows)
    chosen = set(random.Random(seed).sample(ids, round(fraction * len(ids))))
    # genomic order: consecutive containers share the read cache
    return [k for _, _, k in sorted(rows) if k in chosen]


def run_chunk(args) -> list[dict]:
    workdir, out, crIDs, rate, threads, mode = args
    cfg = json.load(open(workdir / "config.json"))
    if mode:
        cfg["consensus_clustering_mode"] = mode
    # the tools' stderr (lamassemble's "using X out of Y sequences") goes to
    # the chunk log, after a CONTAINER marker per attempt
    log_path = out.with_suffix(".log")
    log_f = open(log_path, "a", buffering=1)
    os.dup2(log_f.fileno(), 2)
    logging.basicConfig(
        stream=log_f,
        level=logging.WARNING,
        format="%(asctime)s %(levelname)s %(message)s",
        force=True,
    )
    rows: list[dict] = []
    inner = _PROCESS_CONTAINER  # pool workers are reused: wrap the original

    def timed(**kw):
        rep = min(kw["crs_dict"])
        log_f.write(f"CONTAINER {rep} threads={kw['threads']}\n")
        log_f.flush()
        watched = tool_timeouts._hits
        n0 = len(watched) if watched is not None else 0
        t0 = time.perf_counter()
        res = None
        try:
            res = inner(**kw)
            return res
        finally:
            hits = watched[n0:] if watched is not None else []
            rows.append(
                {
                    "crID": rep,
                    "threads": kw["threads"],
                    "timeout": kw["timeout"],
                    "secs": round(time.perf_counter() - t0, 3),
                    "timed_out": ",".join(hits) or "-",
                    "n_consensus": len(res[0]) if res else -1,
                }
            )

    consensus.process_consensus_container = timed
    consensus.crs_containers_to_consensus(
        samplename=cfg["samplename"],
        input=workdir / "crs_containers.db",
        copy_number_tracks=workdir / "copy_number_tracks.bed.gz",
        output=out,
        lamassemble_mat=cfg["lamassemble_mat"],
        path_alignments=Path(cfg["alignments"]),
        threads=threads,
        buffer_clipped_sequence=500,
        escalation=consensus.parse_escalation(cfg["consensus_escalation"]),
        consensus_method=cfg["consensus_method"],
        reference=Path(cfg["reference"]),
        crIDs=crIDs,
        max_padding_size=cfg["max_padding_size"],
        max_copy_number_threshold=cfg["max_consensus_copy_number"],
        clustering_mode=cfg["consensus_clustering_mode"],
        phasing_flank=cfg["phasing_flank"],
        phasing_fallback=cfg["phasing_fallback"],
        added_error_rate=rate,
    )
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("workdir", type=Path)
    p.add_argument("out_dir", type=Path)
    p.add_argument("--fraction", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--rates", default=",".join(map(str, RATES)))
    p.add_argument("--procs", type=int, default=12)
    p.add_argument("--chunk", type=int, default=10)
    p.add_argument("--threads", type=int, default=12, help="escalation ceiling")
    p.add_argument(
        "--mode", choices=["phased", "legacy"], help="override the clustering mode"
    )
    a = p.parse_args()

    crIDs = sample_crIDs(a.workdir, a.fraction, a.seed)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    (a.out_dir / "crIDs.txt").write_text("\n".join(map(str, crIDs)) + "\n")
    chunks = [crIDs[i : i + a.chunk] for i in range(0, len(crIDs), a.chunk)]
    print(f"{len(crIDs)} containers, {len(chunks)} chunks", flush=True)
    for rate in [float(x) for x in a.rates.split(",")]:
        d = a.out_dir / f"rate_{rate:.2f}"
        d.mkdir(exist_ok=True)
        jobs = [
            (a.workdir, d / f"chunk_{i:03d}.jsonl", c, rate, a.threads, a.mode)
            for i, c in enumerate(chunks)
        ]
        t0 = time.perf_counter()
        with Pool(a.procs) as pool:
            rows = [r for rs in pool.imap_unordered(run_chunk, jobs) for r in rs]
        wall = time.perf_counter() - t0
        rows.sort(key=lambda r: (r["crID"], r["threads"]))
        cols = list(rows[0])
        with open(d / "containers.tsv", "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")
        n_esc = len({r["crID"] for r in rows if r["timed_out"] != "-"})
        print(
            f"rate {rate:.2f}: wall {wall:.0f} s, container secs {sum(r['secs'] for r in rows):.0f}, "
            f"containers with a timeout {n_esc}",
            flush=True,
        )


if __name__ == "__main__":
    main()
