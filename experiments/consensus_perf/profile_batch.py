"""Profile the consensus stage on real containers.

Runs ``crs_containers_to_consensus`` on a set of containers of an existing
svirlpool work dir and reports, per pipeline stage:

* wall time and Python CPU time, exclusive of nested stages,
* CPU time of the external tools the stage launched (``os.wait4`` rusage of
  every reaped child, which includes the child's own reaped descendants, e.g.
  lastal under lamassemble),
* a cProfile of the Python side (``<out>.pstats``).

usage:
  profile_batch.py <workdir> <out_prefix> (--batch N | --crIDs 1,2,3)
                   [--mode legacy|phased] [key=value consensus kwargs ...]
"""

from __future__ import annotations

import argparse
import cProfile
import json
import os
import pstats
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

from svirlpool.localassembly import consensus, read_phasing

# --------------------------------------------------------------------------
# stage timers (exclusive: time inside a nested stage counts to that stage)

STAGES: list[tuple[object, str]] = [
    (consensus, "process_consensus_container"),
    (consensus.read_cache_mod.ReadSequenceCache, "fetch_for_cr"),
    (consensus, "get_read_alignment_intervals_in_cr"),
    (consensus, "trim_reads"),
    (consensus, "summed_indel_distribution"),
    (consensus, "consensus_while_phasing"),
    (read_phasing, "run_ava"),
    (read_phasing, "parse_pair"),
    (read_phasing, "snv_sites"),
    (read_phasing, "recurrent_sites"),
    (read_phasing, "sv_sites"),
    (read_phasing, "pair_weights"),
    (read_phasing, "correlation_cluster"),
    (read_phasing, "refine_clusters"),
    (read_phasing, "low_quality_reads"),
    (read_phasing, "phase_reads"),
    (consensus, "consensus_from_clusters"),
    (consensus, "consensus_while_clustering_with_kmeans"),
    (consensus, "consensus_while_clustering"),
    (consensus, "assemble_consensus"),
    (consensus, "final_consensus"),
    (consensus, "add_unaligned_reads_to_consensuses_inplace"),
    (consensus, "add_cutread_alignments_to_consensus_inplace"),
    (consensus, "parse_ReadAlignmentSignals_from_alignment"),
    (consensus, "create_padding_for_consensus"),
]

wall = defaultdict(float)
pycpu = defaultdict(float)
calls = defaultdict(int)
child_cpu = defaultdict(float)  # (stage, tool) -> s
child_wall = defaultdict(float)
child_n = defaultdict(int)
_stack: list[list] = []  # [name, t_wall_start, t_cpu_start]


def _wrap(owner, name):
    fn = getattr(owner, name)

    def wrapper(*a, **kw):
        now, cpu = time.perf_counter(), time.process_time()
        if _stack:  # pause the parent stage
            top = _stack[-1]
            wall[top[0]] += now - top[1]
            pycpu[top[0]] += cpu - top[2]
        _stack.append([name, now, cpu])
        calls[name] += 1
        try:
            return fn(*a, **kw)
        finally:
            now, cpu = time.perf_counter(), time.process_time()
            top = _stack.pop()
            wall[name] += now - top[1]
            pycpu[name] += cpu - top[2]
            if _stack:
                _stack[-1][1], _stack[-1][2] = now, cpu

    wrapper.__wrapped__ = fn
    setattr(owner, name, wrapper)


for owner, name in STAGES:
    _wrap(owner, name)

# --------------------------------------------------------------------------
# child process accounting

_orig_init = subprocess.Popen.__init__


def _init(self, args, *a, **kw):
    self._perf_t0 = time.perf_counter()
    argv = args if isinstance(args, (list, tuple)) else str(args).split()
    tool = os.path.basename(str(argv[0])) if argv else "?"
    if tool == "samtools" and len(argv) > 1:
        tool = f"samtools {argv[1]}"
    self._perf_tool = tool
    self._perf_stage = _stack[-1][0] if _stack else "-"
    _orig_init(self, args, *a, **kw)


def _try_wait(self, wait_flags):
    try:
        pid, sts, ru = os.wait4(self.pid, wait_flags)
    except ChildProcessError:
        pid, sts = self.pid, 0
        ru = None
    if pid == self.pid and ru is not None:
        key = (self._perf_stage, self._perf_tool)
        child_cpu[key] += ru.ru_utime + ru.ru_stime
        child_wall[key] += time.perf_counter() - self._perf_t0
        child_n[key] += 1
    return (pid, sts)


subprocess.Popen.__init__ = _init
subprocess.Popen._try_wait = _try_wait


# --------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser()
    p.add_argument("workdir", type=Path)
    p.add_argument("out", type=Path)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--batch", type=int)
    g.add_argument("--crIDs")
    p.add_argument("--mode", default="phased")
    p.add_argument("--tmp", type=Path, default=None)
    p.add_argument("extra", nargs="*")
    a = p.parse_args()

    wd = a.workdir
    cfg = json.load(open(wd / "config.json"))
    if a.batch is not None:
        crIDs = consensus._load_crIDs_from_batch_tsv(
            wd / "consensus_batches.tsv", a.batch
        )
    else:
        crIDs = [int(x) for x in a.crIDs.split(",")]
    kw = {
        "samplename": cfg["samplename"],
        "input": wd / "crs_containers.db",
        "copy_number_tracks": wd / "copy_number_tracks.bed.gz",
        "output": a.out.with_suffix(".jsonl"),
        "lamassemble_mat": cfg["lamassemble_mat"],
        "path_alignments": Path(cfg["alignments"]),
        "threads": 1,
        "buffer_clipped_sequence": 500,
        "timeout": 20,
        "consensus_method": cfg["consensus_method"],
        "reference": Path(cfg["reference"]),
        "crIDs": crIDs,
        "tmp_dir_path": a.tmp,
        "max_padding_size": cfg["max_padding_size"],
        "max_copy_number_threshold": cfg["max_consensus_copy_number"],
        "clustering_mode": a.mode,
        "phasing_flank": cfg.get("phasing_flank", 10000),
        "phasing_fallback": cfg.get("phasing_fallback", "single"),
    }
    for x in a.extra:
        k, v = x.split("=", 1)
        kw[k] = eval(v)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    prof = cProfile.Profile()
    t0, c0 = time.perf_counter(), time.process_time()
    ru0 = os.times()
    prof.enable()
    consensus.crs_containers_to_consensus(**kw)
    prof.disable()
    t_total, c_total = time.perf_counter() - t0, time.process_time() - c0
    ru1 = os.times()
    prof.dump_stats(str(a.out.with_suffix(".pstats")))

    rows = []
    stages = sorted(set(wall) | {s for s, _ in child_cpu}, key=lambda s: -wall[s])
    for s in stages:
        ccpu = sum(v for (st, _), v in child_cpu.items() if st == s)
        rows.append((s, calls[s], wall[s], pycpu[s], ccpu))
    out = []
    out.append(
        f"containers {len(crIDs)}  wall {t_total:.1f} s  python cpu {c_total:.1f} s  "
        f"children cpu {ru1.children_user + ru1.children_system - ru0.children_user - ru0.children_system:.1f} s"
    )
    out.append(
        f"{'stage (exclusive)':45s} {'calls':>6s} {'wall':>8s} {'py_cpu':>8s} {'child_cpu':>9s}"
    )
    for s, n, w, c, cc in rows:
        out.append(f"{s:45s} {n:6d} {w:8.1f} {c:8.1f} {cc:9.1f}")
    out.append("")
    out.append(f"{'stage':45s} {'tool':18s} {'n':>6s} {'wall':>8s} {'cpu':>8s}")
    for (s, tool), v in sorted(child_cpu.items(), key=lambda kv: -child_wall[kv[0]]):
        out.append(
            f"{s:45s} {tool:18s} {child_n[(s, tool)]:6d} {child_wall[(s, tool)]:8.1f} {v:8.1f}"
        )
    out.append("")
    txt = "\n".join(out)
    print(txt)
    a.out.with_suffix(".txt").write_text(txt)
    st = pstats.Stats(prof, stream=sys.stdout)
    st.sort_stats("tottime").print_stats(30)


if __name__ == "__main__":
    main()
