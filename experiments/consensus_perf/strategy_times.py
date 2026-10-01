"""Wall-clock vs CPU time of the svirlpool runs and their consensus stage.

  run wall          the svirlpool run job's wall clock (/usr/bin/time)
  run CPU           user + system time of the run and all its children
                    (/usr/bin/time; snakemake's cpu_time misses the nested jobs)
  consensus wall    span of the consensus stage: first batch start to last
                    batch end (batch end = output mtime, start = end - duration)
  consensus sum     sum of the batch jobs' wall times; batches run in parallel,
                    so this is "batch-seconds", not elapsed time

usage: strategy_times.py [variant ...]
"""

import glob
import os
import re
import sys

R = "/home/mayv_c/development/svp_improvements/results"
B = "/home/mayv_c/development/svp_improvements/benchmarks/svirlpool"
variants = sys.argv[1:] or ["cs_accurate", "cs_balanced", "cs_fast"]


def bench_s(path):
    lines = open(path).read().strip().splitlines()
    return float(dict(zip(lines[0].split("\t"), lines[-1].split("\t"), strict=True))["s"])


def time_v(path):
    t = open(path).read()
    user = float(re.search(r"User time \(seconds\): ([\d.]+)", t).group(1))
    sys_ = float(re.search(r"System time \(seconds\): ([\d.]+)", t).group(1))
    wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", t).group(1)
    parts = [float(x) for x in wall.split(":")]
    secs = sum(p * 60**i for i, p in enumerate(reversed(parts)))
    return secs, user + sys_


print("variant      sample  run wall  run CPU  consensus wall  consensus sum  batches")
for v in variants:
    for smp in ("HG002", "HG003", "HG004"):
        tv = f"{B}/run.{v}.20x.{smp}.time.txt"
        if not os.path.exists(tv):
            continue
        wall, cpu = time_v(tv)
        starts, ends, total = [], [], 0.0
        for b in glob.glob(f"{R}/{v}/20x/{smp}/work/benchmarks/consensus/consensus.batch_*.txt"):
            m = re.search(r"consensus\.batch_(\d+)\.(\d+)\.txt$", b)
            out = glob.glob(f"{R}/{v}/20x/{smp}/work/consensus/{m.group(2)}/consensus.batch_{m.group(1)}.jsonl")
            dur = bench_s(b)
            total += dur
            if out:
                end = os.path.getmtime(out[0])
                ends.append(end)
                starts.append(end - dur)
        span = max(ends) - min(starts) if ends else float("nan")
        print(f"{v:12s} {smp}  {wall:8.0f}  {cpu:7.0f}  {span:14.0f}  {total:13.0f}  {len(ends):7d}")
