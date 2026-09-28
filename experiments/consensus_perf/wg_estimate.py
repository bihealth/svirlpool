"""Crude whole-genome extrapolation of an svp_improvements run.

The benchmark covers 22 blocks, 151.7 Mb (5% of hg38, ~5% of the truth SVs),
so every stage that works per region is scaled by genome size / block size.
Per stage: the snakemake benchmarks of the nested svirlpool workflow (wall
seconds, peak RSS). The consensus batches run their tools single-threaded
(except escalations), so their wall time is their CPU time; snakemake's
sampled cpu_time misses the short-lived tools. The run total CPU is from
/usr/bin/time (user + sys over all waited-for descendants).

usage: wg_estimate.py <variant> [...]
"""

import glob
import re
import sys
from collections import defaultdict
from pathlib import Path

R = Path("/home/mayv_c/development/svp_improvements")
BLOCKS_BP = 151_700_000
GENOME_BP = 3_099_734_149  # GRCh38 primary assembly incl. N
SCALE = GENOME_BP / BLOCKS_BP


def run_time(variant: str) -> tuple[float, float, float]:
    """total CPU s, wall s, peak RSS GB of `svirlpool run`"""
    t = (R / f"benchmarks/svirlpool/run.{variant}.20x.HG002.time.txt").read_text()
    user = float(re.search(r"User time \(seconds\): ([\d.]+)", t).group(1))
    sys_ = float(re.search(r"System time \(seconds\): ([\d.]+)", t).group(1))
    w = re.search(
        r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", t
    ).group(1)
    parts = [float(x) for x in w.split(":")]
    wall = sum(p * 60**i for i, p in enumerate(reversed(parts)))
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", t).group(1))
    return user + sys_, wall, rss / 1e6


def stages(variant: str):
    """stage -> [jobs, summed wall s, max wall s, max RSS MB]"""
    out = defaultdict(lambda: [0, 0.0, 0.0, 0.0])
    base = R / f"results/{variant}/20x/HG002/work/benchmarks"
    for fn in glob.glob(f"{base}/**/*.txt", recursive=True):
        name = Path(fn).stem
        stage = "consensus batches" if name.startswith("consensus.batch_") else name
        with open(fn) as f:
            f.readline()
            x = f.readline().split("\t")
        s, rss = float(x[0]), float(x[2]) if x[2] != "-" else 0.0
        o = out[stage]
        o[0] += 1
        o[1] += s
        o[2] = max(o[2], s)
        o[3] = max(o[3], rss)
    return out


print(f"scale: {GENOME_BP / 1e9:.2f} Gb / {BLOCKS_BP / 1e6:.1f} Mb = x{SCALE:.1f}\n")
for v in sys.argv[1:]:
    cpu, wall, rss = run_time(v)
    print(f"== {v}")
    print(
        f"  run on the blocks: CPU {cpu:.0f} s, wall {wall / 60:.1f} min (16 threads), "
        f"peak RSS {rss:.1f} GB"
    )
    st = stages(v)
    print(
        f"  {'stage':50s} {'jobs':>5s} {'wall s':>8s} {'max job s':>9s} {'RSS GB':>7s}"
    )
    for name, (n, s, mx, r) in sorted(st.items(), key=lambda kv: -kv[1][1]):
        print(f"  {name:50s} {n:5d} {s:8.0f} {mx:9.0f} {r / 1e3:7.1f}")
    cons = st["consensus batches"][1]
    print(
        f"  whole genome (x{SCALE:.1f}): CPU {cpu * SCALE / 3600:.0f} h "
        f"(consensus {cons * SCALE / 3600:.0f} h, {st['consensus batches'][0] * SCALE:.0f} batches); "
        f"wall at full use of 16 / 64 cores: "
        f"{cpu * SCALE / 16 / 3600:.1f} / {cpu * SCALE / 64 / 3600:.1f} h"
    )
    print()
