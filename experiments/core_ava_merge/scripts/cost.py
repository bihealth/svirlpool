"""Cost numbers of an arm's sv-calling run (call.log, time.txt).

usage: cost.py ARM [ARM ...]
"""

import ast
import re
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path("/home/mayv_c/development/svp_mendel/arms")


def ts(line):
    return datetime.strptime(line[:23], "%Y-%m-%d %H:%M:%S,%f")


for arm in sys.argv[1:]:
    d = ROOT / arm
    log = (d / "call.log").read_text().splitlines()
    out = {"arm": arm}
    t_start = t_merge0 = t_merge1 = None
    for line in log:
        if not line[:4].isdigit():
            continue
        if t_start is None:
            t_start = ts(line)
        if "Processing chromosomes" in line or "CORE_ALN_MERGE|CORES|" in line and t_merge0 is None:
            t_merge0 = t_merge0 or ts(line)
        if "Merged " in line and "SVcomposites into" in line:
            t_merge1 = ts(line)
            out["merged"] = line.split("INFO - ")[1]
        if "CORE_ALN_MERGE|SUMMARY" in line:
            out["summary"] = ast.literal_eval(line.split("|mode=")[1].split("|", 1)[1])
            t_merge1 = ts(line)
        if "CORE_ALN_MERGE|CORES_LOADED" in line:
            out["cores_loaded"] = line.split("INFO - ")[1]
        if "VERTICAL_MERGE|seconds" in line:
            out["vertical_merge_seconds"] = line.split("=")[-1]
    if t_merge0 and t_merge1:
        out["merge_wall_from_log"] = (t_merge1 - t_merge0).total_seconds()
    tt = (d / "time.txt").read_text() if (d / "time.txt").exists() else ""
    for key in ("Elapsed (wall clock)", "Maximum resident set size", "User time", "System time"):
        m = re.search(re.escape(key) + r".*?: (\S+)\s*$", tt, re.M)
        if m:
            out[key] = m.group(1)
    for k, v in out.items():
        print(f"{k}: {v}")
    print()
