"""Which containers of a svp_improvements variant escalated after a tool
timeout, and the level each one finished at.

usage: escalations.py <variant> [...]
"""

import glob
import re
import sys
from collections import Counter
from datetime import datetime

R = "/home/mayv_c/development/svp_improvements/results"

for v in sys.argv[1:]:
    finished = {}  # container -> (threads, seconds) it finished at, or "unresolved"
    tools = Counter()
    secs = {}  # container -> wall seconds, all levels
    for fn in glob.glob(f"{R}/{v}/20x/HG002/work/consensus/*/consensus.batch_*.log"):
        if fn.endswith(".diag.log"):
            continue
        cur = t0 = None
        for line in open(fn):
            m = re.search(r"representative crID (\d+)\)|Wrote \d+ container results", line)
            if m:
                t = datetime.strptime(line[:23], "%Y-%m-%d %H:%M:%S,%f")
                if cur is not None:
                    secs[cur] = (t - t0).total_seconds()
                cur, t0 = (int(m.group(1)) if m.group(1) else None), t
            m = re.search(
                r"Container (\d+): (.+) timed out \((\d+) thread\(s\), (\d+) s\); "
                r"escalating to (\d+) thread\(s\), (\d+) s",
                line,
            )
            if m:
                finished[int(m.group(1))] = (int(m.group(5)), int(m.group(6)))
                tools.update(m.group(2).split(", "))
                continue
            m = re.search(r"Container (\d+): (.+) timed out at the last level", line)
            if m:
                finished[int(m.group(1))] = "unresolved"
                tools.update(m.group(2).split(", "))
    levels = Counter(str(x) for x in finished.values())
    print(f"== {v}: {len(finished)} container(s) escalated")
    print("  finished at:", dict(sorted(levels.items())))
    print("  timeouts by tool (all levels):", dict(tools))
    esc_s = sum(secs[c] for c in finished if c in secs)
    print(
        f"  wall of the escalated containers: {esc_s:.0f} s of {sum(secs.values()):.0f} s"
    )
    for c, x in sorted(finished.items()):
        print(f"  container {c}: {x}, {secs.get(c, float('nan')):.0f} s")
