"""Sum of the all-vs-all phasing time (minimap2 -a ... phasing_reads -> 'read
phasing: status') over the consensus logs of a svirlpool workdir."""

import glob
import sys
from datetime import datetime

tot, n = 0.0, 0
for f in glob.glob(sys.argv[1] + "/consensus/*/consensus.batch_*.log"):
    if f.endswith(".diag.log"):
        continue
    t0 = None
    for line in open(f):
        if "read_phasing - INFO - minimap2 -a" in line:
            t0 = datetime.strptime(line[:23], "%Y-%m-%d %H:%M:%S,%f")
        elif "read phasing: status" in line and t0 is not None:
            tot += (
                datetime.strptime(line[:23], "%Y-%m-%d %H:%M:%S,%f") - t0
            ).total_seconds()
            n += 1
            t0 = None
print(f"{n} all-vs-all phasings, {tot / 3600:.2f} h")
