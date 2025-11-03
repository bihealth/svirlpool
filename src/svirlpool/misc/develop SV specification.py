# %% script to develop SV specification from BNDs in vcfs
from pathlib import Path

import vcfpy

# %%
vcf_file = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/giab/chm13_benchmark/svirlpool/HG002.33x/GRCh38.vcf.gz"
)
# %%
# 1) find DELs

records = []
breakends = []
for record in vcfpy.Reader.from_path(vcf_file):
    records.append(record)
    if isinstance(record.ALT[0], vcfpy.record.BreakEnd):
        breakends.append(record)

# %%
# filter breakends on chromosome 1
breakends_chr1 = [x for x in breakends if x.CHROM == "chr1"]
# create a bed file from the positions of the breakends plus a 2 kilobases padding
# merge them so that none overlap
# write to a file (output)
import tempfile

tmp_bed = tempfile.NamedTemporaryFile(delete=True, suffix=".bed")

with open(tmp_bed.name, "w") as f:
    for record in breakends_chr1:
        start = record.POS - 2000
        end = record.POS + 2000
        f.write(f"{record.CHROM[3:]}\t{start}\t{end}\n")


output = "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/bnd_test.bed"

import subprocess
from shlex import split

cmd_merge = f"bedtools merge -i {tmp_bed.name}"
with open(output, "w") as f:
    subprocess.run(split(cmd_merge), stdout=f)
# %%
