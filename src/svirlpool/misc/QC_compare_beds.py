# %%
import subprocess
import tempfile
from pathlib import Path
from shlex import split

from . import util

# %%
# set paths
path_crs_bed = Path("/data/hdd/ont-r10.32x/LRSV/full/QC/crs.bed")
path_truth_bed = Path(
    "/home/vinzenz/development/LRSV-detection/development/test/truth.bed"
)

# grep both files for dels and ins and intersect with bedtools
# grep dels in both files and save to two names tmp files
tmp_crs_dels = tempfile.NamedTemporaryFile()
tmp_truth_dels = tempfile.NamedTemporaryFile()
# grep both files and wirite to temp files
cmd_crs_dels = f"grep del {str(path_crs_bed)}"
cmd_truth_dels = f"grep del {str(path_truth_bed)}"
# execute both commands and write to temp files
subprocess.run(split(cmd_crs_dels), stdout=tmp_crs_dels)
subprocess.run(split(cmd_truth_dels), stdout=tmp_truth_dels)
# use bed tools to intersect both files to find true positives
cmd_intersect_positives = (
    f"bedtools intersect -wa -a {tmp_crs_dels.name} -b {tmp_truth_dels.name}"
)
# execute command and save output to df
df_crs_truth_dels_positives = util.execute_to_df(split(cmd_intersect_positives))
# negatives
cmd_intersect_negatives = (
    f"bedtools intersect -v -b {tmp_crs_dels.name} -a {tmp_truth_dels.name}"
)
# execute command and save output to df
df_crs_truth_dels_negatives = util.execute_to_df(split(cmd_intersect_negatives))

# now do the same for ins
tmp_crs_ins = tempfile.NamedTemporaryFile()
tmp_truth_ins = tempfile.NamedTemporaryFile()
# grep both files and wirite to temp files
cmd_crs_ins = f"grep ins {str(path_crs_bed)}"
cmd_truth_ins = f"grep ins {str(path_truth_bed)}"
# execute both commands and write to temp files
subprocess.run(split(cmd_crs_ins), stdout=tmp_crs_ins)
subprocess.run(split(cmd_truth_ins), stdout=tmp_truth_ins)
# use bed tools to intersect both files to find true positives
cmd_intersect_positives = (
    f"bedtools intersect -wa -a {tmp_crs_ins.name} -b {tmp_truth_ins.name}"
)
# execute command and save output to df
df_crs_truth_ins_positives = util.execute_to_df(split(cmd_intersect_positives))
# negatives
cmd_intersect_negatives = (
    f"bedtools intersect -v -b {tmp_crs_ins.name} -a {tmp_truth_ins.name}"
)
# execute command and save output to df
df_crs_truth_ins_negatives = util.execute_to_df(split(cmd_intersect_negatives))

# %%
# do it for the whole bed files regardless of sv type
#  intersect with bedtools
cmd_intersect_positives = (
    f"bedtools intersect -wa -a {str(path_crs_bed)} -b {str(path_truth_bed)}"
)
# execute command and save output to df
df_crs_truth_positives = util.execute_to_df(split(cmd_intersect_positives))
# negatives
cmd_intersect_negatives = (
    f"bedtools intersect -v -b {str(path_crs_bed)} -a {str(path_truth_bed)}"
)
# execute command and save output to df
df_crs_truth_negatives = util.execute_to_df(split(cmd_intersect_negatives))
# %%
