# create test data set
# %%
import subprocess
import tempfile
from pathlib import Path
from shlex import split

import pandas as pd

# %%

path_GiaB_T1_regions = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/benchmarks/HG002_SVs_Tier1_v0.6.bed"
)
path_T1_deletions = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/benchmarks/HG002_SVs_Tier1_v0.6.deletions.bed"
)
min_del_size = 50
sample_size = 300
seed = 123
output = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/del300/targets.bed"
)

# %%
df_GiaB_T1_regions = pd.read_csv(path_GiaB_T1_regions, sep="\t", header=None)
df_Giab_T1_deletions = pd.read_csv(path_T1_deletions, sep="\t", header=None)
bases_in_region = (df_GiaB_T1_regions[2] - df_GiaB_T1_regions[1]).sum()
# randomly pick 300 deletions that match min_del_size
# first subset min_del_size

# %%
mask = df_Giab_T1_deletions[3] > min_del_size
sampled_regions = df_Giab_T1_deletions.loc[mask, :].sample(
    n=sample_size, random_state=seed
)

# %%
# use bedtools intersect to get all T1 regions that overlap with selected deletions
tmp_deletions = tempfile.NamedTemporaryFile(suffix=".bed")
sampled_regions.to_csv(tmp_deletions, sep="\t", header=False, index=False)
# %%
cmd_sort = f"bedtools sort -i {tmp_deletions.name}"
cmd_intersect = (
    f"bedtools intersect -wa -a {str(path_GiaB_T1_regions)} -b {tmp_deletions.name}"
)
with open(output, "w") as f:
    subprocess.run(split(cmd_sort), stdout=f)
    subprocess.run(split(cmd_intersect), stdout=f)

# %%
# sum bases in ouput
df_output = pd.read_csv(output, sep="\t", header=None)
bases_in_output = (df_output[2] - df_output[1]).sum()
# %%
# again for this file:
path_file = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/minidels/minidel.targets.bed"
)
df_file = pd.read_csv(path_file, sep="\t", header=None)
bases_in_file = (df_file[2] - df_file[1]).sum()
# %%
