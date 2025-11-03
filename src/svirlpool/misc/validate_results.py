# this script generates several bash scripts
# %%
from pathlib import Path

# %%
# =============================================================================
#  user input
path_base_script = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/src/scripts/misc/validate_results.sh"
)

path_script = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002HG003/minidels/jobscripts.sh"
)

WDIR = "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002HG003/minidels"
PREFIX = "minidels"
REF = "/fast/work/projects/cubit/18.12/static_data/reference/GRCh37/hs37d5/hs37d5.fa"
BENCHMARK = "benchmark/HG002_SVs_Tier1_v0.6.deletions.vcf.gz"

args0 = [2, 4, 6, 8, 12]
args1 = [12, 24, 42, 84]

# args:
# WDIR=$1
# PREFIX=$2
# REF=$3
# BENCHMARK=$4

# arg_O0=$5
# arg_O1=$6

# argID=$7

# write bash file with script calls to validate_results.sh
with open(path_script, "w") as f:
    print("#!/bin/bash", file=f)
    for arg0 in args0:
        for arg1 in args1:
            print(
                "sbatch",
                str(path_base_script),
                WDIR,
                PREFIX,
                REF,
                BENCHMARK,
                str(arg0),
                str(arg1),
                str(arg0) + "_" + str(arg1),
                sep=" ",
                file=f,
            )

# %%
