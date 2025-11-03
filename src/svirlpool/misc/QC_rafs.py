# %%

from pathlib import Path

import pandas as pd
import plotly.express as px
from tqdm import tqdm

from . import datatypes, util


# %%
def get_n_indels_signals(
    raf: datatypes.ReadAlignmentFragment, min_size: int, max_size: int
) -> int:
    return sum(
        [
            1
            for sv in raf.SV_signals
            if (
                min_size <= sv.size <= max_size
                and sv.sv_type < 2
                and (
                    raf.effective_interval[1]
                    <= sv.ref_start
                    <= raf.effective_interval[2]
                    or raf.effective_interval[1]
                    <= sv.ref_end
                    <= raf.effective_interval[2]
                )
            )
        ]
    )


def get_summed_indels_sizes(
    raf: datatypes.ReadAlignmentFragment, min_size: int, max_size: int
) -> int:
    return sum(
        [
            sv.size
            for sv in raf.SV_signals
            if (
                min_size <= sv.size <= max_size
                and sv.sv_type < 2
                and (
                    raf.effective_interval[1]
                    <= sv.ref_start
                    <= raf.effective_interval[2]
                    or raf.effective_interval[1]
                    <= sv.ref_end
                    <= raf.effective_interval[2]
                )
            )
        ]
    )


def get_median_indel_size(
    raf: datatypes.ReadAlignmentFragment, min_size: int, max_size: int
) -> int:
    return pd.Series(
        [
            sv.size
            for sv in raf.SV_signals
            if (
                min_size <= sv.size <= max_size
                and sv.sv_type < 2
                and (
                    raf.effective_interval[1]
                    <= sv.ref_start
                    <= raf.effective_interval[2]
                    or raf.effective_interval[1]
                    <= sv.ref_end
                    <= raf.effective_interval[2]
                )
            )
        ]
    ).median()


# %%
# load rafs
path_rafs = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/HG002.0.rafs.effective_intervals.tsv.gz"
)
rafs = list(tqdm(util.yield_from_raf(path_rafs)))
path_QC = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/QC"
)
# %%
# build dataframe with columns: effective_size, n_indels, sum_indel_bp, get_median_indel_size
columns = [
    "readname",
    "effective_size",
    "n_indels",
    "sum_indel_bp",
    "median_indel_size",
]
L = []
for raf in tqdm(rafs):
    L.append(
        [
            raf.read_name,
            raf.effective_interval[2] - raf.effective_interval[1],
            get_n_indels_signals(raf, 1, 50),
            get_summed_indels_sizes(raf, 1, 50),
            get_median_indel_size(raf, 1, 50),
        ]
    )
df = pd.DataFrame(L, columns=columns)
df["avg_indel_size"] = df["sum_indel_bp"] / df["n_indels"]
df["n_indels_per_1kbp"] = df["n_indels"] / df["effective_size"] * 1000
# %%
#  plot histograms of n_indels_per_1kbp of df with n_indels > 3 and effective_size >= 10000
fig = px.histogram(
    df.query("effective_size >= 10000").loc[:, ["n_indels"]], x="n_indels", nbins=1000
)
fig.write_html(path_QC / "QC_rafs.n_indels.html")
# %%
