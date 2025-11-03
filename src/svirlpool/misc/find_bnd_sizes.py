# stats of all rafs from giab samples
# load HG002 + parents
# distribution of beak end sizes
# distribution of distances between indels
# %%
from pathlib import Path

import pandas as pd

from ..scripts import util

paths_rafs = {
    "HG002": Path(
        "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/giab/svirlpool/evaluate_bnd_sizes/HG002.rafs.tsv.gz"
    ),
    "HG003": Path(
        "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/giab/svirlpool/evaluate_bnd_sizes/HG003.rafs.tsv.gz"
    ),
    "HG004": Path(
        "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/giab/svirlpool/evaluate_bnd_sizes/HG004.rafs.tsv.gz"
    ),
}

# iterate rafs and fill stats (size of bndl, bndr, distance between indels)
stats = []
for sample, path in paths_rafs.items():
    for raf in util.yield_from_raf(input=path):
        end_last_sv = -1
        for signal in sorted(raf.SV_signals, key=lambda x: x.ref_start):
            if signal.sv_type >= 3:  # BNDL
                stats.append(
                    {"sample": sample, "type": "break end", "size": signal.size}
                )
            else:
                if end_last_sv != -1:
                    stats.append(
                        {
                            "sample": sample,
                            "type": "inter indel distance",
                            "size": signal.ref_start - end_last_sv,
                        }
                    )
                end_last_sv = signal.ref_end

# %%
# create a pandas dataframe
# each row is a value
# the columns are: type, sample, value
# type is one of: BNDL, BNDR, distance
# sample is the sample name
# value is the value
df = pd.DataFrame(stats)

# %%
import matplotlib.pyplot as plt

# create histograms
# rows are samples
# columns are types
fig, axs = plt.subplots(3, 2, figsize=(10, 15))
for i, sample in enumerate(paths_rafs.keys()):
    df_sample = df.query(f"sample == '{sample}'")
    for j, type_ in enumerate(["break end", "inter indel distance"]):
        df_type = df_sample.query(f"type == '{type_}' and size <= 500")
        axs[i, j].hist(
            df_type.loc[:, "size"], bins=100, alpha=0.5, label=sample, color=f"C{j}"
        )
        # axs[i,j].set_yscale("log")
        # axs[i,j].set_xscale("log")
        axs[i, j].set_title(type_)
        axs[i, j].legend()
plt.show()
# %%
