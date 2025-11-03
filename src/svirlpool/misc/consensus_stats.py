# load consensus stats
# glob all the files in consensus/log
# %%
import glob
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns

from . import util

# %%
Path_logs = (
    Path(
        "/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/consensus"
    )
    / "log"
)
QC_dir = Path("/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/QC")
path_crs = Path(
    "/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/HG002HG003.crs"
)
# iterate over all files
# and extract all stats
# example line:
# mcrID: 935 finished in 00:00:43 (43.614378213882446) on 15 reads and created 2 consensuses.
data = []
for path in glob.glob(str(Path_logs / "*.log")):
    with open(path, "r") as f:
        l = next(f).split(" ")
        data.append([int(l[1]), float(l[5][1:-1]), int(l[7]), int(l[11])])
# %%
df = pd.DataFrame(data, columns=["mcrID", "time", "reads", "num_consensus"])
df["log_time"] = np.log10(df["time"])
# add crs stats for mcrIDs
df["size"] = df["mcrID"].map(
    {
        cr.crID: cr.referenceEnd - cr.referenceStart
        for cr in util.yield_from_crs(input=path_crs)
    }
)
df["log_size"] = np.log10(df["size"])
# %%
print(df.corr())
# %%
fig = px.scatter(
    df,
    y="reads",
    x="time",
    hover_data=df.columns,
    color="log_size",
    marginal_x="histogram",
)
fig.write_image(QC_dir / "consensus.racon.scatter.reads.time.png")
fig.write_html(QC_dir / "consensus.racon.scatter.reads.time.html")
# %%
# plot df.corr() as sns heatmap, with title: correlation of consensus stats
plt.figure(figsize=(6, 6))
# correlation of consensus stats (all but not mcrID)
sns.heatmap(df.drop(columns=["mcrID"]).corr(), annot=True)
plt.title("correlation of consensus stats")
plt.savefig(QC_dir / "consensus.racon.heatmap.png")
# %%
print(df.sort_values(by="time", ascending=False).head(10))
# %%
