# %%
import glob
import pathlib

import pandas as pd
import plotly.express as px

# %%
wdir = pathlib.Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/"
)
"truvari/pb-hifi.20x/new_stats.txt"
L = []
for path in glob.glob(str(wdir / "truvari/*/new_stats.txt")):
    df = pd.read_csv(path, sep="\t", header=None)
    tech, cov = path.split("/")[-2].split(".")
    df.columns = ["metric", "score_corrected", "score"]
    df["tech"] = tech
    df["cov"] = int(cov[:-1])
    L.append(df)
df_all = pd.concat(L, axis=0)
# %%
df_all.index = range(df_all.shape[0])
selection = [row.metric in ["F1", "sens", "prec"] for i, row in df_all.iterrows()]
fig = px.line(
    df_all.loc[selection, :].sort_values(by=["cov"]),
    x="cov",
    y="score",
    facet_col="metric",
    color="tech",
    markers=True,
    title="naive metrics",
)
fig.show()

fig2 = px.line(
    df_all.loc[selection, :].sort_values(by=["cov"]),
    x="cov",
    y="score_corrected",
    facet_col="metric",
    color="tech",
    markers=True,
    title="corrected metrics",
)
fig2.show()
