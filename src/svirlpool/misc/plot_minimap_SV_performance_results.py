# %%
import json
from glob import glob
from pathlib import Path

import pandas as pd
import plotly.express as px

# %%
p = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002HG003/minidels/results.txt"
)
summary_files = glob(str(p.parent / "truvari.*" / "summary.json"))

L = []
for summary in summary_files:
    json_data = open(summary)
    data = json.load(json_data)
    json_data.close()
    L.append(
        [
            data["f1"],
            data["precision"],
            data["recall"],
            data["TP-base"],
            data["FP"],
            data["FN"],
        ]
    )
df = pd.DataFrame(L, columns=["F1", "precision", "recall", "TP-base", "FP", "FN"])
df["params"] = [Path(x).parent.name.split(".")[1] for x in summary_files]
# split params into two columns
df[["O0", "O1"]] = df["params"].str.split("_", expand=True)

# %%
# plot 2d scatter of precision and recall
fig = px.scatter(
    df,
    x="precision",
    y="recall",
    hover_data=["params"],
    color="F1",
    title="Precision vs. Recall for different minimap2 parameters",
    color_continuous_scale=px.colors.sequential.matter,
)
fig.update_traces(
    marker=dict(size=8, symbol="diamond", line=dict(width=2, color="DarkSlateGrey")),
    selector=dict(mode="markers"),
)
fig.show()
fig.write_html(str(p.parent / "results_precision_recall_F1.html"))
# make another plot for 'TP-base','FP','FN'
fig2 = px.scatter(
    df,
    x="TP-base",
    y="FP",
    hover_data=["params"],
    color="F1",
    title="TP-base vs. FP for different minimap2 parameters",
    color_continuous_scale=px.colors.sequential.matter,
)
fig2.update_traces(
    marker=dict(size=8, symbol="diamond", line=dict(width=2, color="DarkSlateGrey")),
    selector=dict(mode="markers"),
)
fig2.show()
fig2.write_html(str(p.parent / "results_TP-base_FP_F1.html"))

# %%
