# script to plot all sv sizes and letter compositions
# %%
import logging as log

log.basicConfig(level=log.INFO)
import gzip
import json
from pathlib import Path

# %%
from typing import Generator

import cattrs
import pandas as pd
from plotly import express as px
from plotly import graph_objects as go
from tqdm import tqdm

from ..scripts import datatypes


def load_signals(
    path_signals: Path,
) -> Generator[datatypes.ExtendedSVsignal, None, None]:
    with gzip.open(path_signals, "rt") as f:
        for line in f:
            yield cattrs.structure(json.loads(line.strip()), datatypes.ExtendedSVsignal)


# %% load data ..
SVTYPES = {0: "INS", 1: "DEL", 3: "BNDL", 4: "BNDR"}


def df_from_signals(path_signals: Path) -> pd.DataFrame:
    return pd.DataFrame(
        data=[
            (
                SVTYPES[signal.sv_type],
                signal.size,
                signal.coverage,
                signal.strength,
                signal.repeatID > -1,
                signal.chrID,
            )
            for signal in tqdm(load_signals(path_signals))
        ],
        columns=["sv_type", "size", "GC content", "complexity", "repeat", "distance"],
    )


SAMPLE = "HG002"
path = Path(
    f"/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/giab/alignment_stats/hs37d5/{SAMPLE}.signals.json.gz"
)
df = df_from_signals(path)

# %%
# import pandas as pd
from plotly.subplots import make_subplots

# Define bins and labels for categorizing data
bins = [0, 20, 100, 1000, 10_000, float("inf")]
# generate labels for the bins
labels = []
for i in range(len(bins) - 1):
    labels.append(f"{bins[i]}-{bins[i+1]-1}")
labels[-1] = f"{bins[-2]}+"


# Create a new column in the DataFrame for binning
df["size_range"] = pd.cut(df["size"], bins=bins, labels=labels, right=True)

# Define subplots (2 rows for 'INS' & 'DEL', 4 columns for size ranges)
fig = make_subplots(
    rows=2,
    cols=len(labels),
    subplot_titles=[f"{sv} ({size})" for sv in ["INS", "DEL"] for size in labels],
    shared_xaxes=False,  # Ensuring independent x-axis scales
    shared_yaxes=False,  # Ensuring independent y-axis scales
)

# Map row index for 'INS' and 'DEL'
sv_mapping = {"INS": 1, "DEL": 2}
colors = {
    True: "rgb(255,150,30)",
    False: "rgb(30,150,255)",
}  # Assign colors for stacking

# automatically set x-axis limits based on size ranges
x_ranges = []
for i in range(len(bins) - 1):
    x_ranges.append((bins[i], bins[i + 1] - 1))
# set the last range to df["size"].max()
x_ranges[-1] = (bins[-2], df.query("sv_type == 'INS' | sv_type == 'DEL'")["size"].max())

# Iterate over size ranges and structural variant types to create stacked histograms
for col_idx, (size_range, (x_min, x_max)) in enumerate(zip(labels, x_ranges), start=1):
    for sv_type in ["INS", "DEL"]:
        for repeat_value in [True, False]:  # Create separate bars for True and False
            filtered_data = df[
                (df["size_range"] == size_range)
                & (df["sv_type"] == sv_type)
                & (df["repeat"] == repeat_value)
            ]

            fig.add_trace(
                go.Histogram(
                    x=filtered_data["size"],
                    name=f"{sv_type} {size_range} Repeat={repeat_value}",
                    marker=dict(color=colors[repeat_value], opacity=0.7),
                    histnorm="",
                ),
                row=sv_mapping[sv_type],
                col=col_idx,
            )
        fig.update_xaxes(range=[x_min, x_max], row=sv_mapping[sv_type], col=col_idx)


# Update layout for stacking and better visibility
fig.update_layout(
    title_text="Histogram of Structural Variants by Size Range",
    showlegend=True,
    height=600,
    width=1200,
    barmode="stack",  # Enables stacking of bars
)

fig.show()

# %% include sequence complexity. Plot size vs complexity
fig = px.scatter(
    df.query("(sv_type == 'INS' | sv_type == 'DEL') & size < 5000"),
    x="size",
    y="complexity",
    color="repeat",
    facet_col="sv_type",
    title="Size vs Complexity",
    labels={"size": "Size", "complexity": "Complexity", "repeat": "Repeat"},
)
fig.show()

# %% plot density for 300bp SVs
fig = px.scatter(
    df.query(
        "(sv_type == 'INS' | sv_type == 'DEL') & size >= 250 & size < 350 & distance < 10000 & distance > -1"
    ),
    x="size",
    y="distance",
    color="complexity",
    facet_row="repeat",
    facet_col="sv_type",
    title="Size vs Distance",
    labels={"size": "Size", "distance": "Distance", "repeat": "Repeat"},
)
fig.show()

fig = px.density_heatmap(
    df.query(
        "(sv_type == 'INS' | sv_type == 'DEL') & size >= 250 & size < 350 & distance < 10000 & distance > -1"
    ),
    x="size",
    y="distance",
    facet_col="sv_type",
    facet_row="repeat",
)
fig.show()


# %% plot density for 100 bp to 200 bp INS
fig = px.scatter(
    df.query(
        "(sv_type == 'INS') & size >= 100 & size < 200 & distance < 1000 & distance > -1"
    ),
    x="size",
    y="distance",
    color="complexity",
    facet_row="repeat",
    title="Size insertions vs Distance",
    labels={"size": "Size", "distance": "Distance", "repeat": "Repeat"},
)
fig.show()

fig = px.density_heatmap(
    df.query(
        "(sv_type == 'INS') & size >= 100 & size < 200 & distance < 1000 & distance > -1"
    ),
    x="size",
    y="distance",
    facet_row="repeat",
)
fig.show()

# %% plot density for 100 bp to 200 bp DEL
fig = px.scatter(
    df.query(
        "(sv_type == 'DEL') & size >= 100 & size < 200 & distance < 1000 & distance > -1"
    ),
    x="size",
    y="distance",
    color="complexity",
    facet_row="repeat",
    title="Size deletions vs Distance",
    labels={"size": "Size", "distance": "Distance", "repeat": "Repeat"},
)
fig.show()

fig = px.density_heatmap(
    df.query(
        "(sv_type == 'DEL') & size >= 100 & size < 200 & distance < 1000 & distance > -1"
    ),
    x="size",
    y="distance",
    facet_row="repeat",
)
fig.show()

# %%
# plot distribution of "distance" per sv_type (INS, DEL) with color=repeat

fig = px.histogram(
    df.query("(sv_type == 'INS' | sv_type == 'DEL') & distance >= 1 & distance < 1000"),
    x="distance",
    color="repeat",
    facet_row="sv_type",
)
fig.show()

fig = px.histogram(
    df.query(
        "(sv_type == 'INS' | sv_type == 'DEL') & distance >= 1000 & distance < 6000"
    ),
    x="distance",
    color="repeat",
    facet_row="sv_type",
)
fig.show()

fig = px.histogram(
    df.query(
        "(sv_type == 'INS' | sv_type == 'DEL') & distance >= 6000 & distance < 800000"
    ),
    x="distance",
    color="repeat",
    facet_row="sv_type",
)
fig.show()

fig = px.histogram(
    df.query("(sv_type == 'INS' | sv_type == 'DEL') & distance >= 800000"),
    x="distance",
    color="repeat",
    facet_row="sv_type",
)
fig.show()


# %%


# compute number of row that have repeat == True and distance > 4350 and sv_type == INS
d_INS_tail_rep = df.query(
    "repeat == True & distance > 4350 & sv_type == 'INS' & distance > 0"
).shape[0]
d_INS_tail_nr = df.query(
    "repeat == False & distance > 4350 & sv_type == 'INS' & distance > 0"
).shape[0]
d_INS_head_rep = df.query(
    "repeat == True & distance <= 4350 & sv_type == 'INS' & distance > 0"
).shape[0]
d_INS_head_nr = df.query(
    "repeat == False & distance <= 4350 & sv_type == 'INS' & distance > 0"
).shape[0]

d_DEL_tail_rep = df.query(
    "repeat == True & distance > 4350 & sv_type == 'DEL' & distance > 0"
).shape[0]
d_DEL_tail_nr = df.query(
    "repeat == False & distance > 4350 & sv_type == 'DEL' & distance > 0"
).shape[0]
d_DEL_head_rep = df.query(
    "repeat == True & distance <= 4350 & sv_type == 'DEL' & distance > 0"
).shape[0]
d_DEL_head_nr = df.query(
    "repeat == False & distance <= 4350 & sv_type == 'DEL' & distance > 0"
).shape[0]


# %%

# def run(args,**kwargs):
#     extract_all_signals(
#         path_alignments=args.alignments,
#         path_repeats=args.repeats,
#         path_reference=args.reference,
#         min_signal_size=args.min_signal_size,
#         path_output=args.output)
#     return


# def get_parser():
#     parser = argparse.ArgumentParser(description="Reads crs container objects and create Consensus objects that are written to the output database.")
#     parser.add_argument("-a","--alignments", type=Path, required=True, help="Path to the alignments file.")
#     parser.add_argument("-r","--repeats", type=Path, required=True, help="Path to the repeats file.")
#     parser.add_argument("-f","--reference", type=Path, required=True, help="Path to the reference file.")
#     parser.add_argument("-o","--output", type=Path, required=True, help="Path to the output gzipped json dumped file containing all signals.")
#     parser.add_argument("-s","--min_signal_size", type=int, default=6, required=False, help="Minimum signal size.")
#     parser.add_argument("--logfile", default=None, required=False, help="Path to the logfile.")
#     return parser


# def main():
#     parser = get_parser()
#     args = parser.parse_args()
#     if args.logfile:
#         logfile(str(args.logfile))
#     run(args)
#     return

# if __name__ == "__main__":
#     main()

# %%
