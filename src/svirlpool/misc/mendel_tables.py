# %% load tables
from pathlib import Path

import pandas as pd

paths = {
    "HS1": {
        "svirlpool": Path(
            "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/platinum/svirlpool/chm13/mendelian/mendel.tsv"
        ),
        "sniffles": Path(
            "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/platinum/sniffles/chm13/mendelian/mendel.tsv"
        ),
    },
    "HG38": {
        "svirlpool": Path(
            "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/platinum/svirlpool/GRCh38/mendelian/mendel.tsv"
        ),
        "sniffles": Path(
            "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/platinum/sniffles/GRCh38/mendelian/mendel.tsv"
        ),
    },
}

##% load tables
tables = {}
for ref, tools in paths.items():
    tables[ref] = {}
    for tool, path in tools.items():
        tables[ref][tool] = pd.read_csv(path, sep="\t")
        tables[ref][tool]["Tool"] = tool
        tables[ref][tool]["Ref"] = ref
# concat tables to one table
tables = pd.concat(
    [tables[ref][tool] for ref in tables for tool in tables[ref]], axis=0
)
# %% crete figure
from plotly import graph_objects as go
from plotly.subplots import make_subplots

df = tables.copy()

status_colors = {
    "missing": "rgba(110, 110, 80, 1)",
    "consistent": "rgba(50, 150, 150, 1)",
    "inconsistent": "rgba(150, 0, 0, 1)",
    "incomplete": "rgba(170, 170, 160, 1)",
}
custom_status_order = ["consistent", "incomplete", "missing", "inconsistent"]
df["status"] = pd.Categorical(
    df["status"], categories=custom_status_order, ordered=True
)
df["percentage"] = df["percentage"] * 100.0
df["sample_tool"] = df["sample"].astype(str) + "_" + df["Tool"].astype(str)

# Add extra spacing between categories
# df["sample"] = df["sample"].replace({
#     s: s + " " * (i % 2)  # Add spaces to every other sample name
#     for i, s in enumerate(df["sample"].unique())
# })

df = df.sort_values(["sample", "status"])


# add two columns to df:
# 1) cumulative counts for each sample in the order of status_colors.keys()
# 2) cumulative percentages for each sample in the order of status_colors.keys()
df["cumulative_count"] = (
    df.groupby(["sample", "Ref", "Tool"])["count"].cumsum() - df["count"]
)
df["cumulative_percentage"] = (
    df.groupby(["sample", "Ref", "Tool"])["percentage"].cumsum() - df["percentage"]
)

# Get unique references
unique_refs = df["Ref"].unique()


# Precompute x-axis values for spacing
dict_x_coords = {}
x_labels = []  # To store x-tick labels
x_positions = []  # To store x-tick positions

for i, sample in enumerate(df["sample"].unique()):
    for j, tool in enumerate(df["Tool"].unique()):
        # Assign x-position
        x_pos = i * 8 + j * 3
        dict_x_coords[(sample, tool)] = (x_pos,)

        # Store unique x-tick positions and labels
        if j == 0:  # Only store label once per sample (avoid duplicating per tool)
            x_labels.append(sample)
            x_positions.append(x_pos + 1.5)  # Center label between bars

for metric in ["count", "percentage"]:
    # Create subplot structure (one row per Ref, shared x-axis)
    fig = make_subplots(
        rows=1,
        cols=len(unique_refs),
        shared_xaxes=True,
        subplot_titles=[f"Reference: {ref}" for ref in unique_refs],
    )
    # Add 'count' bars (left y-axis)
    for row_idx, Ref in enumerate(unique_refs, start=1):  # Rows start at 1
        for sample in df["sample"].unique():
            for tool in df["Tool"].unique():
                for status in reversed(custom_status_order):
                    df_subset = df.query(
                        "Ref == @Ref and status == @status and Tool == @tool and sample == @sample"
                    )
                    # if metric == "count":
                    #     print(df_subset[metric])
                    #     print(df_subset["cumulative_count"])
                    print(df_subset.shape)
                    fig.add_trace(
                        go.Bar(
                            x=dict_x_coords[(sample, tool)],
                            y=df_subset[metric],
                            width=2.0,
                            name=None,
                            offsetgroup=sample + tool,
                            marker_color=status_colors[status],
                            base=df_subset[f"cumulative_{metric}"],
                            legendgroup=status,
                            showlegend=False,
                        ),
                        row=1,
                        col=row_idx,
                    )

    # Apply x-tick labels to all subplots
    for i in range(1, len(unique_refs) + 1):
        fig.update_xaxes(
            tickvals=x_positions,  # Position of labels
            ticktext=x_labels,  # Actual sample names
            tickangle=45,
            tickfont=dict(size=14),  # Control font size
            title_text="",
            row=1,
            col=i,
        )

    # Update layout for better spacing
    fig.update_layout(
        height=300 * len(unique_refs),  # Adjust height based on number of subplots
        barmode="stack",  # Keep stacking behavior
        bargap=0.1,  # Space between bars within the same group
        bargroupgap=0.3,  # Space between different groups of bars
        xaxis_title="",
        yaxis_title=metric,
    )

    for status in reversed(custom_status_order):
        color = status_colors[status]
        fig.add_trace(
            go.Bar(
                x=[None],  # Invisible data point
                y=[None],  # Invisible data point
                marker=dict(color=color),
                name=status,  # Legend label
                legendgroup=status,  # Group legends
                showlegend=True,  # Ensure it appears in the legend
            )
        )

    annotations = []
    for i in range(len(x_positions)):
        for j, tool in enumerate(["svirlpool", "sniffles"]):
            annotations.append((i * 8 + j * 3, tool))

    # Add annotations to the figure
    for xpos, ann in annotations:
        for col in range(1, len(unique_refs) + 1):
            fig.add_annotation(
                x=xpos,
                y=9 if metric == "percentage" else 8000,
                text=ann,
                textangle=-90,
                showarrow=False,  # Display arrow
                font=dict(size=14, color="white"),  # Font color
                align="center",
                bgcolor="rgba(255,255,255,0.0)",  # Semi-transparent background
                bordercolor="rgba(255,255,255,0.0)",  # Transparent border
                borderwidth=0,
                row=1,
                col=col,
            )

    fig.show()
    path_svg = Path(
        f"/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/platinum/platinum_mendel_{metric}.svg"
    )
    path_pdf = Path(
        f"/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/platinum/platinum_mendel_{metric}.pdf"
    )
    fig.write_image(path_svg, format="svg")
    fig.write_image(path_pdf, format="pdf")

# %%
