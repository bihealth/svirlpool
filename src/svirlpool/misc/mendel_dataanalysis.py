# %%
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# %%
path_svirlpool_platinum_hs1 = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/platinum/svirlpool/chm13/mendelian/mendel.30x.tsv"
)
path_svirlpool_platinum_hg18 = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/platinum/svirlpool/GRCh38/mendelian/mendel.30x.tsv"
)
path_svirlpool_giab_hg19 = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/giab/svirlpool/hs37d5/mendelian/mendel.30x.tsv"
)
path_svirlpool_giab_hs1 = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/giab/svirlpool/chm13/mendelian/mendel.30x.tsv"
)

path_sniffles_platinum_hs1 = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/platinum/sniffles/chm13/mendelian/mendel.tsv"
)
path_sniffles_platinum_hg18 = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/platinum/sniffles/GRCh38/mendelian/mendel.tsv"
)
path_sniffles_giab_hs1 = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/giab/sniffles/hs1/trio/mendel.tsv"
)
path_sniffles_giab_hg19 = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/giab/sniffles/hs37d5/trio/mendel.tsv"
)

path_output_percentages = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/figures/raw/mendelian_errors_percentages.svg"
)
path_output_counts = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/figures/raw/mendelian_errors_counts.svg"
)


def parse_df(path: Path, tool: str, ref: str) -> pd.DataFrame:
    # parse the dataframe and add columns tool and ref
    df = pd.read_csv(path, sep="\t")
    df["tool"] = tool
    df["ref"] = ref
    return df


# load all data to a dataframe
# each file has columns:
#  sample (str), status (str), count (int), percentage (float)
df = pd.concat(
    [
        parse_df(path_svirlpool_platinum_hs1, "svirlpool", "hs1"),
        parse_df(path_svirlpool_platinum_hg18, "svirlpool", "hg38"),
        parse_df(path_svirlpool_giab_hg19, "svirlpool", "hg19"),
        parse_df(path_svirlpool_giab_hs1, "svirlpool", "hs1"),
        parse_df(path_sniffles_platinum_hs1, "sniffles", "hs1"),
        parse_df(path_sniffles_platinum_hg18, "sniffles", "hg38"),
        parse_df(path_sniffles_giab_hg19, "sniffles", "hg19"),
        parse_df(path_sniffles_giab_hs1, "sniffles", "hs1"),
    ]
)
# re-calculate percentages to avoid rounding errors in the precomputed files
group_totals = df.groupby(["sample", "tool", "ref"])["count"].transform("sum")
df["percentage"] = (df["count"] / group_totals) * 100

# %%
# create a bar plot for percentage and put svirlpool and sniffles side by side
custom_colors = {
    "consistent": "rgb(89, 136, 107)",  # green
    "inconsistent": "rgb(192, 85, 85)",  # red
    "missing": "rgb(255, 200, 92)",  # yellow
    "incomplete": "rgb(255, 248, 193)",  # yellow
}

pattern_map = {
    "consistent": "",  # Solid fill
    "incomplete": "/",  # Diagonal stripe
    "missing": "\\",  # Opposite diagonal stripe
    "inconsistent": "x",  # Crosshatch
}

status_order = ["consistent", "incomplete", "missing", "inconsistent"]

# Step 1: Create composite label
df["sample_tool"] = df["sample"] + " - " + df["tool"]

# Step 2: Build x-axis order with visual spacers (not in the data)
samples = df["sample"].unique()
tools = df["tool"].unique()

spaced_order = []
for sample in samples:
    for tool in tools:
        spaced_order.append(f"{sample} - {tool}")
    spaced_order.append(f"spacer_{sample}")  # spacer just in the axis order

# =============================================================================
#  PLOT PERCENTAGES
# =============================================================================


# Step 3: Plot with visual spacers
fig = px.bar(
    df,
    x="sample_tool",  # only real values here!
    y="percentage",
    color="status",
    pattern_shape="status",
    color_discrete_map=custom_colors,
    barmode="stack",
    facet_col="ref",
    category_orders={"sample_tool": spaced_order, "status": status_order},
    title="Mendelian errors for Svirlpool and Sniffles",
)

# Step 4: Hide spacer x-ticks
fig.update_layout(uniformtext_minsize=8, uniformtext_mode="hide")
fig.update_xaxes(tickvals=[x for x in spaced_order if not x.startswith("spacer_")])
fig.update_layout(
    xaxis_title="Sample",
    yaxis_title="Percentage of Mendelian errors",
    legend_title="Tool",
    title_x=0.5,
)
fig.update_xaxes(
    title_text="Sample",
    tickangle=-90,
)
fig.update_traces(
    textfont_size=12,
    textfont_color="black",
    textposition="outside",
)
fig.update_layout(
    margin=dict(l=20, r=20, t=40, b=20),
    height=600,
    width=1000,
)
# fig.show()
# write fig to svg
fig.write_image(path_output_percentages, format="svg")

# %%

# =============================================================================
#  PLOT COUNTS
# =============================================================================

# Step 3: Plot with visual spacers
fig = px.bar(
    df,
    x="sample_tool",  # only real values here!
    y="count",
    color="status",
    pattern_shape="status",
    color_discrete_map=custom_colors,
    barmode="stack",
    facet_col="ref",
    category_orders={"sample_tool": spaced_order, "status": status_order},
    title="Mendelian errors for Svirlpool and Sniffles",
)

# Step 4: Hide spacer x-ticks
fig.update_layout(uniformtext_minsize=8, uniformtext_mode="hide")
fig.update_xaxes(tickvals=[x for x in spaced_order if not x.startswith("spacer_")])
fig.update_layout(
    xaxis_title="Sample",
    yaxis_title="Counts of Mendelian errors",
    legend_title="Tool",
    title_x=0.5,
)
fig.update_xaxes(
    title_text="Sample",
    tickangle=-90,
)
fig.update_traces(
    textfont_size=12,
    textfont_color="black",
    textposition="outside",
)
fig.update_layout(
    margin=dict(l=20, r=20, t=40, b=20),
    height=600,
    width=1000,
)
# fig.show()
# write fig to svg
fig.write_image(path_output_counts, format="svg")

# %%
# =============================================================================
#  PLOT TABLES
# =============================================================================

tables = []
refs = ["hs1", "hg38", "hg19"]
for ref in refs:
    rows = []
    samples = df.query("ref == @ref")["sample"].unique()
    tools = ["svirlpool", "sniffles"]
    stati = ["consistent", "incomplete", "missing", "inconsistent"]
    for sample in samples:
        for tool in tools:
            row = [sample, tool]
            for status in stati:
                row.append(
                    float(
                        df.query(
                            "sample == @sample and tool == @tool and status == @status and ref == @ref"
                        )["percentage"].values[0]
                    )
                )
            rows.append(row)
    # create a dataframe
    columns = ["sample", "tool"] + list(stati)
    table_df = pd.DataFrame(rows, columns=columns)
    # round all values to 2 decimals and save as string
    table_df[stati] = table_df[stati].map(lambda x: f"{x:.2f}")
    tables.append(table_df)
# %%
# render each table as an svg and format all floats to 2 decimal places
from plotly.subplots import make_subplots

# from plotly.io import to_image
# import plotly.io as pio
# pio.renderers.default = "svg"

for i, table in enumerate(tables):
    path_output = Path(
        f"/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/figures/raw/mendelian_errors_table_{i+1}.svg"
    )

    # Color the column headers based on the status columns
    header_colors = []
    for col in table.columns:
        # If column is one of the status columns, apply corresponding color from custom_colors
        if col in custom_colors:
            header_colors.append(custom_colors[col])
        else:
            # Default color for non-status columns
            header_colors.append("lightgrey")

    fig = make_subplots(rows=1, cols=1)

    # Adjust alignment for specific columns (right-align status columns)
    alignments = ["left"] * len(
        table.columns
    )  # Default to left alignment for all columns
    for col_idx, col in enumerate(table.columns):
        if col in custom_colors:
            alignments[col_idx] = "right"  # Right-align status columns

    fig.add_trace(
        go.Table(
            header=dict(
                values=list(table.columns),
                fill_color=header_colors,  # Apply header colors here
                font=dict(color="black", size=12),
                align="center",
            ),  # Align headers to the center (you can change this if needed)
            cells=dict(
                values=[table[col] for col in table.columns],
                fill_color=["white"] * len(table.columns),  # Keep cells white
                align=alignments,
            ),  # Apply alignment here
            columnwidth=[0.15, 0.13, 0.185, 0.175, 0.175, 0.185],
        )
    )

    fig.update_layout(title_text=f"Table {i+1} - {refs[i]}")
    fig.update_layout(margin=dict(l=20, r=20, t=40, b=20), height=400, width=500)

    fig.write_image(path_output, format="svg")


# %%
