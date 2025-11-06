# %%
import argparse
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go


def create_plot(df: pd.DataFrame, title: str) -> go.Figure:
    status_colors = {
        "consistent": "rgba(93, 156, 89, 200)",
        "inconsistent": "rgba(223, 46, 56, 200)",
    }
    custom_status_order = ["consistent", "inconsistent"]
    df["status"] = pd.Categorical(
        df["status"], categories=custom_status_order, ordered=True
    )
    df = df.sort_values(["file_group", "sample", "status"])
    # percentage is float, but in the 0.0 to 1.0 range, so it needs to be scaled up to 0-100
    df["percentage"] = df["percentage"] * 100.0

    # Create combined x-axis labels (file_group + sample)
    df["x_label"] = df["file_group"].astype(str) + " - " + df["sample"]

    # Make x_label categorical with explicit ordering to preserve file_group order
    x_order = df["x_label"].unique()  # This preserves the order from file_group sorting
    df["x_label"] = pd.Categorical(df["x_label"], categories=x_order, ordered=True)

    # Calculate cumulative percentages for stacking
    df["cumulative_percentage"] = (
        df.groupby("x_label", observed=True)["percentage"].cumsum() - df["percentage"]
    )

    fig = go.Figure()

    # Add percentage bars with count annotations
    for status in reversed(custom_status_order):
        df_subset = df[df["status"] == status]
        fig.add_trace(
            go.Bar(
                x=df_subset["x_label"],
                y=df_subset["percentage"],
                name=f"{status}",
                marker_color=status_colors[status],
                text=[
                    f"{val:.1f}" for val in df_subset["percentage"]
                ],  # Format to 1 decimal place
                textposition="inside",  # Position text inside the bar
                textangle=0,  # Prevent rotation
                textfont={"color": "white", "size": 11},  # Fixed font size
                insidetextanchor="middle",  # Anchor text in the middle
                constraintext="none",  # Don't constrain text size
                legendgroup=status,
                base=df_subset["cumulative_percentage"],  # Stack bars
            )
        )

    # Layout updates
    fig.update_layout(
        title=title,
        barmode="stack",  # Stacked bars
        yaxis={"title": "Percentage (%)", "range": [0, 105]},
        xaxis={"title": "File Group - Sample", "tickangle": -45},
        legend={
            "x": 1.02,
            "y": 1,
            "xanchor": "left",
            "yanchor": "top",
            "font": {"size": 12, "family": "Arial", "color": "black"},
        },
    )
    return fig


# %%
def create_table(df: pd.DataFrame) -> pd.DataFrame:
    """Create a matrix table with samples as rows and file_group-status combinations as columns."""
    # Sort by file_group and sample to match plot order
    df = df.sort_values(["file_group", "sample"])

    # Format percentage to 2 decimal places (multiply by 100 to show as percentages)
    df["percentage_formatted"] = (df["percentage"] * 100.0).round(2)

    # Create combined row identifier (file_group - sample) to match plot x-axis
    df["row_label"] = (
        df["file_group"].cat.categories[df["file_group"].cat.codes].astype(str)
        + " - "
        + df["sample"]
    )

    # Create combined column name for status only (since row already has file_group)
    df["column_name"] = df["status"].cat.categories[df["status"].cat.codes].astype(str)

    # Create ordered list of column names (just statuses)
    statuses = df["status"].cat.categories.tolist()

    # Make column_name categorical with explicit ordering
    df["column_name"] = pd.Categorical(
        df["column_name"], categories=statuses, ordered=True
    )

    # Create row order based on file_group and sample sorting (unique x_labels)
    row_order = (
        df[["file_group", "sample", "row_label"]]
        .drop_duplicates()["row_label"]
        .tolist()
    )

    # Pivot to create matrix: row_label as rows, status as columns
    table = df.pivot_table(
        index="row_label",
        columns="column_name",
        values="percentage_formatted",
        aggfunc="first",
    )

    # Reset index to make row_label a regular column
    table = table.reset_index()

    # Reorder rows to match row_order
    table["row_label"] = pd.Categorical(
        table["row_label"], categories=row_order, ordered=True
    )
    table = table.sort_values("row_label").reset_index(drop=True)

    # Rename row_label column to something more descriptive
    table = table.rename(columns={"row_label": "file_group - sample"})

    return table


# %%
def plot_mendelian(
    inputs: list[Path], names: list[str], output: Path, title: str
) -> pd.DataFrame:
    # make sure the output has the right extension (svg or png)
    if output.suffix not in [".svg", ".png"]:
        raise ValueError(
            f"Output file must have extension .svg or .png. Got {output.suffix} instead."
        )

    # Validate that names match inputs
    if names and len(names) != len(inputs):
        raise ValueError(
            f"Number of names ({len(names)}) must match number of input files ({len(inputs)})"
        )

    # Use provided names or fall back to file stems
    if not names:
        names = [input_file.stem for input_file in inputs]

    all_dfs = []
    for input_file, group_name in zip(inputs, names):
        df = pd.read_csv(input_file, sep="\t")

        # Keep only required columns and ignore others
        required_columns = ["sample", "status", "count", "percentage"]
        df = df[required_columns]

        # Keep only rows with consistent or inconsistent status
        df = df[df["status"].isin(["consistent", "inconsistent"])]

        # Handle percentage conversion - it might be string with "%" or already float
        if df["percentage"].dtype == object:
            df["percentage"] = df["percentage"].str.rstrip("%").astype(float)
        else:
            df["percentage"] = df["percentage"].astype(float)

        # Add file identifier using custom name
        df["file_group"] = group_name
        all_dfs.append(df)

    # Combine all dataframes
    combined_df = pd.concat(all_dfs, ignore_index=True)

    # Make file_group categorical with explicit ordering to preserve input order
    combined_df["file_group"] = pd.Categorical(
        combined_df["file_group"], categories=names, ordered=True
    )

    fig = create_plot(combined_df, title=title)

    # Set DPI to 300 for PNG files
    if output.suffix == ".png":
        fig.write_image(output, scale=3)  # scale=3 gives 300 DPI (default is 100 DPI)
    else:
        fig.write_image(output)

    return combined_df


def run(args, **kwargs):
    combined_df = plot_mendelian(
        inputs=args.input, names=args.names, output=args.output, title=args.title
    )

    # Save table if requested
    if args.table:
        table = create_table(combined_df)
        table.to_csv(args.table, sep="\t", index=False)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="creates a plot of the Mendelian measures from a tsv generated with scripts.mendelian_consistency"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        nargs="+",
        required=True,
        help="Path(s) to the tsv result(s) from scripts.mendelian_consistency.",
    )
    parser.add_argument(
        "-n",
        "--names",
        type=str,
        nargs="+",
        help="Custom names for each input file group. Must match the number of input files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output png svg file.",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Mendelian Violations",
        help="Title of the plot. Default: 'Mendelian Violations'",
    )
    parser.add_argument(
        "--table",
        type=Path,
        help="If provided, will output a combined tsv table of the plotted numbers in tsv format.",
    )
    return parser


def main() -> None:
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
