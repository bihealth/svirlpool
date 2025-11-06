# %%

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
from tqdm import tqdm

from . import datatypes, util


# %%
def get_n_indels_signals(
    raf: datatypes.ReadAlignmentFragment, min_size: int, max_size: int
) -> int:
    return sum([
        1
        for sv in raf.SV_signals
        if (
            min_size <= sv.size <= max_size
            and sv.sv_type < 2
            and (
                raf.effective_interval[1] <= sv.ref_start <= raf.effective_interval[2]
                or raf.effective_interval[1] <= sv.ref_end <= raf.effective_interval[2]
            )
        )
    ])


def get_summed_indels_sizes(
    raf: datatypes.ReadAlignmentFragment, min_size: int, max_size: int
) -> int:
    return sum([
        sv.size
        for sv in raf.SV_signals
        if (
            min_size <= sv.size <= max_size
            and sv.sv_type < 2
            and (
                raf.effective_interval[1] <= sv.ref_start <= raf.effective_interval[2]
                or raf.effective_interval[1] <= sv.ref_end <= raf.effective_interval[2]
            )
        )
    ])


def get_median_indel_size(
    raf: datatypes.ReadAlignmentFragment, min_size: int, max_size: int
) -> int:
    sizes = np.array([
        sv.size
        for sv in raf.SV_signals
        if (
            min_size <= sv.size <= max_size
            and sv.sv_type < 2
            and (
                raf.effective_interval[1] <= sv.ref_start <= raf.effective_interval[2]
                or raf.effective_interval[1] <= sv.ref_end <= raf.effective_interval[2]
            )
        )
    ])
    if len(sizes) == 0:
        return 0
    return np.median(sizes)


# %%
def parse_rafs_to_df(path_rafs) -> pd.DataFrame:
    columns = [
        "readname",
        "effective_size",
        "n_indels",
        "sum_indel_bp",
        "median_indel_size",
    ]
    L = []
    for raf in tqdm(util.yield_from_raf(path_rafs)):
        L.append([
            raf.read_name,
            raf.effective_interval[2] - raf.effective_interval[1],
            get_n_indels_signals(raf, 1, 50),
            get_summed_indels_sizes(raf, 1, 50),
            get_median_indel_size(raf, 1, 50),
        ])
    return pd.DataFrame(L, columns=columns)


def plot_rafs_df_QC_histograms(df: pd.DataFrame, path_QC: Path) -> None:
    df["avg_indel_size"] = df["sum_indel_bp"] / df["n_indels"]
    df["n_indels_per_1kbp"] = df["n_indels"] / df["effective_size"] * 1000
    # plot histograms of n_indels_per_1kbp of df with n_indels > 3 and effective_size >= 10000
    for col in df.columns:
        if col == "readname":
            continue
        fig = px.histogram(df.query("effective_size >= 10000").loc[:, [col]], x=col)
        fig.write_html(path_QC / f"QC_rafs.{col}.html")


def QC_rafs(path_rafs: Path, path_QC: Path) -> None:
    df = parse_rafs_to_df(path_rafs)
    plot_rafs_df_QC_histograms(df, path_QC)


def run(args, **kwargs):
    QC_rafs(path_rafs=args.input, path_QC=args.output)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Convolves over the singal to get the localized signalstrength to be used to create candidate regions."
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="tsv file of candidate regions."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help=r"output directory for the plots of name scheme 'QC_crs.\{colname\}.html.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
