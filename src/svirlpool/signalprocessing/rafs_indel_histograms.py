# %%

import argparse
from pathlib import Path

import plotly.graph_objects as go

from ..util import util

# %%


def get_indel_counts(input: Path | str) -> tuple[dict[str, object], int, int]:
    counts = {"deletions": {}, "insertions": {}, "breakends": 0}
    n_rafs = 0
    sum_bp = 0
    for raf in util.yield_from_raf(input=input):
        n_rafs += 1
        sum_bp += raf.reference_alignment_end - raf.reference_alignment_start
        for sv in raf.SV_signals:
            if sv.sv_type == 0:  # insertion
                if sv.size not in counts["insertions"]:
                    counts["insertions"][sv.size] = 0
                counts["insertions"][sv.size] += 1
            elif sv.sv_type <= 2:  # deletion
                if sv.size not in counts["deletions"]:
                    counts["deletions"][sv.size] = 0
                counts["deletions"][sv.size] += 1
            else:
                counts["breakends"] += 1
    return counts, n_rafs, sum_bp


def histograms_indel_sizes_per_chr(input: Path, figpath: Path) -> None:
    counts, n_rafs, sum_bp = get_indel_counts(input)
    # now we have a dict of the form {"deletions":{size:count,...},"insertions":{size:count,...}}
    # from that we can create a bar plot, where the x coordinate of a bar is the sv size and the y coordinate is the count
    # for both deletions and insertions
    # we only want to plot the range of sizes [1-100]
    fig = go.Figure()
    for sv_type in ["deletions", "insertions"]:
        sizes = []
        counts_ = []
        for size in counts[sv_type]:
            if size > 100:
                continue
            sizes.append(size)
            counts_.append(counts[sv_type][size] / (sum_bp / 1e6))
        fig.add_trace(go.Bar(x=sizes, y=counts_, name=sv_type))
    fig.update_layout(
        barmode="group",
        title=f"Indel sizes frequencies per 1mb (n rafs={n_rafs} sum bp={sum_bp})",
    )

    if str(figpath).endswith(".html"):
        fig.write_html(figpath)
    elif str(figpath).endswith(".png"):
        fig.write_image(figpath)
    else:
        raise ValueError(
            f"Unsupported file format for figure: {figpath}. Supported formats are .html and .png."
        )


def run(args, **kwargs):
    histograms_indel_sizes_per_chr(input=args.input, figpath=args.figure)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Creates histograms of indel sizes for each sv type ins and del."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="tsv.gz file of read alignment fragments (rafs).",
    )
    parser.add_argument(
        "-f",
        "--figure",
        type=Path,
        required=True,
        help="Path to save the figure. Supported formats are .html and .png.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
