# this script gte a bed file of mononucleotide stretches of a given length on a given reference genome
# and intersects it with a signals bed-like file.
# stats are then computed to determine if a sequencing sample is impacted by mononucleotide stretches and if so,
# the signal overlapping the mononucleotides is filtered.

# %%

import argparse
import subprocess
import tempfile
from pathlib import Path
from shlex import split

import pandas as pd
import plotly.express as px
from logzero import logging as log
from tqdm import tqdm

from . import util

# #%%

# # re-build with new signals (contain only INS,DEL,BNDL,BNDR)
# # use bedtools subtract correctly
# # filter any signal that is very close to the annotated mononucleotide stretches
# # and has a similar size to the mononucleotide stretch or is smaller than the mononucleotide stretch

# def filter_mononucleotide_exclusive_deletions(
#         mononucleotides:Path,
#         signals:Path,
#         reference:Path,
#         output:Path,
#         margin:int,
#         overlap_mononucleotide_stretch:float,
#         overlap_signal:float) -> None:
#     """overlap_signal would need to be big (e.g. 0.8) and overlap_mononucleotide_stretch small (e.g. 0.2)."""

#     tmp_dir = tempfile.TemporaryDirectory()

#     # add margin to mononucleotides
#     log.info(f"Adding margin of {margin} bp to mononucleotides")
#     path_slop_mononucleotides = Path(tmp_dir.name) / "mononucleotides_slop.bed"
#     cmd_slop = split(f"bedtools slop -b {str(margin)} -i {str(mononucleotides)} -g {str(reference)}.fai")
#     with open(path_slop_mononucleotides, 'w') as f:
#         subprocess.run(cmd_slop, stdout=f)

#     # convert mononucleotides chr to chrID
#     log.info("Converting mononucleotides chr to chrID")
#     path_chrID_mononucleotides = Path(tmp_dir.name) / "mononucleotides_chrID.bed"
#     util.bed_chr_to_chrID(input=mononucleotides,output=path_chrID_mononucleotides,reference=reference)


#     log.info("Subtracting signals that overlap mononucleotide stretches from signals..")
#     # bedtools subtract to remove any signal that is very close (or overlaps) and has overlaps and has similar size or is smaller than the mononucleotide stretch
#     cmd_subtract = split(f"bedtools subtract -A -f {overlap_signal} -F {overlap_mononucleotide_stretch} -a {str(signals)} -b {str(path_slop_mononucleotides)}")
#     with open(output, 'w') as f:
#         subprocess.run(cmd_subtract, stdout=f)
#     log.info("done")
def plot_QC(df: pd.DataFrame, figure: Path) -> None:
    df["stretch_exclusive"] = 1
    mask_nonexclusive = df.query(
        "(svtype == 'DELR' or svtype == 'DELL') and svsize > stretch"
    ).index
    df.loc[mask_nonexclusive, "stretch_exclusive"] = 0
    fig = px.histogram(
        df.query("svtype == 'DELL' or svtype == 'DELR'"),
        x="stretch",
        facet_row="svtype",
        color="stretch_exclusive",
        barmode="stack",
    )
    fig.write_image(figure.with_suffix(".png"))
    fig.write_html(figure.with_suffix(".html"))


# %%
def filter_mononucleotide_exclusive_deletions(
    mononucleotides: Path,
    signals: Path,
    reference: Path,
    output: Path,
    margin: int,
    figure: Path = Path(""),
) -> None:
    # %%
    tmp_dir = tempfile.TemporaryDirectory()

    # add margin of 2 bp to mononucleotides
    log.info(f"Adding margin of {margin} bp to mononucleotides")
    path_slop_mononucleotides = Path(tmp_dir.name) / "mononucleotides_slop.bed"
    cmd_slop = split(
        f"bedtools slop -b {str(margin)} -i {str(mononucleotides)} -g {str(reference)}.fai"
    )
    with open(path_slop_mononucleotides, "w") as f:
        subprocess.run(cmd_slop, stdout=f)

    # convert mononucleotides chr to chrID
    log.info("Converting mononucleotides chr to chrID")
    path_chrID_mononucleotides = Path(tmp_dir.name) / "mononucleotides_chrID.bed"
    util.bed_chr_to_chrID(
        input=mononucleotides, output=path_chrID_mononucleotides, reference=reference
    )

    # intersect mononucleotides with signals so that a selection of signals is left
    log.info("Intersecting mononucleotides with signals")
    path_signals_on_mononucleotides = (
        Path(tmp_dir.name) / "signals_on_mononucleotides.bed"
    )
    cmd_intersect = split(
        f"bedtools intersect -wo -a {str(signals)} -b {str(path_chrID_mononucleotides)}"
    )
    with open(path_signals_on_mononucleotides, "w") as f:
        subprocess.run(cmd_intersect, stdout=f)
    # %%
    # load intersection to np array
    log.info("Loading intersection to dataframe")
    df_signals_on_mononucleotides = pd.read_csv(
        path_signals_on_mononucleotides,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14],
    )
    df_signals_on_mononucleotides.columns = "chrID start end svtype readstart readend svsize readname sampleID chr forward stretch".split()
    # np_signals_on_mononucleotides = np.loadtxt(path_signals_on_mononucleotides, dtype=str)
    # %%
    # put all dell and delr that are to be kept in a set.
    log.info("Putting all dell and delr that are to be kept in a set")
    df_delr = df_signals_on_mononucleotides.query(
        "svtype == 'DELR' and svsize <= stretch*2"
    ).sort_values(by=["chrID", "start"])
    delr_tuples = set(
        df_delr.apply(lambda x: (x["chrID"], x["start"], x["readname"]), axis=1)
    )
    df_dell = df_signals_on_mononucleotides.query(
        "svtype == 'DELL' and svsize <= stretch*2"
    ).sort_values(by=["chrID", "start"])
    dell_tuples = set(
        df_dell.apply(lambda x: (x["chrID"], x["start"], x["readname"]), axis=1)
    )
    # augmenting DFs
    dell_tuples = dell_tuples.union(
        set(
            df_delr.apply(
                lambda x: (x["chrID"], x["start"] - x["svsize"], x["readname"]), axis=1
            )
        )
    )
    delr_tuples = delr_tuples.union(
        set(
            df_dell.apply(
                lambda x: (x["chrID"], x["start"] + x["svsize"], x["readname"]), axis=1
            )
        )
    )
    # %%
    if figure != Path(""):
        log.info("Plotting QC")
        plot_QC(df=df_signals_on_mononucleotides, figure=figure)
    # %%
    log.info(
        "iterating over signals and dropping DEL signals that are exclusive to mononucleotide stretches"
    )
    removed_signal = 0
    with open(output, "w") as f:
        with open(signals, "r") as g:
            for line in tqdm(g):
                line = line.strip().split("\t")
                chrID, start = int(line[0]), int(line[1])
                svtype, readname = line[3], line[7]
                if svtype not in ("DELL", "DELR"):
                    print(*line, sep="\t", file=f)
                    continue
                t = (chrID, start, readname)
                if svtype == "DELR":
                    if t in delr_tuples:
                        # remove tuple from delr_tuples
                        delr_tuples.remove(t)
                        removed_signal += 1
                        continue
                if svtype == "DELL":
                    if t in dell_tuples:
                        # remove tuple from dell_tuples
                        dell_tuples.remove(t)
                        removed_signal += 1
                        continue
                # if deletion is not stretch exclusive, keep it
                print(*line, sep="\t", file=f)

    log.info(f"removed {removed_signal:,} signals")


# %%


def run(args, **kwargs):
    filter_mononucleotide_exclusive_deletions(
        mononucleotides=args.mononucleotides,
        signals=args.signals,
        reference=args.reference,
        output=args.output,
        margin=args.margin,
        figure=args.figure,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Filters signals for deletions that are mostly covered by mononucleotide stretches."
    )
    parser.add_argument(
        "-m",
        "--mononucleotides",
        type=Path,
        required=True,
        help="Path to mononucleotides bed file.",
    )
    parser.add_argument(
        "-s", "--signals", type=Path, required=True, help="Path to signals bed file."
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to reference fasta file.",
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Path to output file."
    )
    parser.add_argument(
        "-g",
        "--margin",
        type=int,
        required=False,
        default=5,
        help="Margin to add to mononucleotides.",
    )
    parser.add_argument(
        "-f",
        "--figure",
        type=Path,
        required=False,
        default=Path(""),
        help="Path to figure file. can be either png or html. However, both html and png are written.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
