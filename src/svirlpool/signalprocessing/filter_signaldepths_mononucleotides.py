import argparse
import logging
import subprocess
import tempfile
from pathlib import Path
from shlex import split

from ..signalprocessing import alignments_to_rafs
from ..util import util

log = logging.getLogger(__name__)


def filter_mononucleotide_exclusive_deletions(
    mononucleotides: Path,
    signaldepths: Path,
    reference: Path,
    output: Path,
    margin: int,
    threads: int,
) -> None:
    # add margin of 2 bp to mononucleotides
    log.info(f"Adding margin of {margin} bp to mononucleotides")
    tmp_slop_mononucleotides = tempfile.NamedTemporaryFile(
        delete=True, suffix="slop_mononucleotides.bed"
    )
    cmd_slop = split(
        f"bedtools slop -b {str(margin)} -i {str(mononucleotides)} -g {str(reference)}.fai"
    )
    with open(tmp_slop_mononucleotides.name, "w") as f:
        subprocess.run(cmd_slop, stdout=f)

    # convert mononucleotides chr to chrID
    log.info("Converting mononucleotides chr to chrID")
    tmp_chrID_mononucleotides = tempfile.NamedTemporaryFile(
        delete=True, suffix="chrID_mononucleotides.bed"
    )
    util.bed_chr_to_chrID(
        input=mononucleotides,
        output=Path(tmp_chrID_mononucleotides.name),
        reference=reference,
    )

    # intersect mononucleotides with signals so that a selection of signals is left
    tmp_out = tempfile.NamedTemporaryFile(delete=True, suffix=".out.signaldepths.tsv")
    log.info("Intersecting mononucleotides with signals")
    cmd_intersect = split(
        f"bedtools intersect -v -f 0.5 -r -a {str(signaldepths)} -b {tmp_chrID_mononucleotides.name}"
    )
    with open(tmp_out.name, "w") as f:
        subprocess.run(cmd_intersect, stdout=f)

    alignments_to_rafs.compress_and_index_bedlike(
        sort_numerically=True, input=Path(tmp_out.name), output=output, threads=threads
    )


# %%


def run(args, **kwargs):
    filter_mononucleotide_exclusive_deletions(
        mononucleotides=args.mononucleotides,
        signaldepths=args.signaldepths,
        reference=args.reference,
        output=args.output,
        margin=args.margin,
        threads=args.threads,
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
        "-s",
        "--signaldepths",
        type=Path,
        required=True,
        help="Path to bgzipped signaldepths tsv file.",
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
        "-t",
        "--threads",
        type=int,
        required=False,
        default=1,
        help="Number of threads to use for compression.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
