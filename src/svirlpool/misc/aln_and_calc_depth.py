import argparse
import shlex
import subprocess
from pathlib import Path

import pandas as pd

from . import util


def align_and_calc_depth(args, **kwargs):
    # align
    bamout = args.alndir / f"{args.tech}.unknown.alignments.bam"
    util.align_reads_with_minimap(
        args.reference, args.reads, bamout, args.tech_minimap, 32
    )
    # calc depth
    cmd_depth = f"mosdepth {args.tech}.max -t 32 {str(bamout)}"
    subprocess.check_call(shlex.split(cmd_depth))
    avg_depth = float(
        pd.read_csv(f"mosdepth {args.tech}.max.mosdepth.summary.txt", sep="\t").iloc[
            -1, 3
        ]
    )
    # df = execute_to_df(shlex.split(cmd_depth))
    # avg_depth = float(df.iloc[-1,3])
    # rename alignment
    new_bamout = args.alndir / f"{args.tech}.{str(round(avg_depth))}x.alignments.bam"
    cmd_mv = f"mv {bamout} {new_bamout}"
    subprocess.check_call(cmd_mv)
    cmd_mvb = f"mv {bamout}.bai {new_bamout}.bai"
    subprocess.check_call(cmd_mvb)
    print(str(round(avg_depth)))


def get_parser():
    parser = argparse.ArgumentParser(
        description="re-calculates the performance statistics of truvari output"
    )
    parser.add_argument(
        "--reference", type=Path, required=True, help="Path to reference genome"
    )
    parser.add_argument(
        "--reads",
        type=Path,
        required=True,
        help="Path to reads file [fasta,fa,fastq,fq]",
    )
    parser.add_argument(
        "--alndir",
        type=Path,
        required=True,
        help="Path to directory where the alignments are saved",
    )
    parser.add_argument(
        "--tech", type=str, required=True, help="Tech, e.g. ont-r9 or pb-hifi"
    )
    parser.add_argument(
        "--tech_minimap",
        type=str,
        required=True,
        help="Tech to configure minimap2 aligner, e.g. ont or pb-hifi.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    align_and_calc_depth(args)
    return


if __name__ == "__main__":
    main()
