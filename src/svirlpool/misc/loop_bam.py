# %%
# this script lopps through a bam file, alignment by alignment, and exits. It is to measure the time required

import argparse
from pathlib import Path

import pysam


def loop_bam(path: Path) -> None:
    bamfile = pysam.AlignmentFile(path, "rb")
    for read in bamfile:
        pass
    bamfile.close()
    return


def run(args):
    path = Path(args.input)
    loop_bam(path)
    return


def get_parser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="Input bam file", required=True)
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
