import argparse
from pathlib import Path

from . import util


def run(args, **kwargs):
    util.generate_sampledicts(alignments=args.alignments, output=args.output)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Creates dictionaries to handle multiple bam and fastq files from multiple sequencing samples."
    )
    parser.add_argument(
        "-i",
        "--alignments",
        type=Path,
        required=True,
        nargs="+",
        help="Path to all used bam files. The order must be the same as in reads.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to dictionaries in json format.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
