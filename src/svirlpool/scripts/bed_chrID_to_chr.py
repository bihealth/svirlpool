import argparse
from pathlib import Path

from . import util


def run(args, **kwargs):
    util.bed_chrID_to_chr(
        input=args.input, reference=args.reference, output=args.output
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="replaces the chr names in the first column of a tsv (bed) with their chr IDs and prints to std out."
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="Path to bed file"
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to reference genome (only the index '.fai' is used).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(""),
        required=False,
        help="Path to output. Prints to std out if no output is provided.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
