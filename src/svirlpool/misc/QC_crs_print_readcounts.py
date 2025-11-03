import argparse
from pathlib import Path

from scripts.util import load_crs


def print_crs_read_counts(path_crs: Path):
    print("crID all_reads unique_reads".split(" "), sep="\t")
    for cr in load_crs(path_crs):
        fragments_readnames = [svs[7] for svs in cr.sv_signals]
        readnames = set(fragments_readnames)
        print(cr.crID, len(fragments_readnames), len(readnames), sep="\t")


def run(args, **kwargs):
    print_crs_read_counts(path_crs=args.input)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Prints the read counts for each cr to the terminal."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the fastq file of cut reads.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
