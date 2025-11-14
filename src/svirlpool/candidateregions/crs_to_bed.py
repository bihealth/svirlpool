# %%
import argparse
import csv
from pathlib import Path

from ..candidateregions import signalstrength_to_crs

# %%


def print_crs_to_bed(input: Path, output: Path) -> None:
    writer = csv.writer(open(output, "w"), delimiter="\t")
    for cr in signalstrength_to_crs.load_crs_from_db(path_db=input):
        writer.writerow([
            cr.chr,
            cr.referenceStart,
            cr.referenceEnd,
            cr.crID,
            cr.referenceEnd - cr.referenceStart,
        ])
    return None


# argparse block
def run(args, **kwargs):
    print_crs_to_bed(input=args.input, output=args.output)


def get_parser():
    parser = argparse.ArgumentParser(
        description="creates a bed file with columns chr, start, end, crID."
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="Path to crs database file."
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Path to output file."
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
