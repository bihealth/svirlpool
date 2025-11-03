## this script sorts a bed file according to a given reference order
# %%
import argparse
import typing
from pathlib import Path


# %%
def create_ref_dict(reference: Path) -> typing.Dict[int, str]:
    return {
        i: str(line.rstrip().split("\t")[0])
        for i, line in enumerate(open(Path(str(reference) + ".fai"), "r"))
    }


def run(args, **kwargs):
    refdict = create_ref_dict(args.reference)
    refdict_inversed = {v: k for k, v in refdict.items()}
    with open(args.input) as f:
        data = [line.rstrip().split("\t") for line in f]
    data = sorted(data, key=lambda x: [refdict_inversed[x[0]], int(x[1]), int(x[2])])
    if args.output:
        with open(args.output, "w") as f:
            for line in data:
                print("\t".join(line), file=f)
    else:  # copy to terminal
        for line in data:
            print("\t".join(line))


def get_parser():
    parser = argparse.ArgumentParser(
        description="sorts a bed file according to the order of chromosomes that are given in a reference fasta file."
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to reference fasta file. The index of the file (.fai) must be present in the same directory.",
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="Path to input bed file."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=False,
        help="Path to output file. Prints to command line if not given.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
