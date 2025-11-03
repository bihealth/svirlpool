# %%
# mendelian consistency
# this script iterates a VCF file and checks for mendelian consistency among the provided samples.

import argparse
from pathlib import Path

import pandas as pd
import vcfpy

# %%


def inheritance_probabilities():
    # this function gets the genotypes of each variant in the child, father and mother
    # and calculates the probability of each inheritance pattern.
    # The genotypes are coded as 0,1,2 which is the sum of ref=0 and alt=1 alleles per sample genotype.
    # if the child is 0, then the parents can be 0,0 or 0,1 or 1,0 or 1,1.
    # if the child is 1, then the parents can be 0,1 or 1,0 or 1,1.
    # if the child is 2, then the parents can be 1,1 or 1,2 or 2,1 or 2,2.
    # given 0,1,2 for the child, the probabilities of the parent GTs are:
    gt_probs = {
        0: {(0, 0): 0.25, (0, 1): 0.50, (1, 1): 0.25},
        1: {(0, 0): 0.0, (0, 1): 0.5, (1, 1): 0.5},
        2: {(0, 0): 0.0, (1, 1): 0.25, (1, 2): 0.50, (2, 2): 0.25},
    }


def mendelian_errors(
    input: Path,
    father: str = "Father",
    mother: str = "Mother",
    child: str = "Child",
    debug: bool = False,
) -> pd.DataFrame:

    reader = vcfpy.Reader.from_path(input)

    # check in the vcf header if child, father, mother are in the samples
    for sample in [child, father, mother]:
        if sample not in reader.header.samples.names:
            raise ValueError(
                f"Sample {sample} not found in the VCF header. Samples found: {reader.header.samples.names}"
            )

    GT_dict = {}
    for divider in ["/", "|"]:
        for n, num in [(0, "0"), (0, "."), (1, "1"), (1, "2")]:
            for m, num2 in [(0, "0"), (0, "."), (1, "1"), (1, "2")]:
                s = f"{num}{divider}{num2}"
                GT_dict[s] = n + m

    data = []
    for record in reader:
        GTs = {}
        for call in record.calls:
            gt = GT_dict[call.data["GT"]]
            GTs[call.sample] = gt
        # don't save the GTs if all alleles are 0
        if all([GTs[father] == 0, GTs[mother] == 0, GTs[child] == 0]):
            continue
        data.append(GTs)

    df = pd.DataFrame(data)

    error_masks = [
        ((df[father] == 2) | (df[mother] == 2)) & (df[child] == 0),
        ((df[father] == 0) & (df[mother] == 0)) & (df[child] != 0),
        ((df[father] == 0) | (df[mother] == 0)) & (df[child] == 2),
    ]

    total_alleles = df.shape[0]
    total_errors = error_masks[0].sum() + error_masks[1].sum() + error_masks[2].sum()

    print(f"Total alleles: {total_alleles}; total errors: {total_errors}")
    print(f"Error rate: {total_errors/total_alleles:.2f}")

    if debug:
        # print all rows with errors
        for row in df[error_masks[0] | error_masks[1] | error_masks[2]].iterrows():
            print(row)
    return df


# %%


def run(args, **kwargs):
    df = mendelian_errors(
        input=args.input,
        father=args.father,
        mother=args.mother,
        child=args.child,
        debug=args.debug,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Check mendelian consistency in a VCF file and print the basic stats."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="a VCF file with trio information",
    )
    parser.add_argument(
        "-f",
        "--father",
        type=str,
        required=False,
        default="Father",
        help="name of the father sample",
    )
    parser.add_argument(
        "-m",
        "--mother",
        type=str,
        required=False,
        default="Mother",
        help="name of the mother sample",
    )
    parser.add_argument(
        "-c",
        "--child",
        type=str,
        required=False,
        default="Child",
        help="name of the child sample",
    )
    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="print all rows with Mendelian errors",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
