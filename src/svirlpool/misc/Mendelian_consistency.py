# %%
# script to test Mendelian consistency in a vcf file
# testing for exact matches

import argparse
from enum import Enum
from pathlib import PosixPath

import vcfpy

# %%


# parse chr,start and GTs
def parse_variant(variant):
    return (variant.CHROM, variant.POS), [call.data["GT"] for call in variant.calls]


# %%


def possible_inheritances() -> dict[tuple[int, int], set[tuple[int, int]]]:
    p = [(0, 0), (0, 1), (1, 1)]
    all_possible_inheritances = dict[tuple[int, int], set[tuple[int, int]]]()
    for father in p:
        for mother in p:
            # generate all possible genotypes of the son
            all_possible_inheritances[(father, mother)] = set(
                [
                    tuple(sorted([father[a], mother[b]]))
                    for a in range(2)
                    for b in range(2)
                ]
            )
    return all_possible_inheritances


# define an enum for gt inheritance status: consistent, inconsistent, missing
class GTInheritanceStatus(Enum):
    missing = 0
    consistent = 1
    incomplete = 2
    inconsistent = 3


# %%
def is_variant_inconsistent(
    variant: vcfpy.Record, dict_possible_inheritances: dict
) -> GTInheritanceStatus:
    # get the genotypes
    gt_son, gt_father, gt_mother = [call.data["GT"] for call in variant.calls]
    if gt_son == "./.":
        return GTInheritanceStatus.missing
    elif gt_father == "./." and gt_mother == "./.":
        return GTInheritanceStatus.missing
    if "./." not in [gt_son, gt_father, gt_mother]:
        # exact matching
        # parse the genotypes to tuples
        gt_son = tuple(map(int, gt_son.split("/")))
        gt_father = tuple(map(int, gt_father.split("/")))
        gt_mother = tuple(map(int, gt_mother.split("/")))
        # if all genotypes are one of ["0/0", "0/1", "1/1"], then search for exact matches
        if gt_son in dict_possible_inheritances[(gt_father, gt_mother)]:
            return GTInheritanceStatus.consistent
    else:
        if gt_father == "./.":
            gt_son = tuple(map(int, gt_son.split("/")))
            gt_mother = tuple(map(int, gt_mother.split("/")))
            # try all three possibilities where the father has any genotype
            for father in [(0, 0), (0, 1), (1, 1)]:
                if gt_son in dict_possible_inheritances[(father, gt_mother)]:
                    return GTInheritanceStatus.incomplete
        if gt_mother == "./.":  # redundant, but for better understanding
            gt_son = tuple(map(int, gt_son.split("/")))
            gt_father = tuple(map(int, gt_father.split("/")))
            # try all three possibilities where the mother has any genotype
            for mother in [(0, 0), (0, 1), (1, 1)]:
                if gt_son in dict_possible_inheritances[(gt_father, mother)]:
                    return GTInheritanceStatus.incomplete
    return GTInheritanceStatus.inconsistent


# %%


def variants_consistency_stats(path_vcf: PosixPath, output: PosixPath):
    # load the vcf file
    reader = vcfpy.Reader.from_path(path_vcf)
    # load all variants
    variants = [record for record in reader]
    # get all possible inheritances
    dict_possible_inheritances = possible_inheritances()
    # get the number of variants
    n_variants = len(variants)
    # get the number of inconsistent variants
    # counter for each GTInheritanceStatus
    with open(output, "w") as f:
        dict_stats = {
            GTInheritanceStatus.missing: 0,
            GTInheritanceStatus.consistent: 0,
            GTInheritanceStatus.incomplete: 0,
            GTInheritanceStatus.inconsistent: 0,
        }
        for variant in variants:
            status = is_variant_inconsistent(variant, dict_possible_inheritances)
            dict_stats[status] += 1
            # if inconsistent, print the variant CHROM and POS
            print(f"{variant.CHROM}:{variant.POS}\t{status.name}", file=f)
    # print results to terminal
    print(f"Number of variants: {n_variants}")
    for status in dict_stats:
        print(f"{status.name}: {dict_stats[status]}")


def run(args, **kwargs):
    variants_consistency_stats(args.input, args.output)
    return


def get_parser():
    parser = argparse.ArgumentParser(
        description="Build a benchmark for multi-sample SV calling. Output is a json file containing a dict"
    )
    parser.add_argument(
        "-i", "--input", type=PosixPath, required=True, help="Path to vcf(.gz) file"
    )
    parser.add_argument(
        "-o", "--output", type=PosixPath, required=True, help="Path to output tsv file"
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
# %%
