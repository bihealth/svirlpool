# extract pos and svlen from zipped vcf file and write to bed file

# %%
import argparse
from pathlib import Path

import vcfpy


def vcf_deletions_to_bed(path_vcf: Path, path_bed: Path, min_size: int):
    reader = vcfpy.Reader.from_path(path_vcf)
    with open(path_bed, "w") as f:
        for record in reader:
            # add deletions (DEL)
            if "SVTYPE" in record.INFO and record.INFO["SVTYPE"] == "DEL":
                svlen = abs(record.INFO["SVLEN"][0])
                if svlen >= min_size:
                    print(
                        f"{record.CHROM}\t{record.POS}\t{record.POS+svlen}\t{svlen}",
                        file=f,
                    )
            # add insertions (INS)
            if "SVTYPE" in record.INFO and record.INFO["SVTYPE"] == "INS":
                svlen = abs(record.INFO["SVLEN"][0])
                if svlen >= min_size:
                    print(
                        f"{record.CHROM}\t{record.POS}\t{record.POS+1}\t{svlen}", file=f
                    )


# vcf_deletions_to_bed(path_vcf,path_bed,50)
# %%


def run(args, **kwargs):
    vcf_deletions_to_bed(
        path_vcf=args.input, path_bed=args.output, min_size=args.min_size
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Extract insertions and deletions from a VCF file and write to a BED file."
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="(zipped) VCF file"
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Path to BED file."
    )
    parser.add_argument(
        "--min_size",
        type=int,
        default=20,
        help="Minimum size of SVs to be included in BED file.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
