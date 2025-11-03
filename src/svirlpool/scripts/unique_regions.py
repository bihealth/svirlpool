# %%
# this script computes a unique regions file from:
# a reference.fasta (load .fai and generate bedtools .genome file frm it)
# a reference.bed is created from the .fai file.
# a set of bed files of non-unique contributers (e.g. SegDups, lowMappability, tandem-repeats)
# the bed files are merged and slopped by buffer_radius bp.
# the final file is created with bedtools subtract -a reference.bed -b merged_slopped.bed
# %%
import argparse
import subprocess
import tempfile
from pathlib import Path
from shlex import split

from logzero import logger as log

from . import util

# %%


def check_chr_are_in_reference(bedfile: Path, reference: Path) -> set:
    reference_fai = util.create_fai_if_not_exists(reference=reference)
    with open(reference_fai) as f:
        reference_chrs = [line.split("\t")[0] for line in f]
    with open(bedfile) as f:
        bed_chrs = set([line.split("\t")[0] for line in f])
    # return all chromosome names in bedfile that are not in reference
    return bed_chrs - set(reference_chrs)


def merge_bedfiles_and_subtract_from_reference(
    paths_bedfiles: list[Path], output: Path, ref_bed: Path
) -> None:
    # merge bed files and subtract
    cmd_cat = f"cat {' '.join([str(p) for p in paths_bedfiles])}"
    cmd_cut = "cut -f 1-3 "
    cmd_sort = "bedtools sort -i - "
    cmd_merge_beds = "bedtools merge -i - "
    cmd_subtract = f"bedtools subtract -a {ref_bed} -b - "
    cmd_merge_final = "bedtools merge -d 100 -i - "
    with open(output, "w") as f:
        p0 = subprocess.Popen(split(cmd_cat), stdout=subprocess.PIPE)
        p1 = subprocess.Popen(split(cmd_cut), stdin=p0.stdout, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(split(cmd_sort), stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(
            split(cmd_merge_beds), stdin=p2.stdout, stdout=subprocess.PIPE
        )
        p4 = subprocess.Popen(
            split(cmd_subtract), stdin=p3.stdout, stdout=subprocess.PIPE
        )
        p5 = subprocess.Popen(split(cmd_merge_final), stdin=p4.stdout, stdout=f)
        p5.communicate()


# %%
# input
# output = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/hs37d5.unique.bed")
# reference = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/hs37d5.fa")
# paths_bedfiles = [
#     Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/wgEncodeCrgMapabilityAlign100mer.centered.bedgraph"),
#     Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/segdups.merged.bed"),
#     Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/human_hs37d5.trf.slop50.bed"),
#     Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/hs37d5.Ns.bed")
# ]


def create_unique_regions_bed_file(
    paths_bedfiles: list[Path], reference: Path, output: Path, min_size: int
) -> None:
    tmp_ref_bed = tempfile.NamedTemporaryFile(delete=True, suffix=".reference.bed")
    util.fai_to_bed(reference=reference, output=tmp_ref_bed.name)
    # check all chromosome names in all bed files if they are a subset of the reference chr names
    log.info(
        "Checking if all chromosome names in bed files are a subset of the reference.."
    )
    for bedfile in paths_bedfiles:
        exclusive_chrs = check_chr_are_in_reference(
            bedfile=bedfile, reference=reference
        )
        if exclusive_chrs:
            print(
                f"Chromosome names {exclusive_chrs} in {bedfile} are not a subset of the reference."
            )
            raise ValueError(
                f"Chromosome names {exclusive_chrs} in {bedfile} are not a subset of the reference."
            )
    log.info("Merge bed files and subtract from reference to create the output..")
    tmp_unique = tempfile.NamedTemporaryFile(delete=False, suffix=".unique.bed")
    merge_bedfiles_and_subtract_from_reference(
        paths_bedfiles=paths_bedfiles,
        output=Path(tmp_unique.name),
        ref_bed=Path(tmp_ref_bed.name),
    )
    # filter out regions smaller than min_size
    cmd_filter = f"awk '{{if ($3-$2 > {min_size}) print $0}}' {tmp_unique.name}"
    with open(output, "w") as f:
        p0 = subprocess.Popen(split(cmd_filter), stdout=f)
        p0.communicate()


# %%
# subtract merged bed from reference bed


def run(args, **kwargs):
    create_unique_regions_bed_file(
        paths_bedfiles=args.bedfiles,
        reference=args.reference,
        output=args.output,
        min_size=args.min_size,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Create a bed file of unique regions from a reference and a set of bed files of non-unique regions."
    )
    parser.add_argument(
        "-b", "--bedfiles", nargs="+", type=Path, help="path to input bed files"
    )
    parser.add_argument(
        "-r", "--reference", type=Path, help="path to reference fasta file"
    )
    parser.add_argument("-o", "--output", type=Path, help="path to output bed file")
    parser.add_argument(
        "--min-size",
        type=int,
        default=500,
        help="minimum size of unique regions. Smaller regions are discarded.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
