### create non-unique regions file from mappability track
# %%
import argparse
import subprocess
import tempfile
from pathlib import Path
from shlex import split

# %%
from . import util

# %%

# path_mappability = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/wgEncodeCrgMapabilityAlign100mer.bedgraph")
# path_reference = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/hs37d5.fa")
# output = Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/wgEncodeCrgMapabilityAlign50mer.lt_1.slop1000.merged50.bed")


def check_if_all_chromosomes_are_in_refdict(
    path_mappability: Path, path_reference: Path
):
    refdict = util.create_ref_dict(reference=path_reference)
    cmd_cut = f"cut -f1 {path_mappability}"
    cmd_sort = "sort"
    cmd_uniq = "uniq"
    # write output to variable chromosomes and check if all chromosomes are in refdict.values()
    p0 = subprocess.Popen(split(cmd_cut), stdout=subprocess.PIPE)
    p1 = subprocess.Popen(split(cmd_sort), stdin=p0.stdout, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(split(cmd_uniq), stdin=p1.stdout, stdout=subprocess.PIPE)
    chromosomes = p2.communicate()[0].decode().strip().split("\n")
    chromosomes_in_mappability: set = set(chromosomes)
    chromosomes_in_reference: set = set(refdict.values())
    if not chromosomes_in_mappability.issubset(chromosomes_in_reference):
        raise ValueError(
            f"The following chromosomes in the mappability track {str(path_mappability)} are not a subset of the reference chromosomes: {chromosomes_in_mappability-chromosomes_in_reference}. Make sure that all chromosome names in the mappability track are in the reference genome."
        )


def split_large_regions_into_flanking_regions(
    unsplit_regions: Path, output: Path, max_region_size: int, flanking_region_size: int
):
    # the problem arises, that some regions are just too big (greater than 30kb). In that case, they are split
    # into two much smaller (5kb) regions flanking the original region. Any region that is smaller than 30 kb is
    # kept as it is.
    with open(output, "w") as out:
        with open(unsplit_regions) as f:
            for line in f:
                chrom, start, end = line.strip().split("\t")
                start = int(start)
                end = int(end)
                if end - start > max_region_size:
                    end_left = start + flanking_region_size
                    end_right = end - flanking_region_size
                    out.write(f"{chrom}\t{start}\t{end_left}\n")
                    out.write(f"{chrom}\t{end_right}\t{end}\n")
                else:
                    out.write(f"{chrom}\t{start}\t{end}\n")


def create_nonunique_regions(
    path_mappability: Path,
    path_reference: Path,
    k: int,
    output: Path,
    max_region_size: int,
    flanking_region_size: int,
):
    if 2 * flanking_region_size > max_region_size:
        raise (
            f"Flanking region size {flanking_region_size} is too big for the maximum region size {max_region_size}. It can not be bigger than half of the maximum region size."
        )
    if not 30 < int(k) < 300:
        raise ValueError(f"K-mer size {k} is not in the range [30,300].")
    # check if all chr names are in the refdict
    check_if_all_chromosomes_are_in_refdict(
        path_mappability=path_mappability, path_reference=path_reference
    )

    tmp_reference_bed = tempfile.NamedTemporaryFile(
        delete=True, suffix=".reference.bed"
    )
    util.fai_to_bed(reference=path_reference, output=tmp_reference_bed.name)

    cmd_awk = f"awk -v OFS='\\t' '$4==1 {{print $1,$2+{int(k/2)},$3+{int(k/2)}}}' {path_mappability}"
    cmd_subtract = f"bedtools subtract -a {tmp_reference_bed.name} -b -"
    tmp_unsplit_regions = tempfile.NamedTemporaryFile(
        delete=False, suffix=".unsplit.bed"
    )
    with open(tmp_unsplit_regions.name, "w") as f:
        p0 = subprocess.Popen(split(cmd_awk), stdout=subprocess.PIPE)
        p1 = subprocess.Popen(split(cmd_subtract), stdin=p0.stdout, stdout=f)
        p1.communicate()

    # post-processing: split large regions into flanking regions
    split_large_regions_into_flanking_regions(
        unsplit_regions=tmp_unsplit_regions.name,
        output=output,
        max_region_size=max_region_size,
        flanking_region_size=flanking_region_size,
    )
    # the regions should still be in order.


def run(args, **kwargs):
    create_nonunique_regions(
        path_mappability=args.input,
        path_reference=args.reference,
        k=args.k,
        output=args.output,
        max_region_size=args.max_region_size,
        flanking_region_size=args.flanking_region_size,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Creates a file of non-unique regions on a given reference genome and a mappability track (regions of less than 1 mappability) and adds a merged margin to ensure the sufficient size of candidate regions."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to mappability track in bedgraph format with the columns chr,start,end,mappability. Mappability is a value in the range [0,1].",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to reference fasta file.",
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Path to output bed file."
    )
    parser.add_argument(
        "-k", type=int, required=True, help="K-mer size of the mappability track."
    )
    parser.add_argument(
        "--max_region_size",
        type=int,
        default=30000,
        help="Maximum size of a region in the output file. If a region is larger than this, it will be split into two flanking regions.",
    )
    parser.add_argument(
        "--flanking_region_size",
        type=int,
        default=5000,
        help="Size of the flanking regions that are created when a region is split.",
    )

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
