# interval 22 35137826 35142763
# %%
import typing
from pathlib import Path

import pandas as pd
import pysam

# %%
# pseudo-convolve over signals
# save for each CR reads that contribute to ins,del,breakends
# e.g dels:[read0,read1,..]
# save for each CR the connection via deletion to the last CR
# find all connections via breakends
#
# what if a del is partially only available via breakends?
# potential equality check?
#
# %%
# 1) calc all read positions for substrings on table of regions
#   and alignments -> save corellation to CR
# 2) iterate reads and for each read check necessary intervals
#   write subsequences to fastq file with name of CR
# 3) iterate fastq files of CR-related sub fastqs
#   create local assemblies


def get_interval_on_read_in_region(
    a: pysam.libcalignedsegment.AlignedSegment, start: int, end: int
) -> typing.Tuple[int, int]:
    istart = -0
    iend = -1
    positions = a.get_reference_positions(full_length=True)
    for i, p in enumerate(positions):
        if p is None:
            continue
        if p > start:
            istart = max(0, i - 1)
            break
    for j in range(i, len(positions)):
        if positions[j] is None:
            continue
        if positions[j] > end:
            iend = j
            break
    if iend == -1:
        iend = len(positions)
    return istart, iend


def cut_sequence_of_alignmentSegment(
    a: pysam.libcalignedsegment.AlignedSegment, start: int, end: int
) -> str:
    s, e = get_interval_on_read_in_region(a, start, end)
    print(len(a.query_sequence))
    return a.query_sequence[s:e]


# %%
# read bed file as pandas df
path_bed = Path(
    "/home/vinzenz/development/LRSV-detection/development/test/ont-r10/QC/crs.bed"
)
df_bed = pd.read_csv(path_bed, sep="\t", header=None)
output_base_path = Path(
    "/home/vinzenz/development/LRSV-detection/development/test/ont-r10/fastq"
)
# create directory if it does not exist
# somehow buggy: output_base_path.mkdir(parents=True,exist_ok=True)
bamfile = pysam.AlignmentFile(
    "/home/vinzenz/development/LRSV-detection/development/test/ont-r10/minimap2.ont-r10.10x.chr22.bam",
    "rb",
)

# %%
for i, row in df_bed.iterrows():
    gregion = [str(row[0]), int(row[1]), int(row[2])]
    print(gregion)
    region_name = f"{str(gregion[0])}:{str(gregion[1])}-{str(gregion[2])}"
    fname = output_base_path / (region_name + ".fasta")
    # print to fasta, maybe later adjust to fastq
    # which needs to search the read sequence in the original fastq file
    with open(fname, "w") as f:
        alns = list(bamfile.fetch(region=region_name))
        for a in alns:
            print(a.query_name)
            print(">" + a.query_name, file=f)
            print(cut_sequence_of_alignmentSegment(a, row[1], row[2]), file=f)

# read bed file
# for each region in bed file, fetch reads and
# %%
