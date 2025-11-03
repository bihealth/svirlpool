# %%
# iterate a reference genome and find strecthes of adjacent Ns.
# create a bed file with the coordinates of these stretches

import multiprocessing as mp
import typing
from pathlib import Path

from Bio import SeqIO
from tqdm import tqdm

# %%

reference = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/hs37d5.fa"
)
output = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/hs37d5.Ns.bed"
)

# %%


def process_chr(record) -> typing.List[typing.Tuple[str, int, int]]:
    N_regions: typing.List[typing.Tuple[str, int, int]] = []
    N_start = -1
    last_letter = "9"
    for i, letter in enumerate(record.seq):
        if letter == "N" and last_letter != "N":
            N_start = i
        if letter != "N" and last_letter == "N":
            N_regions.append((record.id, N_start, i))
        last_letter = letter
    if last_letter == "N":
        N_regions.append((record.id, N_start, len(record.seq)))
    return N_regions


def regions_with_Ns(
    reference: Path, threads: int
) -> typing.List[typing.Tuple[str, int, int]]:
    N_regions = []
    with mp.Pool(threads) as pool:
        records = list(SeqIO.parse(reference, "fasta"))
        for N_regions_chr in tqdm(pool.imap(process_chr, records), total=len(records)):
            N_regions.extend(N_regions_chr)


# %%
regions = regions_with_Ns(reference=reference, threads=2)
#
# write regions to output as bed file
with open(output, "w") as f:
    for region in regions:
        f.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")
