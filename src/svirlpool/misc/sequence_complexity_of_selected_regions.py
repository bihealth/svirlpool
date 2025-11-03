# %%

import subprocess
import typing
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

# %%


def dna_to_base5(letter: str) -> int:
    if letter == "A":
        return 0
    if letter == "C":
        return 1
    if letter == "G":
        return 2
    if letter == "T":
        return 3
    return 4


def kmer_to_hash(kmer: str, base: int = 3) -> int:
    """Max 20 mers allowed"""
    hashes = np.array([dna_to_base5(letter) for letter in kmer], dtype=int)
    kmer_hash = 0
    for x in hashes:
        kmer_hash = (kmer_hash << base) + x
    return kmer_hash


def sequence_to_kmer_hashes(seq, k: int) -> typing.List[int]:
    return [kmer_to_hash(seq[i : i + k]) for i in range(len(seq) - k + 1)]


def extend_hash(kmer_hash: int, letter_dna5: int, k: int, base: int = 3) -> int:
    return ((kmer_hash << base) + letter_dna5) & 2 ** (k * base) - 1


def sequence_complexity(
    seq, K: typing.List[int] = [3, 5, 7, 9], W: int = 11
) -> np.ndarray:
    if W % 2 != 1:
        raise Exception("window W must be an odd number.")
    if max(K) > 20 or min(K) < 2:
        raise Exception("K must be within [2,20].")

    N = len(seq)
    results = np.zeros(N, dtype=np.float16)
    MAX_UNIQUE_HASHES = np.array([min(W - k + 1, 4**k) for k in K], dtype=int)
    # init
    kmer_hashes = [sequence_to_kmer_hashes(seq=seq[:W], k=k) for k in K]
    current_unique_n_hashes = np.array([len(set(l)) for l in kmer_hashes], dtype=int)
    results[int(W / 2)] = np.product(
        current_unique_n_hashes / MAX_UNIQUE_HASHES, dtype=np.float16
    )
    # loop
    for i in range(int(W / 2) + 1, N - int(W / 2)):
        # print(seq[i-int(W/2):i+int(W/2)+1])
        h = dna_to_base5(seq[i + int(W / 2)])
        # updates on np arrays by random access on cycling position -> no list construction
        for j, k in enumerate(K):
            position = int((i - int(W / 2) - 1) % (W - k + 1))
            last_position = int((i - int(W / 2) - 2) % (W - k + 1))
            old_hash = kmer_hashes[j][last_position]
            # update count of unique elements in k-mers (remove old one?)
            last_hash_unique = int(sum(kmer_hashes[j] == kmer_hashes[j][position])) == 1
            current_unique_n_hashes[j] -= int(last_hash_unique)
            new_hash = extend_hash(kmer_hash=old_hash, letter_dna5=h, k=k)
            kmer_hashes[j][position] = new_hash
            # update count of unique elements in k-mers (add new one?)
            new_hash_unique = int(sum(kmer_hashes[j] == new_hash)) == 1
            current_unique_n_hashes[j] += int(new_hash_unique)
        # calc new score, save
        results[i] = np.product(
            current_unique_n_hashes / MAX_UNIQUE_HASHES, dtype=np.float16
        )
    return results


# path reference
path_reference = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/hs37d5.fa"
)
path_regions = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/worst_offenders.bed"
)
path_bedgraph = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/worst_offenders.sequence_complexity.bedgraph"
)
# %%
regions = pd.read_csv(path_regions, sep="\t", header=None)
# for each region, extraxt the DNA sequence from the reference and calc sequence complexity

dict_sequences = {}
for row in tqdm(regions.iterrows()):
    chr, start, end, crID_mcrID = row[1]
    key = row[0]
    # extract DNA sequence
    cmd_extract_DNA = "samtools faidx {} {}:{}-{} | grep -v '>'".format(
        path_reference, chr, start, end
    )
    dna_sequence = (
        subprocess.check_output(cmd_extract_DNA, shell=True)
        .decode("utf-8")
        .replace("\n", "")
    )
    dict_sequences[key] = dna_sequence

# %%

# calc sequence complexity
dict_results = {}
for key, seq in tqdm(dict_sequences.items()):
    results = sequence_complexity(seq)
    dict_results[key] = results
# %%

# write the results to a bedgraph file.
# that is a simple text file with the following columns:
# chr, start, end, value
# value is the sequence complexity. start and end are generated
# by starting at j=start of region and j+=1 until all values are written.
# chr is the chromosome of the region.
with open(path_bedgraph, "w") as f:
    for row in tqdm(regions.iterrows()):
        key = row[0]
        chr, start, end, crID_mcrID = row[1]
        results = dict_results[key]
        for j, result in enumerate(results):
            f.write(f"{chr}\t{start+j}\t{start+j+1}\t{result}\n")

# %%
