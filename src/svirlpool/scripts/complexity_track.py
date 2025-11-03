# for seq in fasta
# compute rolling complexity
#

# %%

import argparse
import locale
import typing
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from tqdm import tqdm

locale.setlocale(locale.LC_ALL, "")

# %%


def plot_seq_with_intensities(seq, scores: np.ndarray):
    scores /= scores.max()
    figure = plt.figure(figsize=(len(seq) / 5, 1))
    axes = figure.add_subplot(111)
    sns.heatmap(
        data=scores.reshape(1, -1),
        annot=np.array([*seq]).reshape(1, -1),
        fmt="",
        xticklabels=[],
        yticklabels=[],
        cmap="Reds",
        vmin=0,
        vmax=1,
        square=True,
        ax=axes,
        cbar=False,
    )
    plt.show()


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
    for i in tqdm(range(int(W / 2) + 1, N - int(W / 2))):
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


# %%
def calc_sequence_complexity_for_reference(
    input: Path,
    output: Path,
    output_index: Path = Path(""),
    K: typing.List[int] = list(range(2, 11, 1)),
    W: int = 13,
):
    """calculates the sequence complexity with a given window size W for a provided fasta file.

    Args:
        input (Path): Path to fasta file with references.
        output (Path): Path to numpy memmap file (e.g. hg19.complexity.npz).
        output_index (Path, optional): If provided, the index of the memmap will be stored here. Otherwise, it will be at [output].index.
        W (int, optional): Window size must be odd. Defaults to 13.

    Raises:
        Exception: If W is even.
    """
    if max(K) >= W:
        raise Exception(f"W must be greater than max K with W={W} and K={K}")
    if W % 2 != 1:
        raise Exception(f"W={W} must be an odd number.")
    input_fai_path = input.with_suffix(input.suffix + ".fai")
    if output_index == Path(""):
        output_index = output.with_suffix(output.suffix + ".index")

    seqlengths = (
        pd.read_csv(input_fai_path, sep="\t", header=None)
        .iloc[:, 1]
        .astype(int)
        .to_numpy()
    )
    seqlengths = np.concatenate([np.array([0]), seqlengths])
    seqlensums = seqlengths.cumsum()
    N = int(seqlengths.sum())
    findex = []
    M = np.memmap(output, dtype=np.float16, mode="w+", shape=N)
    for i, record in enumerate(SeqIO.parse(input, "fasta")):
        print(
            f"calc complexity on sequence {record.id} with length {len(record.seq):,}"
        )
        findex.append(seqlensums[i])
        M[seqlensums[i] : seqlensums[i + 1]] = sequence_complexity(
            seq=record.seq.upper(), W=W, K=K
        )
        M.flush()
    findex.append(N)
    print("writing index..")
    with open(output_index, "w") as f:
        for value in findex:
            print(value, file=f)


def run(args, **kwargs):
    calc_sequence_complexity_for_reference(
        input=args.input,
        output=args.output,
        output_index=args.output_index,
        W=args.W,
        K=args.K,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="calculates the sequence complexity with a given window size W for a provided fasta file."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to fasta file with references.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to numpy memmap file (e.g. hg19.complexity.npz).",
    )
    parser.add_argument(
        "-n",
        "--output_index",
        type=Path,
        required=False,
        default=Path(""),
        help="If provided, the index of the memmap will be stored here. Otherwise, it will be at [output].index.",
    )
    parser.add_argument(
        "-w",
        "--W",
        type=int,
        required=False,
        default=13,
        help="Window size must be odd. Defaults to 13.",
    )
    parser.add_argument(
        "-k",
        "--K",
        type=typing.List[int],
        required=False,
        default=[2, 3, 4, 5, 10],
        help="All k for the checked k-mers. Choose 1 < k < W+1 < 21",
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
# input=Path('/home/vinzenz/development/LRSV-detection/development/test/refs.fasta')
# output=Path('/home/vinzenz/development/LRSV-detection/development/test/refs.npy')
# calc_sequence_complexity_for_reference(input=input,output=output,K=[2,3,4,5,10])
