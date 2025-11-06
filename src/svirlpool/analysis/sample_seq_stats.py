"""script to plot all sv sizes and letter compositions"""

import argparse
import gzip
import json
import logging as log
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pysam
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from intervaltree import Interval, IntervalTree
from tqdm import tqdm

from ..perform.seqomplexity import compute_complexity
from ..scripts import datatypes

log.basicConfig(level=log.INFO)


def heterogeneity_measure(S: str | Seq) -> float:
    n = len(S)  # Total number of letters
    if n == 0:
        return 0  # Edge case: empty string
    freq = Counter(S)  # Count occurrences of each letter
    m = len(freq)  # Number of unique letters
    contribution = sum(min(1, count / (n / m)) for count in freq.values())
    return contribution


def get_starts_ends(
    t: np.ndarray, x: np.ndarray, reference_start: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # compute a list of tuples of start and end positions for each cigar block
    read_start = x[0] if t[0] == 4 else 0
    x_read = x * ((t == 0) | (t == 1) | (t == 7) | (t == 8))
    x_read = np.cumsum(x_read) + read_start
    x_read_starts = np.array([read_start, *(x_read)[:-1]])
    x_read_ends = x_read
    #
    x_ref = x * ((t == 0) | (t == 2) | (t == 3) | (t == 7) | (t == 8))
    x_ref = np.cumsum(x_ref) + reference_start
    x_ref_starts = np.array([reference_start, *(x_ref)[:-1]])
    x_ref_ends = x_ref
    return x_read_starts, x_read_ends, x_ref_starts, x_ref_ends


def get_indels(
    alignment: pysam.AlignedSegment,
) -> tuple[
    list[tuple[int, int]],
    list[tuple[int, int]],
    list[tuple[int, int]],
    list[tuple[int, int]],
]:
    ref_start = alignment.reference_start
    t, x = zip(*alignment.cigartuples, strict=True)
    x = np.array(x)
    t = np.array(t)
    x_read_starts, x_read_ends, x_ref_starts, x_ref_ends = get_starts_ends(
        t=t, x=x, reference_start=ref_start
    )
    # get deletion ref start, ref end tuples
    deletions_ref = list(zip(x_ref_starts[(t == 2)], x_ref_ends[(t == 2)], strict=True))
    # get deletion read start, read end tuples
    deletions_read = list(
        zip(x_read_starts[(t == 2)], x_read_ends[(t == 2)], strict=True)
    )
    # get insertion ref start, ref end tuples
    insertions_ref = list(
        zip(x_ref_starts[(t == 1)], x_ref_ends[(t == 1)], strict=True)
    )
    # get insertion read start, read end tuples
    insertions_read = list(
        zip(x_read_starts[(t == 1)], x_read_ends[(t == 1)], strict=True)
    )
    return deletions_ref, deletions_read, insertions_ref, insertions_read


def parse_indels(
    aln: pysam.AlignedSegment, min_signal_size: int
) -> tuple[list[datatypes.ExtendedSVsignal], list[datatypes.ExtendedSVsignal]]:
    deletions_ref, deletions_read, insertions_ref, insertions_read = get_indels(aln)
    deletions = [
        datatypes.ExtendedSVsignal(
            ref_start=int(dell),
            ref_end=int(delr),
            read_start=int(rdell),
            read_end=int(rdelr),
            size=int(abs(delr - dell)),
            sv_type=int(1),
            chr=aln.reference_name,
            chrID=aln.reference_id,
            coverage=0,
            readname=aln.query_name,
            samplename="",
            forward=aln.is_forward,
        )
        for ((dell, delr), (rdell, rdelr)) in zip(
            deletions_ref, deletions_read, strict=True
        )
        if abs(delr - dell) >= min_signal_size
    ]
    insertions = [
        datatypes.ExtendedSVsignal(
            ref_start=int(insl),
            ref_end=int(insr),
            read_start=int(rinsl),
            read_end=int(rinsr),
            size=int(abs(rinsr - rinsl)),
            sv_type=int(0),
            chr=aln.reference_name,
            chrID=aln.reference_id,
            coverage=0,
            readname=aln.query_name,
            samplename="",
            forward=aln.is_forward,
        )
        for ((insl, insr), (rinsl, rinsr)) in zip(
            insertions_ref, insertions_read, strict=True
        )
        if abs(rinsr - rinsl) >= min_signal_size
    ]
    return insertions, deletions


def parse_bnds(
    aln: pysam.AlignedSegment, min_signal_size: int
) -> list[datatypes.ExtendedSVsignal]:
    # bnd left
    results: list[datatypes.ExtendedSVsignal] = []
    if aln.cigartuples[0][0] in (4, 5) and aln.cigartuples[0][1] >= min_signal_size:
        results.append(
            datatypes.ExtendedSVsignal(
                ref_start=aln.reference_start,
                ref_end=aln.reference_start,
                read_start=0,
                read_end=aln.cigartuples[0][1],
                size=aln.cigartuples[0][1],
                sv_type=3,
                chr=aln.reference_name,
                chrID=aln.reference_id,
                coverage=0,
                readname=aln.query_name,
                samplename="",
                forward=aln.is_forward,
            )
        )
    if aln.cigartuples[-1][0] in (4, 5) and aln.cigartuples[-1][1] >= min_signal_size:
        results.append(
            datatypes.ExtendedSVsignal(
                ref_start=aln.reference_end,
                ref_end=aln.reference_end,
                read_start=aln.query_length
                - (aln.cigartuples[-1][1] if aln.cigartuples[-1][0] == 5 else 0),
                read_end=aln.query_length,
                size=aln.cigartuples[-1][1],
                sv_type=4,
                chr=aln.reference_name,
                chrID=aln.reference_id,
                coverage=0,
                readname=aln.query_name,
                samplename="",
                forward=aln.is_forward,
            )
        )
    return results


def parse_signals(
    aln: pysam.AlignedSegment, min_signal_size: int
) -> list[datatypes.ExtendedSVsignal]:
    insertions, deletions = parse_indels(aln, min_signal_size)
    breakends = parse_bnds(aln, min_signal_size)
    return [*insertions, *deletions, *breakends]


def get_alt_sequence_for_insertion(
    ins: datatypes.ExtendedSVsignal, aln: pysam.AlignedSegment
) -> str:
    if ins.sv_type != 0:
        raise (ValueError("Not an insertion signal"))
    return aln.query_sequence[ins.read_start : ins.read_end]


# from typing import Generator
# def yield_all_signals(path_alignments:Path, min_signal_size:int) -> Generator[tuple[datatypes.ExtendedSVsignal, str], None, None]:
#     with pysam.AlignmentFile(path_alignments, "rb") as f:
#         for aln in tqdm(f):
#             if not aln.is_unmapped:
#                 for signal in parse_signals(aln=aln, min_signal_size=min_signal_size):
#                     seq = get_alt_sequence_for_insertion(ins=signal, aln=aln) if signal.sv_type == 0 else ""
#                     yield (signal, seq)


def add_ref_seq_to_deletions(
    signals: list[tuple[datatypes.ExtendedSVsignal, str]], reference: Path
) -> list[tuple[datatypes.ExtendedSVsignal, str]]:
    results = []
    with tempfile.TemporaryDirectory() as tmpdir:
        path_bed = Path(tmpdir) / "dels.bed"
        # write all deletions to a bed file
        with open(path_bed, "w") as bed:
            for signal, seq in signals:
                if signal.sv_type == 1:
                    print(f"{signal.chr}:{signal.ref_start}-{signal.ref_end}", file=bed)
        # use samtools faidx to extract all sequences defined in 'bed' from the reference sequence
        path_fa = Path(tmpdir) / "dels.fa"
        with open(path_fa, "w") as fa:
            subprocess.run(["samtools", "faidx", reference, "-r", path_bed], stdout=fa)
        # now create two iterators. The first iterates over path_fa, the second iterates over signals
        # if signal is a deletion, get the next sequence from path_fa and add it to the signal
        # is signal is not a deletion, continue
        it_signals = iter(signals)
        it_fa = SeqIO.parse(path_fa, "fasta")
        with tqdm(total=len(signals)) as pbar:
            for signal, seq in it_signals:
                if signal.sv_type == 1:
                    results.append((signal, next(it_fa).seq))
                else:
                    results.append((signal, seq))
                pbar.update(1)
    return results


def create_intervaltrees(path_repeats: Path) -> dict[str, IntervalTree]:
    # construct dict of interval tree from repeats
    it_repeats: dict[str, IntervalTree] = {}
    with open(path_repeats, "r") as f:
        for i, line in enumerate(f):
            chr, start, end = line.strip().split("\t")
            if chr not in it_repeats:
                it_repeats[chr] = IntervalTree()
            it_repeats[chr].add(Interval(int(start), int(end), i))
    return it_repeats
    # # iterate over signals, check if signal is in interval tree, if yes, add repeat information
    # for signal,dna in tqdm(signals):
    #     interval_signal = Interval(signal.ref_start, signal.ref_end) if signal.sv_type == 1 else Interval(signal.ref_start-1, signal.ref_start+1)
    #     if not signal.chr in it_repeats:
    #         raise ValueError(f"Chromosome {signal.chr} not found in repeat intervals")
    #     reps = it_repeats[signal.chr][interval_signal]
    #     if len(reps) > 0:
    #         signal.repeatID = reps[0].data


KVARIANTS = list(range(13, 2, -2))


def sequence_complexity(seq: str | Seq) -> float:
    n = len(seq)
    if n < 3:
        return 1.0
    m = min(n, 15)
    m = m if m % 2 == 1 else m - 1
    kvars = KVARIANTS
    if m < 15:
        # limit kvars so that kvars < m
        kvars = list(range(m, 2, -2))
    result = float(np.mean(compute_complexity(m, kvars, str(seq))))
    return result


def test_sequence_complexity() -> None:
    # generate a set of strings, empty, short, long, with different complexities
    import random
    import string

    test_strings = [
        "",
        "A",
        "AC",
        "ACG",
        "ACGT",
        "ACGT" * 10,
        "".join(random.choices(string.ascii_uppercase + string.digits, k=10**3)),
        "".join(random.choices(string.ascii_uppercase + string.digits, k=10**4)),
        "".join(random.choices(string.ascii_uppercase + string.digits, k=10**5)),
        "".join(random.choices(string.ascii_uppercase + string.digits, k=10**6)),
        "".join(random.choices(string.ascii_uppercase + string.digits, k=10**7)),
    ]
    for test_string in test_strings:
        print(f"Testing string with length {str(len(test_string))}")
        print(f"Complexity: {sequence_complexity(test_string)}")


# test_sequence_complexity()

# # %%

# DATA_ALIGNMENTS = {
#     "dummy":Path("/home/vinzenz/development/LRSV-detection/development/HG/trio/HG002/consensus.bam"),
#     "HG002":Path("/data/hdd/HG002/ont-r10.32x/minimap2.ont-r10.32x.bam"),
#     "HG003":Path("/data/hdd/HG003/ont-r10.32x/minimap2.ont-r10.32x.bam"),
#     "HG004":Path("/data/hdd/HG004/ont-r10.32x/minimap2.ont-r10.32x.bam")
# }
# REFERENCE = Path("/home/vinzenz/development/LRSV-detection/development/reference/hs37d5/hs37d5.fa")
# REPEATS = Path("/home/vinzenz/development/LRSV-detection/development/reference/hs37d5/human_hs37d5.trf.bed")
# MIN_SV_SIZE = 6
# SAMPLE = "HG002"

# ALL_SIGNALS_PATH = Path(f"/home/vinzenz/development/LRSV-detection/development/HG/trio/{SAMPLE}.all_signals.json.gz")


# #%%


# # signal.strength is the complexity of the sequence
# # signal.coverage is the GC content of the sequence
# # signal.repeatID is the repeat ID of the sequence
# # signal.sv_type is the type of the signal (0: insertion, 1: deletion, 3: breakend left, 4: breakend right)
# # signal.chrID can be used to store the minimum distance to the next SV signal

# log.info("parsing repeats..")
# intervaltrees = create_intervaltrees(REPEATS)
# all_signals : list[datatypes.ExtendedSVsignal] = []
# dict_region_del_idx:dict[str:int] = dict() # this dict links a region string (chr:start-end) to the index of the deletion in all_signals
# log.info("parsing SV signals..")
# with tempfile.TemporaryDirectory() as tmpdir:
#     path_regions = Path(tmpdir) / "dels.regions"
#     with open(path_regions, "w") as regions:
#         idx_signal = 0
#         with pysam.AlignmentFile(DATA_ALIGNMENTS["HG002"], "rb") as f:
#             for aln in tqdm(f):
#                 if not aln.is_unmapped:
#                     for signal in parse_signals(aln=aln, min_signal_size=MIN_SV_SIZE):
#                         if signal.sv_type == 0: # insertion
#                             seq = get_alt_sequence_for_insertion(ins=signal, aln=aln)
#                             signal.coverage = int(round(SeqUtils.gc_fraction(seq) * 100))
#                             signal.strength = sequence_complexity(seq)

#                         elif signal.sv_type == 1: # deletion
#                             region_string = f"{signal.chr}:{signal.ref_start}-{signal.ref_end}"
#                             print(region_string, file=regions)
#                             dict_region_del_idx[region_string] = idx_signal
#                         else: # break end
#                             pass
#                         signal_interval = Interval(signal.ref_start, signal.ref_start+1) if signal.sv_type != 1 else Interval(signal.ref_start, signal.ref_end)
#                         if signal.chr in intervaltrees:
#                             reps = intervaltrees[signal.chr][signal_interval.begin:signal_interval.end]
#                             if len(reps) > 0:
#                                 signal.repeatID = list(reps)[0].data if len(reps) > 0 else -1
#                         all_signals.append(signal)
#                         idx_signal += 1
#     # extract sequences from reference
#     log.info(f"extracting sequences from reference..")
#     path_fa = Path(tmpdir) / "dels.fa"
#     with open(path_fa, "w") as fa:
#         subprocess.run(["samtools", "faidx", REFERENCE, "-r", path_regions], stdout=fa)
#     log.info(f"annotating GC content and heterogeneity from extracted sequences..")
#     for seqrec in tqdm(SeqIO.parse(path_fa, "fasta")):
#         signal = all_signals[dict_region_del_idx[seqrec.id]]
#         signal.coverage = int(round(SeqUtils.gc_fraction(seqrec.seq) * 100))
#         signal.strength = sequence_complexity(seqrec.seq)

# log.info(f"sorting signals..")
# # iterate sorted signals (sorted by chrID, start, end)
# all_signals = sorted(all_signals, key=lambda x: (x.chrID, x.ref_start, x.ref_end))

# # re-use chrID to store the minimum distance to the next SV signal or the signal before
# # but only if their chr is the same, otherwise its 0
# log.info(f"annotating distances to next SV signal..")
# for i in tqdm(range(1,len(all_signals)-1)):
#     dist_to_prior = all_signals[i].ref_start - all_signals[i-1].ref_end if all_signals[i].chrID == all_signals[i-1].chrID else -1
#     dist_to_next = all_signals[i+1].ref_start - all_signals[i].ref_end if all_signals[i].chrID == all_signals[i+1].chrID else -1
#     if dist_to_next >= 0 and dist_to_prior >= 0:
#         all_signals[i].chrID = min(dist_to_next, dist_to_prior)
#     elif dist_to_next >= 0:
#         all_signals[i].chrID = dist_to_next
#     elif dist_to_prior >= 0:
#         all_signals[i].chrID = dist_to_prior
#     else:
#         all_signals[i].chrID = -1

# #%% save all signals to a compressed file using json and gzip


# with gzip.open(ALL_SIGNALS_PATH, "wt") as f:
#     json.dump([signal.unstructure() for signal in all_signals], f)

# # %% create a dataframe with columns: sy_type, size, GC content, complexity, repeatID
# # goes to other script
# SVTYPES = {0:"INS", 1:"DEL", 3:"BNDL", 4:"BNDR"}

# import pandas as pd
# df = pd.DataFrame(
#     data=[
#         (SVTYPES[signal.sv_type], signal.size, signal.coverage, signal.strength, signal.repeatID > -1, signal.chrID)
#         for signal in all_signals
#     ],
#     columns=["sv_type", "size", "GC content", "complexity", "repeat", "distance"]
# )

# DATAFRAME_PATH = Path(f"/home/vinzenz/development/LRSV-detection/development/HG/trio/{SAMPLE}.dataframe.tsv.gz")
# with gzip.open(DATAFRAME_PATH, "wt") as f:
#     df.to_csv(f, sep="\t", index=False)
# %%


def extract_all_signals(
    path_alignments: Path,
    path_repeats: Path,
    path_reference: Path,
    min_signal_size: int,
    path_output: Path,
) -> None:
    log.info("parsing repeats..")
    intervaltrees = create_intervaltrees(path_repeats)
    all_signals: list[datatypes.ExtendedSVsignal] = []
    dict_region_del_idx: dict[
        str, int
    ] = {}  # this dict links a region string (chr:start-end) to the index of the deletion in all_signals
    log.info("parsing SV signals..")
    with tempfile.TemporaryDirectory() as tmpdir:
        path_regions = Path(tmpdir) / "dels.regions"
        with open(path_regions, "w") as regions:
            idx_signal = 0
            with pysam.AlignmentFile(path_alignments, "rb") as f:
                for aln in tqdm(f):
                    if aln.is_unmapped or aln.is_secondary:
                        continue
                    else:
                        for signal in parse_signals(
                            aln=aln, min_signal_size=min_signal_size
                        ):
                            if signal.sv_type == 0:  # insertion
                                seq = get_alt_sequence_for_insertion(
                                    ins=signal, aln=aln
                                )
                                signal.coverage = int(
                                    round(SeqUtils.gc_fraction(seq) * 100)
                                )
                                signal.strength = sequence_complexity(seq)

                            elif signal.sv_type == 1:  # deletion
                                region_string = (
                                    f"{signal.chr}:{signal.ref_start}-{signal.ref_end}"
                                )
                                print(region_string, file=regions)
                                dict_region_del_idx[region_string] = idx_signal
                            else:  # break end
                                pass
                            signal_interval = (
                                Interval(signal.ref_start, signal.ref_start + 1)
                                if signal.sv_type != 1
                                else Interval(signal.ref_start, signal.ref_end)
                            )
                            if signal.chr in intervaltrees:
                                reps = intervaltrees[signal.chr][
                                    signal_interval.begin : signal_interval.end
                                ]
                                if len(reps) > 0:
                                    signal.repeatID = (
                                        list(reps)[0].data if len(reps) > 0 else -1
                                    )
                            all_signals.append(signal)
                            idx_signal += 1
        # extract sequences from reference
        log.info("extracting sequences from reference..")
        path_fa = Path(tmpdir) / "dels.fa"
        with open(path_fa, "w") as fa:
            subprocess.run(
                ["samtools", "faidx", str(path_reference), "-r", str(path_regions)],
                stdout=fa,
            )
        log.info("annotating GC content and heterogeneity from extracted sequences..")
        for seqrec in tqdm(SeqIO.parse(path_fa, "fasta")):
            signal = all_signals[dict_region_del_idx[seqrec.id]]
            signal.coverage = int(round(SeqUtils.gc_fraction(seqrec.seq) * 100))
            signal.strength = sequence_complexity(seqrec.seq)

    log.info("sorting signals..")
    # iterate sorted signals (sorted by chrID, start, end)
    all_signals = sorted(all_signals, key=lambda x: (x.chrID, x.ref_start, x.ref_end))

    # re-use chrID to store the minimum distance to the next SV signal or the signal before
    # but only if their chr is the same, otherwise its 0
    log.info("annotating distances to next SV signal..")
    for i in tqdm(range(1, len(all_signals) - 1)):
        dist_to_prior = (
            all_signals[i].ref_start - all_signals[i - 1].ref_end
            if all_signals[i].chrID == all_signals[i - 1].chrID
            else -1
        )
        dist_to_next = (
            all_signals[i + 1].ref_start - all_signals[i].ref_end
            if all_signals[i].chrID == all_signals[i + 1].chrID
            else -1
        )
        if dist_to_next >= 0 and dist_to_prior >= 0:
            all_signals[i].chrID = min(dist_to_next, dist_to_prior)
        elif dist_to_next >= 0:
            all_signals[i].chrID = dist_to_next
        elif dist_to_prior >= 0:
            all_signals[i].chrID = dist_to_prior
        else:
            all_signals[i].chrID = -1
    log.info(f"writing signals to {path_output}")
    with gzip.open(path_output, "wt") as f:
        for signal in all_signals:
            print(json.dumps(signal.unstructure()), file=f)


def run(args, **kwargs) -> None:
    extract_all_signals(
        path_alignments=args.alignments,
        path_repeats=args.repeats,
        path_reference=args.reference,
        min_signal_size=args.min_signal_size,
        path_output=args.output,
    )
    return


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Reads crs container objects and create Consensus objects that are written to the output database."
    )
    parser.add_argument(
        "-a",
        "--alignments",
        type=Path,
        required=True,
        help="Path to the alignments file.",
    )
    parser.add_argument(
        "-r", "--repeats", type=Path, required=True, help="Path to the repeats file."
    )
    parser.add_argument(
        "-f",
        "--reference",
        type=Path,
        required=True,
        help="Path to the reference file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output gzipped json dumped file containing all signals.",
    )
    parser.add_argument(
        "-s",
        "--min_signal_size",
        type=int,
        default=6,
        required=False,
        help="Minimum signal size.",
    )
    return parser


def main() -> None:
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
