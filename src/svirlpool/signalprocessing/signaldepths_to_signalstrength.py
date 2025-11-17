# %%

import argparse
import csv
import gzip
import json
import logging
import multiprocessing
import shlex
import subprocess
import sys
import tempfile
from bisect import insort_left
from pathlib import Path

import numpy as np
import numpy.typing as npt

from ..signalprocessing import alignments_to_rafs
from ..util import datatypes, util

# logger
log = logging.getLogger(__name__)


# creates a score between 0 and 1 for the similarity of two signal sizes.
def size_score(difference: int, min_sv_length: int, max_sv_length) -> float:
    if difference > max_sv_length:
        return 0.0
    if difference < min_sv_length:
        return 1.0
    # return float(1/(difference*(100/min_sv_length)))
    return float(min_sv_length / (difference * 100))


def size_similarity_score(a: npt.NDArray[np.int32], b: npt.NDArray[np.int32]) -> float:
    """computes a score between 0 and 1 for the similarity of two sizes. The function is symmetric and convex."""
    min_sv_size = 0.01 * min(a[4], b[4])
    max_sv_size = 0.10 * min(a[4], b[4])
    score = size_score(abs(a[4] - b[4]), min_sv_size, max_sv_size)
    return score


def distance_collapsed_sv_signals(
    a: npt.NDArray[np.int32], b: npt.NDArray[np.int32], repeat_factor: float = 10.0
) -> np.int32:
    # if a and b are on the same chromosome (item 0)
    # then return the distance between the two signals
    if a[0] == b[0]:
        a, b = (b, a) if b[1] < a[1] else (a, b)
        # if a and b have an overlap in the interval of their last two elements, then return 0
        try:
            if b[-2] == a[-2] and a[-2] >= 0:  # both share same repeat
                # return np.int32(0)
                return np.int32(max(0, b[1] - a[2]) / repeat_factor)
        except TypeError:
            log.warning(f"TypeError: {str(a)} {str(b)}")
        return np.int32(max(0, b[1] - a[2]))
    else:
        return np.int32(2**31 - 1)


def genomic_distance_score(d: int, flatness: float = 15.0) -> float:
    """Flatness factor: Higher flatness means that the distance score is more tolerant to distance."""
    if d == 0:
        return 1.0
    return 1 / ((d + flatness) / flatness)


def delta(
    subj: npt.NDArray[np.int32], obj: npt.NDArray[np.int32], flatness: float
) -> float:
    """returns 1 if no distance between two loci, else x >= 0"""
    distance = np.uint32(2**31 - 1)
    size_score = 0.0
    if subj[0] != obj[0] or subj[7] == obj[7]:
        return 0  # closeness is 0 if they are on different chromosomes or are the same read
    # rewrite this to use DELL and DELR instead of only DEL
    # INS is now code 0, DELL is code 1, DELR is code 2, BNDL is code 3, BNDR is code 4
    match (subj[3], obj[3]):
        case (0, 0):  # INS vs INS
            distance = distance_collapsed_sv_signals(subj, obj)
            size_score = size_similarity_score(subj, obj)
        case (0, 1):  # INS vs DELL
            pass
        case (0, 2):  # INS vs DELR
            pass
        case (0, 3):  # INS vs BNDL
            distance = distance_collapsed_sv_signals(subj, obj)
        case (0, 4):  # INS vs BNDR
            distance = distance_collapsed_sv_signals(subj, obj)
        case (1, 0):  # DELL vs INS
            pass
        case (1, 1):  # DELL vs DELL
            distance = distance_collapsed_sv_signals(subj, obj)
            size_score = size_similarity_score(subj, obj)
        case (1, 2):  # DELL vs DELR
            pass
        case (1, 3):  # DELL vs BNDL
            pass
        case (1, 4):  # DELL vs BNDR
            distance = distance_collapsed_sv_signals(subj, obj)
        case (2, 0):  # DELR vs INS
            pass
        case (2, 1):  # DELR vs DELL
            pass
        case (2, 2):  # DELR vs DELR
            distance = distance_collapsed_sv_signals(subj, obj)
            size_score = size_similarity_score(subj, obj)
        case (2, 3):  # DELR vs BNDL
            distance = distance_collapsed_sv_signals(subj, obj)
        case (2, 4):  # DELR vs BNDR
            pass
        case (3, 0):  # BNDL vs INS
            pass
        case (3, 1):  # BNDL vs DELL
            pass
        case (3, 2):  # BNDL vs DELR
            pass
        case (3, 3):  # BNDL vs BNDL
            distance = distance_collapsed_sv_signals(subj, obj)
        case (3, 4):  # BNDL vs BNDR
            pass
        case (4, 0):  # BNDR vs INS
            pass
        case (4, 1):  # BNDR vs DELL
            pass
        case (4, 2):  # BNDR vs DELR
            pass
        case (4, 3):  # BNDR vs BNDL
            pass
        case (4, 4):  # BNDR vs BNDR
            distance = distance_collapsed_sv_signals(subj, obj)
    if distance < 0:
        return 1.0
    elif distance == 0:
        return 1.0
    else:
        distance_score = genomic_distance_score(distance, flatness)
        return max(distance_score, size_score)


def pseudo_convolve(
    np_signals: npt.NDArray[np.int32],
    kernelsize: int,
    flatness: float,
    warning_collapsed_region_size: int = 2_000,
    high_density_threshold: int = 200,
    high_density_score: float = 10.0,
) -> np.ndarray:
    """convolves over a numpy array of signals and returns the sum of weighted distances
    by a kernel function (radius) which is dynamically expanded if it overlaps any repeats.
    the flatness factor describes how tolerant the distance is measured, e.g. 4.0 means that
    two SV signals with 4 bp distance still have 50% spatial distance score.

        Args:
            np_signals: An array with each row of: refID, start, end, svsize, svID, repeatIDstart, repeatIDend
            kernelsize: int. The radius of the kernel function
            flatness: float. The flatness factor
            high_density_threshold: int. If region contains more signals than this, assign high_density_score
            high_density_score: float. Score to assign to high-density regions
    """
    values = np.zeros(np_signals.shape[0], dtype=float)
    j, k = 0, 0
    last_warned_index = 0
    i = 0

    while i < np_signals.shape[0]:
        # extend k to the right
        while (
            k < np_signals.shape[0]
            and distance_collapsed_sv_signals(np_signals[i, :], np_signals[k, :])
            <= kernelsize
        ):
            k += 1
        # extend j to the left
        while (
            j < i
            and distance_collapsed_sv_signals(np_signals[j, :], np_signals[i, :])
            > kernelsize
        ):
            j += 1

        region_size = k - j

        # Check for high-density region first
        if region_size > high_density_threshold:
            # Assign high density score to entire region at once
            values[j:k] = high_density_score
            if k > last_warned_index:
                last_warned_index = k
                log.info(
                    f"High-density region {str(j)} to {str(k)} ({str(np_signals[j, :4])} to {str(np_signals[k - 1, :4])}) with {str(region_size)} signals. Assigned score {high_density_score}."
                )
            # Skip to end of this region
            i = k
            j = k
            continue

        # Normal processing for lower-density regions
        if region_size > warning_collapsed_region_size:
            if k > last_warned_index:
                last_warned_index = k
                log.warning(
                    f"region {str(j)} to {str(k)} ({str(np_signals[j, :4])} to {str(np_signals[k - 1, :4])}) is very large ({str(region_size)}). (>{str(warning_collapsed_region_size)})."
                )

        for l in range(j, k):  # noqa: E741
            if l == i:
                continue
            next_val = delta(subj=np_signals[i], obj=np_signals[l], flatness=flatness)
            values[i] += next_val

        i += 1

    return values


def process_func(kwargs) -> str:
    return pseudo_convolve(**kwargs)


def slop_bed(input: Path, reference: Path, output: Path, bp: int) -> None:
    # create a .genome file (two columns: chr, length) from reference.fa.fai
    # first check if .fai is present
    path_fai = Path(str(reference) + ".fai")
    if not path_fai.exists():
        raise FileNotFoundError(
            f"Could not find {str(path_fai)}. Make sure to generate an index for your reference with 'samtools faidx'"
        )
    # create a .genome file
    tmp_genome = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".genome")
    with open(path_fai, "r") as f:
        with open(tmp_genome.name, "w") as g:
            for line in f:
                line_split = line.strip().split("\t")
                print(line_split[0], line_split[1], sep="\t", file=g)
    cmd_bedtools_slop = (
        f"bedtools slop -i {str(input)} -g {str(tmp_genome.name)} -b {bp}"
    )
    with open(output, "w") as f:
        subprocess.check_call(shlex.split(cmd_bedtools_slop), stdout=f)
    tmp_genome.close()


def parse_signaldepths_to_npsignals(input: Path) -> npt.NDArray[np.int32]:
    """Reduces signaldepths to numpy array. Readnames are relplaced with hashes and the following columns are dropped: read_start, read_end, sampleID, forward, chr. \
    The resulting numpy array has the following columns: chrID, start, end, svType, SVsize, depth, repeatID, read_hash."""
    np_signals = np.loadtxt(
        input, usecols=(0, 1, 2, 3, 6, 11, 12), delimiter="\t", dtype=int
    )
    # build dictionary of readname to hash
    readnames = np.loadtxt(input, usecols=(7), delimiter="\t", dtype=str)
    # create a dictionary of readnames to hash, where hash is an incremental index with observation
    readnames_dict = {readname: k for k, readname in enumerate(np.unique(readnames))}
    # np_strdata = pd.read_csv(tmp_signals_indexed,sep='\t',usecols=(7,8,9),header=None,dtype={7:str,8:int,9:str}).to_numpy()
    # compute hash value of each sampleID
    # readHashes = pd.read_csv(tmp_signals_indexed,sep='\t',usecols=[7],header=None,dtype=str).iloc[:,0].apply(lambda x: hash(x))
    # add readHashes to np_signals
    readHashes = np.array([readnames_dict[x] for x in readnames])
    return np.hstack((np_signals, np.array(readHashes).reshape(-1, 1)))


def annotate_signals_with_repeats(
    repeats: Path, signals: Path, output: Path, tmp_dir_path: Path | None = None
) -> None:
    # assume reapeats has chrIDs in first column
    # assume repeats is non-overlapping
    # assumes that repeats have an ID in the last (4th) column
    tmp_intersected = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path,
        delete=False if tmp_dir_path else True,
        suffix=".tmp_intersected.signals.bed",
    )
    cmd_intersection = f"bedtools intersect -wao -a {str(signals)} -b {str(repeats)}"
    # intersection creates 5 new columns: chr_repeat, start_repeat, end_repeat, repeatID, overlap
    with open(tmp_intersected.name, "w") as f:
        subprocess.check_call(shlex.split(cmd_intersection), stdout=f)
    # now drop all new columns except the repeatID
    with open(tmp_intersected.name, "r") as f:
        with open(output, "w") as g:
            for line in f:
                line_split = line.strip().split("\t")
                repeatID = line_split[-2]
                old_line = line_split[:-5]
                if repeatID == ".":
                    repeatID = "-1"
                print(*old_line, repeatID, sep="\t", file=g)


# def replace_sv_events_with_codes(input:Path,output:Path) -> None:
#     cmd_sed = f"sed 's/INS/0/g; s/DELL/1/g; s/DELR/2/g; s/BNDL/3/g; s/BNDR/4/g' {str(input)}"
#     with open(output,'w') as f:
#         subprocess.check_call(shlex.split(cmd_sed),stdout=f)


def split_DEL_to_DELL_and_DELR(
    path_signaldepths: Path, path_output: Path, threads: int
) -> None:
    cache = []
    last_chrID = None
    tmp_output = tempfile.NamedTemporaryFile(suffix=".tmp.split_dels.bed", delete=True)
    with open(tmp_output.name, "w") as out:
        for line in gzip.open(path_signaldepths, "rt"):
            current_signal = line.strip().split("\t")
            if last_chrID and last_chrID != current_signal[0]:
                # save all items in cache to output
                for cached_signal in cache:
                    print(*cached_signal, sep="\t", file=out)
                cache = []
            while len(cache) > 0 and int(current_signal[1]) > int(cache[0][1]):
                cached_signal = cache.pop(0)
                print(*cached_signal, sep="\t", file=out)
            if current_signal[3] == 1:
                # create DELL
                dell = current_signal.copy()
                dell[2] = str(int(dell[1]) + 1)
                dell[3] = str(1)
                print(*dell, sep="\t", file=out)
                delr = current_signal.copy()
                delr[1] = str(int(delr[2]) - 1)
                delr[3] = str(2)
                insort_left(a=cache, x=delr, key=lambda x: (int(x[0]), int(x[1])))
                # cache.append(delr)
                # cache = sorted(cache, key=lambda x: int(x[1]))
            else:
                print(*current_signal, sep="\t", file=out)
            last_chrID = current_signal[0]
        if len(cache) > 0:
            for cached_signal in cache:
                print(*cached_signal, sep="\t", file=out)
    alignments_to_rafs.compress_and_index_bedlike(
        sort_numerically=True,
        input=tmp_output.name,
        output=path_output,
        threads=threads,
    )


def signaldepths_to_signalstrength(
    reference: Path,
    signaldepths: Path,
    repeats: Path,
    output: Path,
    flatness_distance: float,
    threads: int,
    kernel_signal_radius: int,
    slop_repeats_radius: int,
    tmp_dir_path: Path | None = None,
) -> None:
    csv.field_size_limit(sys.maxsize)
    # --- load and prepare data --- #

    # 1. transform repeats to contain chrID instead of chr and add a column with the index
    # save to temporary file

    # load depth track
    # load last column to numpy array
    # flatness_noise = flatness_signal * (kernel_noise_radius / kernel_signal_radius)
    # skips 5,6,7 = readID,sampleID,chr
    # np_signals = pd.read_csv(tmp_signals_indexed,sep='\t',usecols=(0,1,2,3,6,7,11,12),header=None,low_memory=False,dtype=object).to_numpy()

    # repeats might contain more than 3 columns, however, there can only be 3 columns.
    # copy only the first three columns to the tmp_repeats file
    tmp_repeats = tempfile.NamedTemporaryFile(
        suffix=".repeats.bed", dir=tmp_dir_path, delete=False if tmp_dir_path else True
    )
    with open(repeats, "r") as f:
        with open(tmp_repeats.name, "w") as g:
            for line in f:
                line_split = line.strip().split("\t")
                print(*line_split[:3], sep="\t", file=g)

    log.info("slop repeats and replace chr in first column with chrID..")

    tmp_repeats_slopped = tempfile.NamedTemporaryFile(
        suffix=".repeats.slop.bed",
        dir=tmp_dir_path,
        delete=False if tmp_dir_path else True,
    )
    slop_bed(
        input=tmp_repeats.name,
        reference=reference,
        output=tmp_repeats_slopped.name,
        bp=slop_repeats_radius,
    )

    tmp_repeats_chrID = tempfile.NamedTemporaryFile(
        suffix=".chrID.bed", dir=tmp_dir_path, delete=False if tmp_dir_path else True
    )
    # repeats need to have chrID instead of chr in chr column to fit with signals, adds an ID as the last column
    util.bed_chr_to_chrID(
        input=tmp_repeats_slopped.name,
        reference=reference,
        output=tmp_repeats_chrID.name,
        add_id=True,
    )

    # split DELs to DELLs and DELRs
    tmp_signaldepths_split_dels = tempfile.NamedTemporaryFile(
        suffix=".split_dels.bed",
        dir=tmp_dir_path,
        delete=False if tmp_dir_path else True,
    )
    split_DEL_to_DELL_and_DELR(
        path_output=tmp_signaldepths_split_dels.name,
        path_signaldepths=signaldepths,
        threads=threads,
    )

    log.info("annotate signals with repeats..")
    tmp_signals_with_repeats = tempfile.NamedTemporaryFile(
        suffix=".signals_with_repeats.bed",
        dir=tmp_dir_path,
        delete=False if tmp_dir_path else True,
    )
    annotate_signals_with_repeats(
        tmp_dir_path=tmp_dir_path,
        output=tmp_signals_with_repeats.name,
        signals=signaldepths,
        repeats=tmp_repeats_chrID.name,
    )

    log.info("parse signaldepths to numpy array..")
    np_signals = parse_signaldepths_to_npsignals(tmp_signals_with_repeats.name)

    log.info("convolve signal..")
    # parallelize this with multiprocessing. split np_signals into chunks (each chunk is one chromosome, row 0 is chrID)
    # and run pseudo_convolve on each chunk. then concatenate the results.
    # create a set or job arguments for pseudo_convolve
    jobs_args = [
        {
            "np_signals": np_signals[np_signals[:, 0] == refID],
            "kernelsize": kernel_signal_radius,
            "flatness": flatness_distance,
        }
        for i, refID in enumerate(set(np_signals[:, 0]))
    ]
    if threads > 1:
        with multiprocessing.Pool(threads) as pool:
            values_signal = np.concatenate(
                pool.map(process_func, jobs_args, chunksize=1)
            )
    else:
        values_signal = np.concatenate([
            pseudo_convolve(**kwargs) for kwargs in jobs_args
        ])
    # add values_signal to signaldepths. Open signals and to each row add the value_signal and write to output
    log.info("write output..")
    tmp_output = tempfile.NamedTemporaryFile(delete=True, suffix=".tsv")
    with open(tmp_signals_with_repeats.name, "r") as g:
        with open(tmp_output.name, "w") as f:
            writer = csv.writer(f, delimiter="\t")
            for i, line in enumerate(g):
                line_split = line.strip().split("\t")
                eSVs: datatypes.ExtendedSVsignal = (
                    parse_signalstrengths_to_extendedSVsignal(line_split)
                )
                eSVs.strength = float(values_signal[i])
                chrID, start, end = eSVs.chrID, eSVs.ref_start, eSVs.ref_end
                writer.writerow([chrID, start, end, json.dumps(eSVs.unstructure())])
                # print(*line_split[:3],json.dumps(eSVs.unstructure()),sep='\t',file=f)
    # compress output
    alignments_to_rafs.compress_and_index_bedlike(
        input=tmp_output.name, output=output, threads=threads, sort_numerically=True
    )


def parse_signalstrengths_to_extendedSVsignal(
    _l: list[object],
) -> datatypes.ExtendedSVsignal:
    # 0      1          2        3        4           5         6     7         8           9        10   11        12        13
    # chrID, ref_start, ref_end, sv_type, read_start, read_end, size, readname, samplename, forward, chr, coverage, repeatID, strength
    return datatypes.ExtendedSVsignal(
        chr=str(_l[10]),
        chrID=int(_l[0]),
        ref_start=int(_l[1]),
        ref_end=int(_l[2]),
        sv_type=int(_l[3]),
        read_start=int(_l[4]),
        read_end=int(_l[5]),
        size=int(_l[6]),
        readname=str(_l[7]),
        samplename=str(_l[8]),
        forward=int(_l[9]),
        coverage=int(_l[11]) if len(_l) > 11 else 0,
        repeatID=int(_l[12]) if len(_l) > 12 else -1,
        strength=float(_l[13]) if len(_l) > 13 else 0.0,
    )


def run(args, **kwargs):
    signaldepths_to_signalstrength(
        signaldepths=args.input,
        repeats=args.repeats,
        reference=args.reference,
        output=args.output,
        flatness_distance=args.flatness_distance,
        kernel_signal_radius=args.kernel_signal_radius,
        threads=args.threads,
        slop_repeats_radius=args.slop_repeats_radius,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Convolves over the singal to get the localized signalstrength to be used to create candidate regions."
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="signaldepth file."
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to reference in uncompressed .fa(sta) format.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to bed-like file that augments signaldepths by signalstrength information.",
    )
    parser.add_argument(
        "-R",
        "--repeats",
        type=Path,
        required=True,
        help="Path to bed file with repeats.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=8,
        help="Number of threads to use.",
    )
    parser.add_argument(
        "-s",
        "--flatness_distance",
        type=float,
        required=False,
        default=15.0,
        help="Factor that shapes the flatness of the distance function bewteen genomic loci.",
    )
    parser.add_argument(
        "-k",
        "--kernel_signal_radius",
        type=int,
        required=False,
        default=1000,
        help="Radius of kernel for signal.",
    )
    parser.add_argument(
        "--slop_repeats_radius",
        type=int,
        required=False,
        default=50,
        help="Radius of slop for repeats.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
