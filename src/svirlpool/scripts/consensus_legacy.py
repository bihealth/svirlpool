# %%

import argparse
import json
import shlex
import sqlite3
import subprocess
import tempfile
import time
import typing
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from random import sample as randomsample
from shlex import split

import cattrs
import logzero
import numpy as np
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from logzero import logger as log

from . import alignments_to_rafs, datatypes, util


def yield_crs_per_region(
    region: tuple[int, int, int], filename: Path
) -> typing.Generator[datatypes.CandidateRegion, None, None]:
    cmd_tabix = f"tabix -0 -f {str(filename)} {region[0]}:{region[1]}-{region[2]}"
    process = subprocess.Popen(shlex.split(cmd_tabix), stdout=subprocess.PIPE)
    for line in process.stdout:
        raw = (
            line.decode("utf-8")
            .strip()
            .split("\t")[-1]
            .strip()
            .replace('""', '"')[1:-1]
        )
        yield cattrs.structure(json.loads(raw), datatypes.CandidateRegion)


def regions_from_intervals(intervals: np.ndarray) -> list[tuple[int, int, int]]:
    return [(int(row[5]), int(row[6]), int(row[7])) for row in intervals]


def load_crs_per_regions(
    regions: list[tuple[int, int, int]], path_crs: Path
) -> list[datatypes.CandidateRegion]:
    # load all candidate regions that are in the intervals
    return list(
        set(
            [
                cr
                for region in regions
                for cr in yield_crs_per_region(region=region, filename=path_crs)
            ]
        )
    )


def extract_signals(rass: list[datatypes.ReadAlignmentSignals]) -> np.ndarray:
    """Create a numpy array with columns: ref_pos, ref_end,"""
    readnames = {raf.read_name: i for i, raf in enumerate(rass)}
    signals = np.array(
        [
            np.array(
                [
                    sv.ref_start,
                    sv.ref_end,
                    sv.size,
                    sv.sv_type,
                    k,
                    readnames[ras.read_name],
                ],
                dtype=np.int32,
            )
            for k, ras in enumerate(rass)
            for sv in ras.SV_signals
        ],
        dtype=np.int32,
    )
    if signals.size == 0:
        return np.array([])
    signals = signals[signals[:, 0].argsort()]
    return signals


def sv_neighbors(
    svSignal: datatypes.SVsignal, signals: np.ndarray, radius: int
) -> np.ndarray:
    # check if signals is empty. If so, return an empty array
    if signals.size == 0:
        return np.array([])
    if svSignal.sv_type == 0:  # insertion
        # check if there are signals of the same type in the vicinity
        start, end = svSignal.ref_start - radius, svSignal.ref_start + radius
        return signals[
            (signals[:, 3] == 0) & (signals[:, 0] >= start) & (signals[:, 0] <= end)
        ]
    if svSignal.sv_type == 1:  # deletion
        # check if there are signals of the same type in the vicinity
        start_left, end_left = svSignal.ref_start - radius, svSignal.ref_start + radius
        start_right, end_right = svSignal.ref_end - radius, svSignal.ref_end + radius
        mask_left = (
            (signals[:, 3] == 1)
            & (signals[:, 0] >= start_left)
            & (signals[:, 0] <= end_left)
        )
        mask_right = (
            (signals[:, 3] == 1)
            & (signals[:, 0] >= start_right)
            & (signals[:, 0] <= end_right)
        )
        return signals[mask_left | mask_right]
    if svSignal.sv_type == 3 or svSignal.sv_type == 4:  # break end
        # check if there are signals of the same type in the vicinity
        start, end = svSignal.ref_start - radius, svSignal.ref_start + radius
        return signals[
            ((signals[:, 3] == 3) | (signals[:, 3] == 4))
            & (signals[:, 0] >= start)
            & (signals[:, 0] <= end)
        ]
    raise ValueError("sv must be of type SVsignal and must have sv_type 0, 1, 3 or 4.")


# =============================================================================
#  clustering
# =============================================================================

import numpy as np
import pandas as pd
import plotly.express as px
from sklearn.cluster import AgglomerativeClustering


def row_to_annotated_readname(row: np.ndarray) -> str:
    return f"{row[8]}.{row[9]}"


def annotated_readname_to_readname_and_sampleID(
    annotated_readname: str,
) -> tuple[str, int]:
    x = annotated_readname.split(".")
    return x[0], int(x[1])


def sum_dels_per_read(
    cr: datatypes.CandidateRegion, readnames: list[str]
) -> tuple[list[str], dict[str, int]]:
    sum_dict: dict[str, int] = dict.fromkeys(readnames, 0)
    for row in cr.sv_signals:
        if row[5] == 1:
            sum_dict[row_to_annotated_readname(row=row)] += int(row[6])
    return readnames, sum_dict


def sum_ins_per_read(
    cr: datatypes.CandidateRegion, readnames: list[str]
) -> tuple[list[str], dict[str, int]]:
    sum_dict: dict[str, int] = dict.fromkeys(readnames, 0)
    for row in cr.sv_signals:
        if row[5] == 0:
            sum_dict[row_to_annotated_readname(row=row)] += int(row[6])
    return readnames, sum_dict


def n_bnds_per_read(
    cr: datatypes.CandidateRegion, readnames: list[str]
) -> tuple[list[str], dict[str, int]]:
    sum_dict: dict[str, int] = dict.fromkeys(readnames, 0)
    for row in cr.sv_signals:
        if row[5] >= 3:
            sum_dict[row_to_annotated_readname(row=row)] += 100
    return readnames, sum_dict


def cr_to_measures_matrix(cr: datatypes.CandidateRegion) -> np.ndarray:
    # returns a ndarray with columns: deletions, insertions, BNDs for each readname
    readnames = sorted(
        set([row_to_annotated_readname(row=row) for row in cr.sv_signals])
    )
    _, sums_dels = sum_dels_per_read(cr=cr, readnames=readnames)
    _, sums_inss = sum_ins_per_read(cr=cr, readnames=readnames)
    readnames, sums_bnds = n_bnds_per_read(cr=cr, readnames=readnames)
    matrix = np.array(
        [list(sums_dels.values()), list(sums_inss.values()), list(sums_bnds.values())]
    ).T
    return readnames, matrix


def ordered_differences(X: np.ndarray) -> np.ndarray:
    if len(X) < 2:
        return np.array([0.0])
    return np.diff(np.sort(X))


def evaluate_clustering(X: np.ndarray, labels: np.ndarray) -> np.ndarray:
    # returns a matrix of the form [[sum_ins,sum_dels,sum_bnds] for each cluster]
    cols = X.shape[1]
    labels = np.array(labels, dtype=int)
    clusters = sorted(set(labels))
    values = np.zeros((len(clusters), cols))
    for cluster in clusters:
        x_selected = X[labels == cluster]
        if x_selected.shape[0] < 2:
            continue
        for i in range(cols):
            value = np.max(ordered_differences(x_selected[:, i]))
            values[cluster, i] = value
    return values


def print_clustering(cr, labels):
    _, matrix = cr_to_measures_matrix(cr)
    df = pd.DataFrame(matrix, columns=["insertions", "deletions", "BNDs"])
    df["cluster"] = labels
    fig = px.scatter(
        df,
        x="deletions",
        y="insertions",
        color="cluster",
        title=f"crID: {cr.crID}, n_clusters: {len(set(labels))}",
    )
    fig.show()


def find_best_clustering(
    matrix: datatypes.CandidateRegion, max_clusters: int = 10
) -> dict:
    labels = np.zeros(len(matrix), dtype=int)
    cluster_evaluation = evaluate_clustering(X=matrix, labels=labels)
    value = sum(
        [np.max(ordered_differences(matrix[:, i])) for i in range(matrix.shape[1])]
    )
    last_clustering = {
        "value": value,
        "n_clusters": 1,
        "labels": labels,
        "cluster_evaluation": cluster_evaluation,
    }
    n_distinct_rows = len(set([tuple(row) for row in matrix]))
    if n_distinct_rows < 3 or max_clusters < 2:
        value = (
            0
            if len(matrix) == 1
            else sum(
                [
                    np.max(ordered_differences(matrix[:, i]))
                    for i in range(matrix.shape[1])
                ]
            )
        )
        n_clusters = 1
        labels = np.zeros(len(matrix), dtype=int)
        return {
            "value": value,
            "n_clusters": n_clusters,
            "labels": labels,
            "cluster_evaluation": cluster_evaluation,
        }
    max_clusters = min(max_clusters, n_distinct_rows // 3 + 1)
    for n_clusters in range(2, max_clusters + 1):
        clustering = AgglomerativeClustering(n_clusters=n_clusters).fit(matrix)
        cluster_evaluation = evaluate_clustering(X=matrix, labels=clustering.labels_)
        value: float = int(np.sum(cluster_evaluation))
        if value > last_clustering["value"]:
            return last_clustering
        last_clustering = {
            "value": value,
            "n_clusters": n_clusters,
            "labels": np.array(clustering.labels_, dtype=int),
            "cluster_evaluation": cluster_evaluation.tolist(),
        }
    return last_clustering


# =============================================================================
#  clustering
# =============================================================================


# concave distance function
def closeness_function(x, k: int = 10):
    v = abs(x)
    # raise error if v-k < 2
    maxrange = 1000
    cutoff = 0.05
    if v <= k:
        return 1.0
    if v > maxrange:
        return 0.0
    else:
        if np.log2(v - k) - cutoff <= 0.0:
            return 1.0
        else:
            return 1 / np.log2(v - k) - cutoff


def probability_of_sv_presence_with_neighbors(
    svSignal: datatypes.SVsignal,
    neighbors: np.ndarray,
    N_samples: int,
    size_factor: float = 0.2,
    radius: int = 10,
) -> float:
    # a high prob is given, if many samples are represented in the neighbors
    # and the sv is very close to its neighbors
    # and the size of the signal is very similar to the neighbors (for insertions and deletions)
    # check if neighbors is empty. If so, return 0.0
    if neighbors.shape[0] == 0:
        return 0.0
    if svSignal.sv_type == 0:  # insertion
        margin = max(5, size_factor * svSignal.size)
        mask_similar_size = abs(neighbors[:, 2] - svSignal.size) <= margin
        # aplly distance function to score close signals higher than distant signals
        # pick one signal per sample that fits best
        neighbors_selected = neighbors[mask_similar_size]
        proximities = np.array(
            [
                closeness_function(x, radius)
                for x in neighbors_selected[:, 0] - svSignal.ref_start
            ]
        )
        # pick the best matching signal for each sample
        sum_signal = 0.0
        for rn_ID in np.unique(neighbors_selected[:, 5]):
            mask_sample = neighbors_selected[:, 5] == rn_ID
            sum_signal += np.max(proximities[mask_sample])
        return sum_signal / N_samples
    if svSignal.sv_type == 1:  # deletion
        # similar to insertion, but the signal is split into two parts
        # the left and the right side of the deletion
        margin = max(5, size_factor * svSignal.size)
        mask_similar_size = abs(neighbors[:, 2] - svSignal.size) <= margin
        # aplly distance function to score close signals higher than distant signals
        # pick one signal per sample that fits best
        neighbors_selected = neighbors[mask_similar_size]
        proximities_left = np.array(
            [
                closeness_function(x, radius)
                for x in neighbors_selected[:, 0] - svSignal.ref_start
            ]
        )
        proximities_right = np.array(
            [
                closeness_function(x, radius)
                for x in neighbors_selected[:, 1] - (svSignal.ref_start + svSignal.size)
            ]
        )
        # pick the best matching signal for each sample
        sum_signal = 0.0
        for rn_ID in np.unique(neighbors_selected[:, 5]):
            mask_sample = neighbors_selected[:, 5] == rn_ID
            sum_signal += np.max(proximities_left[mask_sample]) * 0.5
            sum_signal += np.max(proximities_right[mask_sample]) * 0.5
        return sum_signal / N_samples
    if svSignal.sv_type == 3 or svSignal.sv_type == 4:  # break end
        proximities = np.array(
            [
                closeness_function(x, radius)
                for x in neighbors[:, 0] - svSignal.ref_start
            ]
        )
        # pick the best matching signal for each sample
        sum_signal = 0.0
        for rn_ID in np.unique(neighbors[:, 5]):
            mask_sample = neighbors[:, 5] == rn_ID
            sum_signal += np.max(proximities[mask_sample])
        return sum_signal / N_samples
    raise ValueError(
        "sv must be of type AlignmentInsertion, AlignmentDeletion or AlignmentBreakend."
    )


def probability_of_sv_with_similar_svs(
    svSignal: datatypes.SVsignal,
    signals: np.ndarray,
    N_samples: int,
    radius: int = 10_000,
    size_tolerance: float = 0.05,
) -> float:
    # for each sample the best matching row in signals is found (matching by type and size)
    # the probability is the number of samples that have a similar signal in the vicinity.
    # It can't be calculated for breakends
    # works only fo SVs that are greater than a given size
    # check if signals is empty. If so, return 0.0
    if signals.shape[0] == 0:
        return 0.0
    if svSignal.sv_type == 0:  # insertion
        margin = max(5, size_tolerance * svSignal.size)
        # subset signals to signals of the same type and signals that max. differ in size by margin
        mask = (signals[:, 3] == 0) & (abs(signals[:, 2] - svSignal.size) <= margin)
        candidates = signals[mask]
        # apply closeness_function to the size difference of the candidates
        similarities = np.array(
            [closeness_function(x, radius) for x in candidates[:, 2] - svSignal.size]
        )
        # pick the best matching signal for each sample
        sum_signal = 0.0
        for rn_ID in np.unique(candidates[:, 5]):
            mask_sample = candidates[:, 5] == rn_ID
            sum_signal += np.max(similarities[mask_sample])
        return sum_signal / N_samples
    if svSignal.sv_type == 1:  # deletion
        margin = max(5, size_tolerance * svSignal.size)
        # subset signals to signals of the same type and signals that max. differ in size by margin
        mask = (signals[:, 3] == 1) & (abs(signals[:, 2] - svSignal.size) <= margin)
        candidates = signals[mask]
        # apply closeness_function to the size difference of the candidates
        similarities = np.array(
            [closeness_function(x, radius) for x in candidates[:, 2] - svSignal.size]
        )
        # pick the best matching signal for each sample
        sum_signal = 0.0
        for rn_ID in np.unique(candidates[:, 5]):
            mask_sample = candidates[:, 5] == rn_ID
            sum_signal += np.max(similarities[mask_sample])
        return sum_signal / N_samples
    raise ValueError(
        "sv must be of type SVsignal and must have sv_type 0 (insertion) or 1 (deletion)."
    )


# function to calculate both scores from probability_of_sv_presence_with_neighbors and probability_of_sv_with_similar_svs
# and picks the maximum of both
def probability_of_sv(
    svSignal: datatypes.SVsignal,
    signals: np.ndarray,
    N_samples: int,
    bias_for_locality: float = 0.5,
    size_factor: float = 0.2,
    size_tolerance: float = 0.05,
    radius_close: int = 10,
    radius_neighbors: int = 1000,
    radius_far: int = 10_000,
) -> float:
    """Computes the probability of the presence of a structural variant.
    Args:
        sv (_type_): structural variant of the type AlignmentInsertion, AlignmentDeletion or AlignmentBreakend
        signals (np.ndarray): array of signals with columns: [position_start,position_end,size,type,raf_index,read_index]
        N_samples (int): number of samples in all rafs
        size_factor (float, optional): The maximum difference in size of two SV signals when matching locally. Defaults to 0.05.
        radius_close (int, optional): radius to consider SV signals in locality sensitive approximation. Defaults to 10.
        radius_far (int, optional): radius to consider SV signals in similarity sensitive approximation. Defaults to 10_000.

    Returns:
        float: probability of the presence of the SV
    """
    # check input. is signals empty? is N_samples 0?
    if signals.size == 0:
        return 0.0
    if N_samples == 0:
        return 0.0
    neighbors = sv_neighbors(
        svSignal=svSignal, radius=radius_neighbors, signals=signals
    )
    score_locality = probability_of_sv_presence_with_neighbors(
        svSignal=svSignal,
        neighbors=neighbors,
        N_samples=N_samples,
        size_factor=size_factor,
        radius=radius_close,
    )
    # can compute similarity score only for insertions and deletions
    if svSignal.sv_type == 3 or svSignal.sv_type == 4:
        return score_locality
    score_similarity = probability_of_sv_with_similar_svs(
        svSignal=svSignal,
        signals=signals,
        N_samples=N_samples,
        size_tolerance=size_tolerance,
        radius=radius_far,
    )
    # weighted mean
    return (
        bias_for_locality * score_locality
        + (1.0 - bias_for_locality) * score_similarity
    )
    # weighted maximum
    # return 2*max(bias_for_locality*score_locality,(1.0-bias_for_locality)*score_similarity)


# def sum_dels_per_read(rass:list[datatypes.ReadAlignmentSignals]) -> dict[str,int]:
#     readnames = set([ras.read_name for ras in rass])
#     sum_dict:dict[str,int] = {readname:0 for readname in readnames}
#     for ras in rass:
#         for sv in ras.SV_signals:
#             if sv.sv_type == 1:
#                 sum_dict[ras.read_name] += sv.size
#     return sum_dict

# def sum_ins_per_read(rass:list[datatypes.ReadAlignmentSignals]) -> dict[str,int]:
#     readnames = set([ras.read_name for ras in rass])
#     sum_dict:dict[str,int] = {readname:0 for readname in readnames}
#     for ras in rass:
#         for sv in ras.SV_signals:
#             if sv.sv_type == 0:
#                 sum_dict[ras.read_name] += sv.size
#     return sum_dict

# def n_bnds_per_read(rass:list[datatypes.ReadAlignmentSignals]) -> dict[str,int]:
#     readnames = set([ras.read_name for ras in rass])
#     sum_dict:dict[str,int] = {readname:0 for readname in readnames}
#     for ras in rass:
#         for sv in ras.SV_signals:
#             if sv.sv_type >= 3:
#                 sum_dict[ras.read_name] += 100
#     return sum_dict


def calc_ras_scores(rass: typing.List[datatypes.ReadAlignmentSignals]) -> np.ndarray:
    ras_scores = np.zeros(len(rass))
    for i, ras in enumerate(rass):
        for svSignal in ras.SV_signals:
            weight = svSignal.size if svSignal.sv_type < 3 else 100
            ras_scores[i] += weight
    return ras_scores


# def calc_ras_scores(
#         rass:typing.List[datatypes.ReadAlignmentSignals],
#         size_factor:float=0.05,
#         radius_close:int=10,
#         radius_neighbors:int=100,
#         radius_far:int=10_000,
#         bias:float=0.5) -> np.ndarray:
#     # calc a score for each raf. The score of a raf is the sum of probabilities of each sv signal in the raf
#     N_samples = len(set([ras.read_name for ras in rass]))
#     signals = extract_signals(rass=rass)
#     # generate regions with repeats / low complexity
#     # For each
#     ras_scores = np.zeros(len(rass))
#     for i,ras in enumerate(rass):
#         for svSignal in ras.SV_signals:
#             weight = svSignal.size if svSignal.sv_type < 3 else 100
#             ras_scores[i] += probability_of_sv(
#                 svSignal=svSignal,
#                 signals=signals,
#                 N_samples=N_samples,
#                 bias_for_locality=bias,
#                 size_factor=size_factor,
#                 radius_close=radius_close,
#                 radius_neighbors=radius_neighbors,
#                 radius_far=radius_far) * weight
#     return ras_scores

# %%


def cut_reads(
    mcrID: int, path_reads_db: Path, path_intervals_db: Path
) -> tuple[dict[str, SeqRecord], list]:
    con_intervals = sqlite3.connect(
        "file:" + str(path_intervals_db) + "?mode=ro", uri=True
    )
    cur_intervals = con_intervals.cursor()
    cur_intervals.execute(f"SELECT * FROM intervals WHERE mcrID = {mcrID}")
    intervals = np.array(cur_intervals.fetchall(), dtype=object)
    # connect intervals. Each read should have only one interval
    # if a read has multiple intervals, then take the min and max position on the read
    # sort intervals by readname
    readnames = set(intervals[:, 0])
    n_interval_readnames = len(readnames)
    # columns:
    # qname,crID,mcrID,rstart,rend,cr.refID,crStart,crENd,sampleID
    con_reads = sqlite3.connect("file:" + str(path_reads_db) + "?mode=ro", uri=True)
    cur_reads = con_reads.cursor()
    reads = {
        row[0]: row[1:]
        for row in cur_reads.execute(
            """
            SELECT * FROM reads
            WHERE readname IN ({0})""".format(
                ", ".join([f"'{readname}'" for readname in readnames])
            )
        ).fetchall()
    }
    cur_reads.close()

    # filter all readnames whose readname are not in reads
    readnames = set(reads.keys())
    if n_interval_readnames - len(readnames) > 0:
        raise ValueError(
            f"Filtered {n_interval_readnames - len(readnames)} readnames from intervals, because they are not in the reads."
        )
        # log.warning(f"Filtered {n_interval_readnames - len(readnames)} readnames from intervals, because they are not in the reads.")
    intervals_merged = []
    for readname in readnames:
        intervals_read = intervals[intervals[:, 0] == readname]
        if len(intervals_read) == 1:
            intervals_merged.append(list(intervals_read[0]))
            continue
        # if multiple intervals, then take the min and max position on the read
        intervals_merged.append(
            [
                readname,  # readname
                intervals_read[0, 1],  # crID
                intervals_read[0, 2],  # mcrID
                min(intervals_read[:, 3]),  # read start
                max(intervals_read[:, 4]),  # read end
                intervals_read[0, 5],  # refID
                min(intervals_read[:, 6]),  # ref start
                max(intervals_read[:, 7]),  # ref end
                intervals_read[0, 8],  # sampleID
                intervals_read[0, 9],  # forward
                sum(intervals_read[:, 10]),  # qlen
                min(intervals_read[:, 11]),  # qstart
                max(intervals_read[:, 12]),  # qend
                min(intervals_read[:, 13]),  # rstart
                max(intervals_read[:, 14]),
            ]  # rend
        )

    # create reads dict
    # iterate intervals and cut all reads accordingly. Save the subreads to a tmp fastq file
    cutreads = dict()
    empty_reads = set()

    for row in intervals_merged:
        readname = row[0]
        try:
            dna = reads[readname][1]
        except KeyError:
            raise ValueError(f"Readname {readname} not found in reads.")
        qual = reads[readname][2]
        start, end = row[3], row[4]
        sampleID = row[8]
        forward = bool(row[9])
        ext_readname = f"{readname}.{mcrID}"
        rec = SeqRecord(
            id=ext_readname,
            name=ext_readname,
            description=f"mcrID:{mcrID},sampleID:{sampleID}",
            seq=Seq(dna[start:end]),
            letter_annotations={"phred_quality": qual[start:end]} if qual else None,
        )
        # check if read is empty:
        if len(rec.seq) == 0:
            empty_reads.add(readname)
            continue
        cutreads[ext_readname] = rec
    # filter intervals with empty reads
    intervals_merged = [row for row in intervals_merged if row[0] not in empty_reads]
    return cutreads, intervals_merged


# def score_breakend(breakend:datatypes.AlignmentBreakend,
#                    ignored_length:int=30,
#                    max_size:int=10_000,
#                    factor:float=10.0) -> float:
#     val = float(min(max(ignored_length,breakend.clippedLength),max_size))
#     # transform value in a concave function?
#     return float(np.log2(val))*factor


# def score_raf_svs(
#             a:datatypes.ReadAlignmentFragment,
#             min_signal_size:int=30,
#             min_bnd_size:int=100,
#             margin_size:int=100,
#             max_tail_size:int=10_000,
#             consensus_length:int=0,
#             bnd_factor:float=10.0) -> float:
#     """computes the sum of the maximum number of bases of dels, ins and clipped tail"""
#     sum_del = sum([deletion.ref_length if deletion.ref_length >= min_signal_size else 0 for deletion in a.deletions]) if len(a.deletions) > 0 else 0
#     sum_ins = sum([insertion.read_length if insertion.read_length >= min_signal_size else 0 for insertion in a.insertions]) if len(a.insertions) > 0 else 0
#     bnd_scores = [0.0,0.0]
#     if a.breakendLeft and a.breakendLeft.referencePosition >= margin_size:
#         bnd_scores[0]=score_breakend(breakend=a.breakend_left,ignored_length=min_bnd_size,max_size=max_tail_size,factor=bnd_factor)
#     if a.breakendRight and a.breakendRight.referencePosition <= consensus_length - margin_size:
#         bnd_scores[1]=score_breakend(breakend=a.breakend_right,ignored_length=min_bnd_size,max_size=max_tail_size,factor=bnd_factor)
#     return float(sum_del + sum_ins + bnd_scores[0] + bnd_scores[1])


def parse_ReadAlignmentSignals_from_alignment(
    alignment: pysam.AlignedSegment,
    sampleID: int,
    min_signal_size: int,
    min_bnd_size: int,
) -> datatypes.ReadAlignmentSignals:
    ref_start, ref_end, read_start, read_end = alignments_to_rafs.get_start_end(
        alignment
    )
    sv_signals = alignments_to_rafs.parse_SVsignals_from_alignment(
        alignment=alignment,
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=read_start,
        read_end=read_end,
        min_signal_size=min_signal_size,
        min_bnd_size=min_bnd_size,
    )
    return datatypes.ReadAlignmentSignals(
        sampleID=int(sampleID),
        read_name=str(alignment.query_name),
        reference_name=str(alignment.reference_name),
        alignment_forward=bool(not alignment.is_reverse),
        SV_signals=sorted(sv_signals),
    )


# pick first read DNA sequence and write to tmp fasta file
# index fasta file with samtools index
# write all reads to tmp fastq file
# align all other reads to first read with minimap2 (util.align_reads_with_minimap())
# load as rafs (alignments_to_rafs.alignment_to_alignment_segment)
# score rafs (score_raf_svs)
def alignments_reads_to_longest(
    longest: SeqRecord,
    path_reads: Path,
    path_alignments: Path,
    threads: int,
    tmp_dir_path: Path | None = None,
) -> None:
    """Takes a SeqRecord and a path to a fastq file of reads and returns a dict of {read_name: score}
    The score is the sum of the maximum number of bases of dels, ins and clipped tail"""
    # keep_tmps_path = Path(keep_tmps_path)
    # if keep_tmp_path is != Path('') and if it does not exist, then create it
    # if keep_tmps_path != Path('') and not keep_tmps_path.exists():
    #     keep_tmps_path.mkdir(parents=True,exist_ok=True)
    # write longest read to tmp fasta file
    tmp_longest = tempfile.NamedTemporaryFile(
        suffix=".longest.fasta",
        dir=tmp_dir_path,
        delete=False if tmp_dir_path else True,
    )
    log.info("writing longest read to tmp fasta file")
    with open(tmp_longest.name, "w") as f:
        SeqIO.write(longest, f, "fasta")
    # index fasta file with samtools index
    subprocess.run(split(f"samtools faidx {tmp_longest.name}"), check=True)
    log.info("aligning reads to longest read")
    util.align_reads_with_minimap(
        reads=path_reads,
        reference=tmp_longest.name,
        bamout=path_alignments,
        tech="map-ont",
        aln_args=" --secondary=no -Y --sam-hit-only -U 10,500 -H",
        threads=threads,
    )
    # remove fasta index again
    return


def get_handicap(
    rass: typing.List[datatypes.ReadAlignmentSignals], handicap_cutoff: float = 0.9
) -> float:
    # adjust the scores if the per sample summed ins or summed dels have a high variance
    # (variance in the sense of difference between two neighbors in a sorted list)
    ins_sums = np.array(
        [sum([i.size for i in ras.SV_signals if i.sv_type == 0]) for ras in rass]
    )
    del_sums = np.array(
        [sum([d.size for d in ras.SV_signals if d.sv_type == 1]) for ras in rass]
    )
    sum_indel_bp = np.sort(ins_sums + del_sums)
    # print('\n',ins_sums,del_sums,'\n',sorted(np.diff(sum_indel_bp)),'\n')
    # calculate the handicap
    # its the mean value of the handicap_cutoff% lowest values of the diff of the sorted sum_indel_bp
    # instead, take the mean of the 3 reads with the least sum of ins and dels
    np.mean(sorted(np.diff(sum_indel_bp))[: int(handicap_cutoff * len(sum_indel_bp))])
    handicap = np.mean(
        sorted(np.diff(sum_indel_bp))[: int(handicap_cutoff * len(sum_indel_bp))]
    )
    # adjust the threshold accoring to some measure on the diff of the sorted sum_indel_bp
    return handicap


def score_ras_from_alignments(
    alignments: Path, min_signal_size: int, min_bnd_size: int
) -> tuple[dict[str, float], float]:
    # def parse_ReadAlignmentFragment_from_alignment(
    #         alignment:pysam.AlignedSegment,
    #         sampleID:int,
    #         min_signal_size:int,
    #         min_bnd_size:int
    #         ) -> datatypes.ReadAlignmentFragment:

    # score alignments
    # get rafs
    log.info("converting alignments to rafs")
    rass = [
        parse_ReadAlignmentSignals_from_alignment(
            alignment=a,
            sampleID=a.query_name.split(".")[1],
            min_bnd_size=min_bnd_size,
            min_signal_size=min_signal_size,
        )
        for a in pysam.AlignmentFile(alignments)
    ]

    # score rafs (read_name might occur multiple times)
    # re-write raf_scores. It should first create statistics over all sv signal
    # then generate weights for each sv signal, based on a maximum score
    # of 1) a distance score and 2) a size score.
    log.info("scoring rafs")

    handicap = get_handicap(rass=rass, handicap_cutoff=0.9)

    # return dict of {read_name: score}
    # always pick the highest scoring raf for each read_name
    rass_scores = calc_ras_scores(rass=rass)
    rafs_scores_dict = defaultdict(float)
    for ras, score in zip(rass, rass_scores):
        if ras.read_name in rafs_scores_dict:
            if rafs_scores_dict[ras.read_name] < score:
                rafs_scores_dict[ras.read_name] = score
        else:
            rafs_scores_dict[ras.read_name] = score
    log.debug(
        f"done scoring rafs with {len(set([ras.read_name for ras in rass]))} reads and {len(rass)} ReadAlignmentSignals objects."
    )
    log.debug(f"Resulting: {rafs_scores_dict}.")
    log.debug(f"handicap: {handicap}")
    return rafs_scores_dict, handicap


def randomsample_bam_to_sam(
    input_bamfile: Path, output_samfile: Path, n_samples: int
) -> list[str]:
    # read alignments from bam file
    alignments = list(pysam.AlignmentFile(input_bamfile))
    readnames = set([a.query_name for a in alignments])
    # drop all alignments that have the same name as the reference
    alignments = [a for a in alignments if a.query_name != a.reference_name]
    # random sample n_samples from alignments if n_samples < len(alignments)
    if n_samples < len(readnames):
        readnames = set(randomsample(list(readnames), n_samples))
    # filter all alignments that have a readname in readnames
    alignments = [a for a in alignments if a.query_name in readnames]
    # write to samfile
    with pysam.AlignmentFile(
        output_samfile, "w", template=pysam.AlignmentFile(input_bamfile, "rb")
    ) as f:
        for a in alignments:
            f.write(a)
    return readnames


# for loci with few repeated subsequences
def make_consensus_with_racon(
    longest_read: SeqRecord,
    alignments: Path,
    chosen_reads: list[SeqRecord],
    name: str,
    output: Path,
    threads: int,
    maxreads: int = 20,
    tmp_dir_path: Path | None = None,
) -> int:
    log.info(f"making consensus with racon for {name} with {len(chosen_reads)} reads")
    if longest_read.id not in set([r.id for r in chosen_reads]):
        raise ValueError("chosen_reads must contain the longest read.")
    if len(chosen_reads) == 1:
        # can't make a consensus from one read. instead,just write the one read sequence to the output.
        log.warning("only one read present in chosen_reads. Writing it to the output")
        output_record = deepcopy(longest_read)
        output_record.id = name
        output_record.name = name
        output_record.description = f"cluster_n: {len(chosen_reads)}"
        with open(output, "w") as f:
            SeqIO.write(output_record, f, "fasta")
        return len(longest_read.seq)

    tmp_reads_sam = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path,
        suffix=f".longest.{longest_read.id}.sam",
        delete=False if tmp_dir_path else True,
    )

    sampled_reads = randomsample_bam_to_sam(
        input_bamfile=alignments, output_samfile=tmp_reads_sam.name, n_samples=maxreads
    )

    tmp_reads_fasta = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, suffix=".reads.fasta", delete=False if tmp_dir_path else True
    )
    with open(tmp_reads_fasta.name, "w") as f:
        for read in chosen_reads:
            if read.id in sampled_reads:
                SeqIO.write(read, f, "fasta")

    tmp_longest_fasta = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path,
        suffix=".longest.fasta",
        delete=False if tmp_dir_path else True,
    )
    with open(tmp_longest_fasta.name, "w") as f:
        SeqIO.write(longest_read, f, "fasta")

    tmp_output = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path,
        suffix=".consensus.fasta",
        delete=False if tmp_dir_path else True,
    )

    try:
        cmd_racon = f"racon -t {threads} {str(tmp_reads_fasta.name)} {str(tmp_reads_sam.name)} {str(tmp_longest_fasta.name)}"
        with open(tmp_output.name, "w") as f:
            subprocess.check_call(split(cmd_racon), stdout=f)
        # check if output has content
        n_lines = len(open(tmp_output.name).readlines())
        if n_lines == 0:
            log.warning(
                f"racon could not generate consensus from {len(chosen_reads)} reads and uses the input longest read as output"
            )
            with open(tmp_output.name, "w") as f:
                longest_read.description = "cluster_n: 1"
                SeqIO.write(longest_read, f, "fasta")
    except:
        # just write longest to tmp output
        log.warning(
            f"racon could not generate consensus from {len(chosen_reads)} reads and uses the input longest read as output"
        )
        with open(tmp_output.name, "w") as f:
            SeqIO.write(longest_read, f, "fasta")
    with open(tmp_output.name, "r") as f:
        try:
            contig = SeqIO.read(f, "fasta")
            contig.id = name
            contig.name = name
            contig.description = f"cluster_n: {len(chosen_reads)}"
        except:
            log.warning(f"No sequence found in {str(tmp_output.name)}")
    # write contig to output
    with open(output, "w") as f:
        SeqIO.write(contig, f, "fasta")
    log.info(
        f"consensus for {name} created with racon and written to file {str(output)}"
    )
    # read the consensus sequence and return its length
    return len(contig.seq)


# in case a locus is hard to assemble
def make_consensus_with_lamassemble(
    chosen_reads: list[SeqRecord],
    output: Path,
    name: str,
    min_coverage: float,
    threads: int,
    lamassemble_mat_path: Path,
    format="fastq",
):
    # chosen_reads is of the form [SeqRecord]
    tmp_dir = tempfile.TemporaryDirectory()
    path_tmp_dir = Path(tmp_dir.name)
    # first write reads to tmp fastq file
    path_reads = path_tmp_dir / f"reads.{format}"
    SeqIO.write(chosen_reads, path_reads, format)
    cmd_consensus = split(
        f"lamassemble \
        -n {name} \
        -s {str(min_coverage)}% \
        -P {str(threads)} \
        {str(lamassemble_mat_path)} \
        {str(path_reads)}"
    )
    consensus_file_handle = open(output, mode="w")
    subprocess.check_call(cmd_consensus, stdout=consensus_file_handle)


# def merge_consensus_sequences(
#             name:str,
#             reads:dict,
#             selection:np.ndarray,
#             consensus_meta:dict,
#             path_consensus:Path,
#             lamassemble_mat_path:Path,
#             min_coverage:float=0.1,
#             threads:int=1,
#             method='racon') -> int:
#     # merges consensus sequences of a selection of clusters
#     # selection is a np.ndarray of cluster IDs that are merged
#     # reads is of form {id:SeqRecord}
#     # consensus_meta is of the form:
#     # {ind: [representative,nreads,[readIDs],path_consensus,length,scores_against_other_consensuses]}
#     # all ind are given in selection
#     # ---
#     # assemble with lamassemble: make_consensus_with_lamassemble
#     # create list of all reads
#     if method not in ['lamassemble','racon']:
#         raise ValueError("method must be either 'lamassemble' or 'racon'.")
#     merged_reads = []
#     for ind in selection:
#         merged_reads += [reads[readID] for readID in consensus_meta[ind][2]]
#     if len(merged_reads) == 1:
#         consensus = SeqRecord(
#             seq=merged_reads[0].seq,
#             id=name,
#             description="")
#         # write this to path_consensus
#         with open(path_consensus, "w") as f:
#             SeqIO.write(consensus, f, "fasta")
#             return 1
#     if method == 'lamassemble':
#         make_consensus_with_lamassemble(
#             chosen_reads=merged_reads,
#             name=name,
#             output=path_consensus,
#             min_coverage=min_coverage,
#             threads=threads,
#             lamassemble_mat_path=lamassemble_mat_path)
#     if method == 'racon':
#         merge_consensuses_with_racon(
#             name=name,
#             consensus_meta=consensus_meta,
#             path_consensus=path_consensus,
#             reads=reads,
#             selection=selection,
#             threads=threads
#         )
#     return len(merged_reads)

# def make_consensus(
#             name:str,
#             longest_read_id:str,
#             reads:dict,
#             partition_readIDs:list[str],
#             path_consensus:Path,
#             lamassemble_mat_path:Path,
#             min_coverage:float=0.1,
#             threads:int=1,
#             min_consensus_size:int=500,
#             method:str="racon") -> typing.Tuple[Path,int]:
#     # dict_partitions is of form {id_logest:[ids]}
#     # reads is of form {id:SeqRecord}
#     """creates an indexed consensus.{mcrID}.fasta in the provided tmp_dir_path with lamassemble"""
#     # if only one read (the longest one) is passed, then just write longest read seq to path_consensus with name name
#     if method not in ['lamassemble','racon']:
#         raise ValueError("method must be either 'lamassemble' or 'racon'.")
#     if len(partition_readIDs) == 0:
#         with open(path_consensus, "w") as f:
#             # create new record with name name
#             consensus_record = SeqRecord(reads[longest_read_id].seq, id=name, description="singleton consensus from read "+longest_read_id)
#             SeqIO.write(consensus_record, f, "fasta")
#         cmd_index = split(f"samtools faidx {str(path_consensus)}")
#         subprocess.check_call(cmd_index)
#         consensus_length = len(reads[longest_read_id].seq)
#         return consensus_length

#     # if more than one read is passed, then create consensus sequence with chosen tool, e.g. lamassemble
#     if method == "lamassemble":
#         make_consensus_with_lamassemble(
#             name=name,
#             output=path_consensus,
#             chosen_reads=[reads[k] for k in [*partition_readIDs,longest_read_id]],
#             threads=threads,
#             min_coverage=min_coverage,
#             lamassemble_mat_path=lamassemble_mat_path)

#     if method == "racon":
#         chosen_reads = [reads[k] for k in partition_readIDs]
#         make_consensus_with_racon(
#             longest_read=reads[longest_read_id],
#             alignments=alignments,
#             chosen_reads=chosen_reads,

#             name=name,
#             output=path_consensus,
#             threads=threads)


# # read consensus_path to list of lines
# consensus = SeqIO.read(path_consensus, "fasta")
# consensus_length = len(consensus.seq)
# if consensus_length < min_consensus_size:
#     log.info(f"Consensus length is less than {min_consensus_size} bp. Returning 0.")
#     return 0
# cmd_index = split(f"samtools faidx {str(path_consensus)}")
# subprocess.check_call(cmd_index)
# consensus_length = int(open(str(path_consensus)+'.fai').readline().split("\t")[1])
# return consensus_length

# %%


def sort_reads_by_size(reads, intervals) -> typing.Dict[str, SeqRecord]:
    # coverage_sizes is of the form readID:(ref_coverage_size,read_covera_size)
    coverage_sizes = {
        f"{row[0]}.{row[2]}": (int(row[4]) - int(row[3]), int(row[12]) - int(row[11]))
        for row in intervals
    }
    # sort by coverage size
    sorted_read_ids = {
        k: reads[k]
        for k, v in sorted(
            coverage_sizes.items(), key=lambda item: item[1][0], reverse=True
        )
    }
    # sort reads according to sorted_read_ids
    return {k: reads[k] for k in sorted_read_ids}


# TODO: consider informative parts of the reads. Scoring can only be achieved if the aligned reads cover the regions
# of the longest read. If the longest read is not covered by the aligned reads, then the score is the length of the aligned read.
def create_read_partitions(
    pool: dict[str, int],
    reads: dict[str, SeqRecord],
    reads_fastq: Path,
    basepath_alignments: Path,
    min_signal_size: int,
    min_bnd_size: int,
    reads_partitioning_threshold: float,
    threads: int = 2,
    clusterID_start: int = 0,
    use_handicap: bool = True,
    tmp_dir_path: Path | None = None,
) -> tuple[dict[str, list[str]], dict[str, Path]]:
    """use_handicap adds a tolerance to the maximum pairwise difference scores that are used to partition the reads. set use_handicap=False\
to have a very strict cutoff."""
    dict_partitions = defaultdict(list)
    log.info(
        f"partitioning reads by alignment to longest read of remaining {len(pool)} reads"
    )
    # pick reads with score < threshold_init_alignment_signal
    all_dict_scores: dict[str, dict[str, float]] = {}
    alignment_paths: dict[str, Path] = (
        {}
    )  # save the paths of the read-to-longest alignments
    # -------------------------------------------------------------------------
    #  partitioning reads: align reads to ongest read in pool, then remove all reads with score < threshold
    #  from pool and add them to dict_partitions
    iteration = -1
    clusterID = clusterID_start
    while len(pool) > 1:
        iteration += 1
        longest_read: SeqRecord = reads[max(pool, key=pool.get)]
        alignment_paths[longest_read.id] = (
            f"{basepath_alignments}.cluster-{clusterID}.bam"
        )
        clusterID += 1
        alignments_reads_to_longest(
            longest=longest_read,
            path_alignments=alignment_paths[longest_read.id],
            path_reads=reads_fastq,
            threads=threads,
            tmp_dir_path=tmp_dir_path,
        )
        dict_scores, handicap = score_ras_from_alignments(
            alignments=alignment_paths[longest_read.id],
            min_bnd_size=min_bnd_size,
            min_signal_size=min_signal_size,
        )
        log.debug(
            f"iteration {iteration}: {len(dict_scores)} reads scored. dict_scores: {dict_scores}"
        )
        all_dict_scores[longest_read.id] = dict_scores
        picked_reads = {
            read
            for read, score in dict_scores.items()
            if score < reads_partitioning_threshold + handicap * int(use_handicap)
        }
        if longest_read.id not in picked_reads:
            picked_reads.add(longest_read.id)
        # update pool by removing picked reads
        pool = {
            read: length for read, length in pool.items() if read not in picked_reads
        }
        # update dict_partitions by adding picked reads
        dict_partitions[longest_read.id] = list(picked_reads)
    # if one is remaining in pool, add it to dict_partitions
    if len(pool) == 1:
        log.info("pool has one remaining read. Adding it to dict_partitions")
        longest_read = reads[list(pool.keys())[0]]
        dict_partitions[longest_read.id] = [longest_read.id]
        alignment_paths[longest_read.id] = (
            f"{basepath_alignments}.cluster-{clusterID}.bam"  # careful, this file never exists!
        )
        clusterID += 1
    dict_partitions = {
        readname: sorted(partition) for readname, partition in dict_partitions.items()
    }
    log.info("done partitioning reads")
    return dict_partitions, alignment_paths


def create_read_partitions_with_initial_signal_clusters(
    crs: list[datatypes.CandidateRegion],
    pool: dict[str, int],
    reads: dict[str, SeqRecord],
    basepath_alignments: Path,
    min_signal_size: int,
    min_bnd_size: int,
    use_handicap: bool,
    reads_partitioning_threshold: float,
    threads: int = 2,
    tmp_dir_path: Path | None = None,
) -> tuple[dict[str, list[str]], dict[str, Path]]:
    measures_results: list[tuple[list[str], np.ndarray]] = [
        cr_to_measures_matrix(cr) for cr in crs
    ]
    all_readnames = np.concatenate([mr[0] for mr in measures_results], axis=0)
    all_matrix = np.concatenate([mr[1] for mr in measures_results], axis=0)
    # if the signals are dominated by INS or DELS, attempt a clustering with create_read_partitions_with_inital_signal_clusters
    n_samples = len(set([sv_signal[9] for cr in crs for sv_signal in cr.sv_signals]))
    max_clusters = 4 + int(np.ceil(np.log2(n_samples) * n_samples / 2))
    best_clustering: dict = find_best_clustering(
        matrix=all_matrix, max_clusters=max_clusters
    )
    clusterID_start = 0
    collected_results = []
    for clusterID in set(best_clustering["labels"]):
        # create a pool for each cluster
        labels = best_clustering["labels"] == clusterID
        readnames = all_readnames[labels].tolist()
        sub_pool = {
            readname: pool[readname] for readname in readnames if readname in pool
        }
        # check if pool is empty. If so, continue
        if len(sub_pool) < 1:
            continue
        adjusted_reads_partitioning_threshold = max(
            reads_partitioning_threshold,
            sum(best_clustering["cluster_evaluation"][clusterID]),
        )
        #
        tmp_reads_fastq = tempfile.NamedTemporaryFile(
            dir=tmp_dir_path,
            suffix=f".reads.cluster-{clusterID}.fastq",
            delete=False if tmp_dir_path else True,
        )
        # write reads to tmp fastq file
        reads_cluster = {
            readname: read for readname, read in reads.items() if readname in sub_pool
        }
        SeqIO.write(
            sequences=reads_cluster.values(),
            handle=open(tmp_reads_fastq.name, "w"),
            format="fastq",
        )

        dict_partitions, alignment_paths = create_read_partitions(
            pool=sub_pool,
            reads=reads_cluster,
            reads_fastq=Path(tmp_reads_fastq.name),
            basepath_alignments=basepath_alignments,
            min_signal_size=min_signal_size,
            min_bnd_size=min_bnd_size,
            reads_partitioning_threshold=adjusted_reads_partitioning_threshold,
            threads=threads,
            use_handicap=use_handicap,
            clusterID_start=clusterID_start,
            tmp_dir_path=tmp_dir_path,
        )
        collected_results.append((dict_partitions, alignment_paths))
        clusterID_start += len(dict_partitions)
    # merge collected results to dict_partitions and alignment_paths
    all_dict_partitions = {
        readname: sorted(partition)
        for result in collected_results
        for readname, partition in result[0].items()
    }
    all_alignment_paths = {
        readname: path
        for result in collected_results
        for readname, path in result[1].items()
    }
    return all_dict_partitions, all_alignment_paths



# %%
# =============================================================================


# write function to get the minimum and maximum extends from intervals
# for each chromosome, find the minimum start and maximum end of all intervals
def find_regions_from_intervals(
    intervals: list,
) -> typing.Dict[str, typing.Tuple[int, int]]:
    regions = dict()
    np_intervals = np.array([row[1:] for row in intervals], dtype=int)
    for refID in np.unique(np_intervals[:, 4]):
        # get all intervals of refID
        intervals_ref = np_intervals[np_intervals[:, 4] == refID]
        # get minimum start and maximum end
        region = intervals_ref[:, 5].min(), intervals_ref[:, 6].max()
        regions[refID] = region
    return {(k, *interval) for k, interval in regions.items()}


# write a function to find a minimal set of cut reads that cover as much of the
# regions as possible


# assuming reads have been sorted by size
def write_longest_read(
    path_output_fasta, reads, mcrID, threads, path_output_bam, start_time, logfile
):
    logzero.logfile(logfile, loglevel=logzero.DEBUG)
    reads_keys: typing.Dict[str,] = list(reads.keys())
    if len(reads_keys) == 0:
        raise (
            f"No reads found for mcrID {mcrID} and output bam file {path_output_bam}."
        )
    keys_iter = iter(reads_keys)
    with open(path_output_fasta, "w") as f:
        readID = next(keys_iter)
        read = reads[readID]
        read.description = "failed singleton consensus from read " + read.id
        read.id = f"consensus.{str(mcrID)}.0.failed"
        read.name = f"consensus.{str(mcrID)}.0.failed"
        SeqIO.write(read, f, "fasta")
        # and align the read to itself
    util.align_reads_with_minimap(
        reads=path_output_fasta,
        reference=path_output_fasta,
        bamout=path_output_bam,
        threads=threads,
        aln_args="--secondary=no --sam-hit-only -U 10,500",
        tech="map-ont",
    )


def consensus(
    mcrID: int,
    path_crs: Path,
    path_reads_db: Path,
    path_intervals_db: Path,
    path_output_fasta: Path,
    path_output_bam: Path,
    min_signal_size: int,
    min_bnd_size: int,
    reads_partitioning_threshold: float,
    threads: int,
    use_handicap: bool,
    use_preclustering: bool,
    min_cluster_size: int,
    reads_min_size: int = 500,
    logfile: None | Path = None,
    tmp_dir_path: Path | None = None,
) -> None:
    if tmp_dir_path:
        tmp_dir_path = Path(tmp_dir_path)
    # -------------------------------------------------------------------------
    tmp_wdir = tempfile.TemporaryDirectory()
    tmp_wdir_path = Path(tmp_wdir.name)
    # -------------------------------------------------------------------------
    logzero.logfile(logfile, loglevel=logzero.DEBUG, mode="w")
    start_time = time.time()
    log.info("cut reads")

    reads, intervals = cut_reads(
        mcrID=mcrID, path_reads_db=path_reads_db, path_intervals_db=path_intervals_db
    )
    # sampleIDs:list[int] = sorted(list(set([i[8] for i in intervals])))

    # -------------------------------------------------------------------------
    # find correct covered region of read
    log.info("sort reads")
    reads = sort_reads_by_size(reads=reads, intervals=intervals)

    # filter reads that are shorter than 500 bases
    reads_filtered = {k: v for k, v in reads.items() if len(v.seq) >= reads_min_size}
    reads_lengths = [str(len(v.seq)) for _, v in reads_filtered.items()]
    log.info(f"all lengths of filtered cut reads to be used:{','.join(reads_lengths)}")

    # check if all reads were filtered. if so, raise an error
    if len(reads_filtered) == 0:
        write_longest_read(
            path_output_fasta=path_output_fasta,
            reads=reads,
            mcrID=mcrID,
            threads=threads,
            path_output_bam=path_output_bam,
            start_time=start_time,
            logfile=logfile,
        )
        log.warn(
            f"no reads longer than {reads_min_size} bases could be found for mcrID {mcrID}."
        )
        return
    else:
        reads = reads_filtered

    # pool is a dict of the form {read_name:length}
    # the key must be changed. reads should be sorted according to the greatest group of other reads that are similar to it
    log.info("creating pool")
    # create pool and drop any cut reads that are shorter than 500 bases
    pool = {r.id: len(r.seq) for r in reads.values() if len(r.seq) >= reads_min_size}
    # if show_alignments:
    #     # print all chosen reads and their lengths
    #     print("pool:")
    #     for read in pool.keys():
    #         print(pool[read],read)
    # if the pool is empty, then just write the longest read to output
    if len(pool) == 0:
        write_longest_read(
            path_output_fasta=path_output_fasta,
            reads=reads,
            mcrID=mcrID,
            threads=threads,
            path_output_bam=path_output_bam,
            start_time=start_time,
            logfile=logfile,
        )
        return
    # log_data["preparations_time"] = time.time()-start_time
    # -------------------------------------------------------------------------
    #  read partitioning
    # each partition is written to this dict in the form {representative_read_id: [read_0, read_1, ...]}
    # create partitions of reads until there are no more partitions than samples+1
    regions = find_regions_from_intervals(intervals)
    crs = load_crs_per_regions(regions=regions, path_crs=path_crs)
    time_read_partitioning = time.time()
    dict_partitions = None
    basepath_alignments = tmp_wdir_path / "reads-to-longest-alignments"

    # write tmp fastq file of sorted reads
    reads_fastq = tempfile.NamedTemporaryFile(
        prefix="reads.",
        suffix=f".{mcrID}.fastq",
        delete=False if tmp_dir_path else True,
        dir=tmp_dir_path,
    )
    SeqIO.write(reads.values(), reads_fastq.name, "fastq")

    if use_preclustering:
        dict_partitions, alignment_paths = (
            create_read_partitions_with_initial_signal_clusters(
                crs=crs,
                pool=pool,
                reads=reads,
                basepath_alignments=basepath_alignments,
                min_bnd_size=min_bnd_size,
                min_signal_size=min_signal_size,
                reads_partitioning_threshold=reads_partitioning_threshold,
                tmp_dir_path=tmp_dir_path,
                threads=threads,
                use_handicap=use_handicap,
            )
        )
    else:
        dict_partitions, alignment_paths = create_read_partitions(
            pool=pool,
            reads=reads,
            basepath_alignments=basepath_alignments,
            min_bnd_size=min_bnd_size,
            reads_fastq=reads_fastq.name,
            min_signal_size=min_signal_size,
            reads_partitioning_threshold=reads_partitioning_threshold,
            tmp_dir_path=tmp_dir_path,
            threads=threads,
            use_handicap=use_handicap,
            clusterID_start=0,
        )
    # k = longestID
    consensus_meta: dict[int, list] = {
        i: [k, len(read_ids), read_ids]
        for i, (k, read_ids) in enumerate(dict_partitions.items())
    }
    n_partitions = len(consensus_meta)

    if dict_partitions is None:
        raise ValueError(
            f"no partitions could be created for mcrID {mcrID}. Pool={str(pool)}; consensus_meta={str(consensus_meta)}"
        )
    time_read_partitioning = time.time() - time_read_partitioning
    # log_data["read_partitioning_time"] = time.time()-time_read_partitioning

    # -------------------------------------------------------------------------
    #  create consensus sequences of all partitions
    log.info("creating consensus sequences of all partitions")
    time_consensus_creation = time.time()
    for i, (r0_id, read_ids) in enumerate(dict_partitions.items()):
        # write reads of partition to tmp fastq file
        tmp_reads_of_pool = tempfile.NamedTemporaryFile(
            prefix="pooled_reads.",
            suffix=f".cluster-{i}.fastq",
            delete=False if tmp_dir_path else True,
            dir=tmp_dir_path,
        )
        with open(tmp_reads_of_pool.name, "w") as f:
            SeqIO.write([reads[read_id] for read_id in read_ids], f, "fastq")
        path_consensus = tmp_wdir_path / f"consensus.{i}.fasta"
        consensus_meta[i].append(path_consensus)

        chosen_reads = [reads[k] for k in dict_partitions[r0_id]]
        consensus_length = make_consensus_with_racon(
            name=f"consensus.{mcrID}.{i}",
            alignments=alignment_paths[r0_id],
            chosen_reads=chosen_reads,
            longest_read=reads[r0_id],
            output=path_consensus,
            maxreads=20,
            threads=threads,
            tmp_dir_path=tmp_dir_path,
        )
        # get length of consensus sequence by reading the fasta file with Biopython
        # consensus_length = len(next(SeqIO.parse(path_consensus, "fasta")).seq)
        consensus_meta[i].append(consensus_length)
    time_consensus_creation = time.time() - time_consensus_creation
    # log_data["consensus_creation_time"] = time.time()-time_consensus_creation

    # filter consensus_meta. delete any consensus with length < reads_min_size
    consensus_meta = {k: v for k, v in consensus_meta.items() if v[4] >= reads_min_size}

    # [consensus_meta[i][3] for i in consensus_meta.keys()]
    # if consensus meta is empty, then just write the longest read to output
    if len(consensus_meta) == 0:
        log.warning(f"no consensus sequence could be created for mcrID {mcrID}.")
        write_longest_read(
            path_output_fasta=path_output_fasta,
            reads=reads,
            mcrID=mcrID,
            threads=threads,
            path_output_bam=path_output_bam,
            start_time=start_time,
            logfile=logfile,
        )
        return
    # -------------------------------------------------------------------------
    #  cat the created consensus sequences to one file
    cmd_cat = ["cat", *[consensus_meta[i][3] for i in consensus_meta.keys()]]
    tmp_all_consensus = tempfile.NamedTemporaryFile(
        prefix="all_consensus.",
        suffix=".fasta",
        delete=False if tmp_dir_path else True,
        dir=tmp_dir_path,
    )
    with open(tmp_all_consensus.name, "w") as f:
        subprocess.check_call(cmd_cat, stdout=f)

    new_consensuses = [consensus_meta[i][3] for i in range(len(consensus_meta))]
    # check every consensus fasta file if the header is present.
    # the first line should start with ">consensus"
    for fp in new_consensuses:
        with open(fp, "r") as f:
            if not f.readline().startswith(">consensus"):
                raise ValueError(
                    f"""consensus fasta file {str(fp)} does not have the correct header. mcrID = {mcrID}.
new_consensuses = {str(new_consensuses)}"""
                )

    # -------------------------------------------------------------------------
    #  all cut reads to consensus alignments go to: f"consensus/{dir_id}/consensus.{mcrID}.bam"
    log.info("aligning all reads of candidate region to all final consensus sequences")
    time_align_reads = time.time()
    # cat all consensus sequences and write to output file
    # align reads to consensuses
    cmd_cat = ["cat", *[str(path) for path in new_consensuses]]
    with open(path_output_fasta, "w") as f:
        subprocess.check_call(cmd_cat, stdout=f)
    cmd_index_fasta = split(f"samtools faidx {str(path_output_fasta)}")
    subprocess.check_call(cmd_index_fasta)
    tmp_cutread_to_consensus_bam = tempfile.NamedTemporaryFile(
        prefix="cutreads_to_consensus.",
        suffix=".bam",
        delete=False if tmp_dir_path else True,
        dir=tmp_dir_path,
    )

    util.align_reads_with_minimap(
        reads=Path(reads_fastq.name),
        reference=path_output_fasta,
        bamout=tmp_cutread_to_consensus_bam.name,
        tech="map-ont",
        threads=threads,
        aln_args=" --secondary=no --sam-hit-only -U 10,500",
    )
    # filter alignments to keep only primary alignments
    cmd_filter = f"samtools view -b -F 256 {str(tmp_cutread_to_consensus_bam.name)}"
    with open(path_output_bam, "w") as f:
        subprocess.check_call(split(cmd_filter), stdout=f)
    # index bam file
    cmd_index = split(f"samtools index {str(path_output_bam)}")
    subprocess.check_call(cmd_index)
    # print elapsed time, formatted like hh:mm:ss
    time_align_reads = time.time() - time_align_reads
    elapsed_time = time.time() - start_time

    # ===================================================================================================================
    #   not working properly right now
    n_final_consensuses = len(new_consensuses)
    # log.info(f"filter the consensus sequences to only those that are present in the cut reads to consensus alignments")
    # # filter all consensus sequences that are not present as a reference in the cut reads to consensus alignments
    # consensus_sequences_present_in_cutreads_alignments:dict[str,int] = dict()
    # with pysam.AlignmentFile(path_output_bam, "rb") as bamfile:
    #     consensus_sequences_present_in_cutreads_alignments = {refname:0 for refname in bamfile.references}
    #     for aln in bamfile.fetch():
    #         consensus_sequences_present_in_cutreads_alignments[aln.reference_name] += 1
    # # write a tmp file with all consensus sequences that are present in the cut reads to consensus alignments
    # # and write all to the fasta output file
    # tmp_filtered_consensus_output = tempfile.NamedTemporaryFile(prefix='filtered_consensus.',suffix=".fasta",delete=False if tmp_dir_path else True,dir=tmp_dir_path)
    # n_final_consensuses = 0
    # assert min_cluster_size > 0, "min_cluster_size must be greater than 0."
    # with open(tmp_filtered_consensus_output.name, "w") as f:
    #     with open(path_output_fasta, "r") as fasta:
    #         for record in SeqIO.parse(fasta, "fasta"):
    #             if record.id in consensus_sequences_present_in_cutreads_alignments and consensus_sequences_present_in_cutreads_alignments[record.id] >= min_cluster_size:
    #                 SeqIO.write(record, f, "fasta")
    #                 n_final_consensuses += 1
    # # replace path_output_fasta with tmp_filtered_consensus_output
    # shutil.move(tmp_filtered_consensus_output.name, path_output_fasta)
    # # index the filtered consensus output
    # cmd_index_fasta = split(f"samtools faidx {str(path_output_fasta)}")
    # subprocess.check_call(cmd_index_fasta)
    # ================================================================================================================

    # log_data["mcrID"] = mcrID
    # log_data["elapsed_time"] = elapsed_time
    # log_data["n_reads"] = len(reads)
    # log_data["n_consensuses"] = len(new_consensuses)
    # log_data["consensus_lengths"] = [consensus_meta[i][4] for i in consensus_meta.keys()]
    logstring = f"mcrID: {mcrID}; \
time: {elapsed_time}; \
reads: {len(reads)}; \
consensuses: {n_final_consensuses}; \
final_alignment_time: {time_align_reads}; \
time_consensus_creation: {time_consensus_creation}; \
time_read_partitioning: {time_read_partitioning}; \
n_partitions: {n_partitions}"
    log.info(logstring)


# %%


def run(args, **kwargs):
    consensus(
        mcrID=args.mcrID,
        path_crs=args.crs,
        path_reads_db=args.reads_db,
        path_intervals_db=args.intervals_db,
        path_output_fasta=args.output_fasta,
        path_output_bam=args.output_bam,
        min_signal_size=args.min_signal_size,
        reads_min_size=args.min_reads_size,
        min_bnd_size=args.min_bnd_size,
        reads_partitioning_threshold=args.reads_partitioning_threshold,
        threads=args.threads,
        tmp_dir_path=args.tmp_dir_path,
        logfile=args.logfile,
        use_handicap=args.no_handicap,
        use_preclustering=args.use_preclustering,
        min_cluster_size=args.min_cluster_size,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Create consensus sequences of cut reads in candidate region."
    )
    parser.add_argument(
        "-m", "--mcrID", type=int, required=True, help="candidate region ID"
    )
    parser.add_argument(
        "-c",
        "--crs",
        type=Path,
        required=True,
        help="path to (bgzipped) candidate regions tsv file",
    )
    parser.add_argument(
        "-r",
        "--reads-db",
        type=Path,
        required=True,
        help="path to sqlite3 database of reads",
    )
    parser.add_argument(
        "-i",
        "--intervals-db",
        type=Path,
        required=True,
        help="path to sqlite3 database of read intervals",
    )
    parser.add_argument(
        "-f",
        "--output-fasta",
        type=Path,
        required=True,
        help="path to output fasta file that will contain the consensus sequences",
    )
    parser.add_argument(
        "-b",
        "--output-bam",
        type=Path,
        required=True,
        help="path to output bam file that will contain the alignments of all reads of the candidate region to the consensus sequences",
    )
    parser.add_argument(
        "--min-signal-size",
        type=int,
        default=12,
        help="minimum signal size of deletions and insertions to be seen.",
    )
    parser.add_argument(
        "--min-reads-size",
        type=int,
        default=500,
        help="minimum size of reads to be considered for consensus creation.",
    )
    parser.add_argument(
        "--min-bnd-size",
        type=int,
        default=300,
        help="minimum size of a breakpoint to be considered.",
    )
    parser.add_argument(
        "--reads-partitioning-threshold",
        type=float,
        default=6.0,
        help="threshold of maximum summed indels of initial alignments of reads to the longest read of their cluster.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=3,
        help="number of threads to use for minimap2 and lamassemble.",
    )
    parser.add_argument(
        "--logfile", type=Path, default=None, required=False, help="path to logfile"
    )
    parser.add_argument(
        "--tmp-dir-path",
        default=None,
        required=False,
        help="path to temporary directory that does not delete any temporary files. Use for debugging only.",
    )
    parser.add_argument(
        "--no-handicap",
        action="store_false",
        help="use no handicap for the reads_partitioning_threshold.",
    )
    parser.add_argument(
        "--use-preclustering",
        action="store_true",
        help="use preclustering to partition reads.",
    )
    parser.add_argument(
        "--min-cluster-size",
        type=int,
        default=2,
        help="minimum number of cut reads in a cluster to be considered for a final consensus sequence.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    # if args.logfile:
    #     logzero.logfile(str(args.logfile))
    run(args)
    return


if __name__ == "__main__":
    main()

# # %%
# # dev scoring function for msas, given a representative
# # load msa. Its fasta-like, meaning: The line starting with '>' is the name of the sequence. The following line is the sequence.
# # the sequence has all four DNA letters and '-' for gaps.
# # the sequence is parsed to a numpy array of shape (n_seqs, n_bases)
# # DNA letters are encoded as 1,2,3,4 for A,C,G,T and 0 for gaps

# # first define a data class for a msa alignment object. It has a string:name and a np.array:seq of type uint8
# from dataclasses import dataclass
# @dataclass
# class read_alignment:
#     name: str
#     seq: np.ndarray

# def parse_msa(
#     path_msa:Path,
#     alphabet:dict={'-':0,'A':1,'C':2,'G':3,'T':4}) -> typing.List[read_alignment]:
#     with open(path_msa, "r") as f:
#         for line in f:
#             if line.startswith('>'):
#                 name = line[1:].strip()
#                 seq = np.array([alphabet[base] for base in next(f).strip()],dtype=np.uint8)
#                 yield read_alignment(name,seq)

# #%%

# id:int=6541
# path_msa = Path(f'/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/test/{id}.msa')

# msas = list(parse_msa(path_msa))
# # %%
# from scripts.signaldepth_to_crs import find_stretches
# # given a representative, score all reads in the msa against this representative
# # a weights array is computed from the array of sequences. It counts the number of non-gaps at each position.
# # and divides the number of non-gaps by the total number of sequences. This is the weight of each position.
# weights = np.sum([msa.seq != 0 for msa in msas],axis=0) / len(msas)
# # maybe don't use the gaps
# min_sv_size = 8
# # compute the number of mismatches at each position for each read compared to the representative.
# # pick the first read as representative for development purposes
# representative_name = msas[0].name
# # get a score for each read compared to the representative, but scale each position by the weight
# scores = []
# scores_no_weights = []
# for msa in msas:
#     if msa.name == representative_name:
#         continue
#     matches = msa.seq == msas[0].seq
#     # compute stretches of matches and mismatches. first and last are ignored
#     stretches = np.array(find_stretches(matches)[1:-1],dtype=np.uint16)
#     # subset stretches to rows with 0 in the third column and a difference of column 2 and column 1 > min_sv_size
#     big_gaps = stretches[(stretches[:,2] == 0) & (stretches[:,1] - stretches[:,0] > min_sv_size)]
#     # and sum them up
#     sum_big_gaps = np.sum(big_gaps[:,1] - big_gaps[:,0])
#     scores.append(sum_big_gaps)

# # score = np.sum((msa.seq != msas[0].seq) * weights)
# # # and compute the score without the weights
# # score_no_weights = np.sum(msa.seq != msas[0].seq)
# # scores.append(score)
# # scores_no_weights.append(score_no_weights)
# # %%
# import matplotlib.pyplot as plt
# plt.figure()
# plt.plot(sorted(scores))
# # %%
