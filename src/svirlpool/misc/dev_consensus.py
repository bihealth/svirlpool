# %%
import re
import subprocess
import tempfile
import typing
from math import ceil
from pathlib import Path
from shlex import split

import numpy as np
import pysam
import scipy.spatial.distance as ssd
from Bio import SeqIO
from logzero import logger as log
from numpy import typing as npt
from scipy.cluster.hierarchy import cut_tree, ward

from . import alignments_to_raf, datatypes, util

# %%


# compute raf
def rafs_from_alns(
    alns: typing.List[pysam.libcalignedsegment.AlignedSegment], min_signal_size: int = 6
) -> typing.List[datatypes.ReadAlignmentFragment]:
    rafs = []
    for aln in alns:
        if not aln.is_unmapped:
            has_signal, raf = alignments_to_raf.compute_raf(
                a=aln,
                sampleID=aln.query_name,
                min_signal_size=min_signal_size,
                min_mapq=0,
                min_clipped_length=50,
            )
            # all rafs are stored, so that perfect matches can be computed
            rafs.append(raf)
        else:
            log.info(f"Skipping unmapped read {aln.query_name}")
    return rafs


# create ava mapping
def create_ava_mapping(path_reads: Path, path_ava: Path) -> None:
    cmd_ava = f"minimap2 -a -x ava-ont -X --sam-hit-only {path_reads} {path_reads}"
    cmd_view = "samtools view -b"
    cmd_sort = "samtools sort"
    log.info(f"Running {cmd_ava}")
    cmd_index = f"samtools index {path_ava}"
    with open(path_ava, "wb") as f:
        p0 = subprocess.Popen(split(cmd_ava), stdout=subprocess.PIPE)
        p1 = subprocess.Popen(split(cmd_view), stdin=p0.stdout, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(split(cmd_sort), stdin=p1.stdout, stdout=f)
        p2.communicate()
        subprocess.check_call(split(cmd_index))


def compute_ava_rafs(path_ava, min_signal_size: int):
    ava_alns = list(pysam.AlignmentFile(path_ava, "rb"))
    ava_rafs = rafs_from_alns(ava_alns, min_signal_size)
    return ava_rafs


def create_read_dicts(path_reads):
    dict_index_readIDs = {
        i: read.id for i, read in enumerate(SeqIO.parse(path_reads, "fastq"))
    }
    dict_readIDs_index = {v: k for k, v in dict_index_readIDs.items()}
    return dict_index_readIDs, dict_readIDs_index


def score_raf_aligned_bases(a: datatypes.ReadAlignmentFragment) -> int:
    """computes the number of aligned bases"""
    return a.referenceAlignmentEnd - a.referenceAlignmentStart


def score_raf_svs(
    a: datatypes.ReadAlignmentFragment,
    min_signal_size: int = 30,
    min_bnd_size: int = 100,
    weight_bnd: int = 100,
) -> float:
    """computes the sum of the maximum number of bases of dels, ins and clipped tail"""
    sum_del = (
        sum(
            [
                deletion.ref_length if deletion.ref_length > min_signal_size else 0
                for deletion in a.deletions
            ]
        )
        if len(a.deletions) > 0
        else 0
    )
    sum_ins = (
        sum(
            [
                insertion.read_length if insertion.read_length > min_signal_size else 0
                for insertion in a.insertions
            ]
        )
        if len(a.insertions) > 0
        else 0
    )
    bnd = (
        a.breakendLeft.clippedLength > min_bnd_size if a.breakendLeft else False
    ) or (a.breakendRight.clippedLength > min_bnd_size if a.breakendRight else False)
    return float(sum_del + sum_ins + int(bnd) * weight_bnd)


def raf_distance_matrix(
    dict_readIDs_index, rafs, max_sv_density: float = 0.1, score_breakend: int = 200
):
    seen_pairs = set()
    M_sizes = np.zeros((len(dict_readIDs_index), len(dict_readIDs_index)), dtype=int)
    M_check = M_sizes.copy()
    M_svs = M_sizes.copy()
    M_svs.dtype = float
    for raf in sorted(rafs, key=lambda x: score_raf_aligned_bases(x), reverse=True):
        pair = tuple(sorted([dict_readIDs_index[raf.sampleID], int(raf.referenceID)]))
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        # iterate through rafs sorted by size on ref, biggest first.
        i, j = pair
        M_sizes[i, j] += score_raf_aligned_bases(raf)
        M_svs[i, j] += score_raf_svs(raf, min_bnd_size=score_breakend)
        M_check[i, j] += 1
    # reduce to [0,1] matrix
    M_sv_density = np.divide(
        M_svs, M_sizes, out=np.zeros_like(M_svs), where=M_sizes != 0
    )
    # exclude scores where SV signal is too much to be noise
    M = M_sizes * (M_sv_density < max_sv_density)
    M = M / M.max()
    # distance matrix
    M = 1 - M
    # make matrix symmetric
    M = (M + M.T) * 0.5
    # set diagonal to 0
    M = M * (1 - np.eye(M.shape[0]))
    return M


def initial_clustering(N_clusters: int, M: npt.ArrayLike) -> npt.ArrayLike:
    condensed_distance_matrix = ssd.squareform(M)
    linkage = ward(condensed_distance_matrix)
    return cut_tree(linkage, N_clusters).flatten()


def select_reads(
    fastq: Path, selected_readIDs: typing.List[str], fastq_out: Path
) -> None:
    tmp_readIDs = tempfile.NamedTemporaryFile(mode="w", delete=True, suffix=".txt")
    with open(tmp_readIDs.name, "w") as f:
        for readname in selected_readIDs:
            print(readname, file=f)
    cmd_subset = f"seqtk subseq {fastq} {tmp_readIDs.name}"
    with open(fastq_out, "wb") as f:
        subprocess.check_call(split(cmd_subset), stdout=f)


def first_clustering(path_reads: Path, n_clusters: int, min_signal_size: int):
    # tmp path to ava bam file
    tmp_ava = tempfile.NamedTemporaryFile(mode="w", delete=True, suffix=".bam")
    create_ava_mapping(path_reads, tmp_ava.name)
    ava_rafs = compute_ava_rafs(tmp_ava.name, min_signal_size)
    dict_index_readIDs, dict_readIDs_index = create_read_dicts(path_reads)
    ava_distance_matrix = raf_distance_matrix(dict_readIDs_index, ava_rafs)
    clustering = initial_clustering(n_clusters, ava_distance_matrix)
    return clustering, dict_index_readIDs


def make_consensus(
    name: str,
    path_reads: Path,
    path_consensus: Path,
    lamassemble_mat_path: Path,
    min_coverage: float = 0.1,
    threads: int = 1,
) -> typing.Tuple[Path, int]:
    """creates an indexed consensus.{mcrID}.fasta in the provided tmp_dir_path with lamassemble

    Args:
        crID (int): crID to be given as the name to the consensus sequence
        path_reads (Path): path to subset of reads that carry a signal and were selected from the original read sequences
        tmp_dir_path (Path): path to tmp or debugging dir
        lamassemble_mat_path (Path): path to lamassemble parameters matrix, e.g. lamassemble/train/promethion.mat
        min_coverage (float, optional): _description_. Defaults to 0.1.

    Returns:
        typing.Tuple[Path,int]: (path to consensus file, length of assembled sequence)
    """
    cmd_consensus = split(
        f"lamassemble \
        -n {name} \
        -s {str(min_coverage)}% \
        -P {str(threads)} \
        {str(lamassemble_mat_path)} \
        {str(path_reads)}"
    )
    consensus_file_handle = open(path_consensus, mode="w")
    subprocess.check_call(cmd_consensus, stdout=consensus_file_handle)
    # read consensus_path to list of lines
    h = open(str(path_consensus))
    next(h)
    consensus_length = len(next(h))
    if consensus_length < 30:
        return 0
    cmd_index = split(f"samtools faidx {str(path_consensus)}")
    subprocess.check_call(cmd_index)
    consensus_length = int(open(str(path_consensus) + ".fai").readline().split("\t")[1])
    return consensus_length


def make_consensuses_from_clusters(
    clustering: npt.ArrayLike,
    path_reads: Path,
    dict_index_readIDs: dict,
    path_base: Path,
    lamassemble_mat_path: Path,
    N: int,
) -> typing.Dict[int, Path]:
    """creates consensus sequences from clustered reads"""
    consensus_paths = {}
    for cluster in set(clustering):
        consensus_paths[cluster] = path_base / f"consensus.{cluster}.fasta"
        selected_readIDs = np.array(list(dict_index_readIDs.values()))[
            clustering == cluster
        ]
        tmp_fastq = path_base / f"reads.{cluster}.fastq"
        select_reads(path_reads, selected_readIDs, tmp_fastq.name)
        if (
            make_consensus(
                name=f"consensus.{N}.{cluster}",
                lamassemble_mat_path=lamassemble_mat_path,
                path_consensus=consensus_paths[cluster],
                path_reads=tmp_fastq.name,
            )
            == 0
        ):
            # write the first read's sequence to the consensus file
            with open(consensus_paths[cluster], "w") as f:
                seqrec = next(SeqIO.parse(tmp_fastq.name, "fastq"))
                seqrec.id = f"consensus.{N}.{cluster}"
                SeqIO.write(seqrec, f, "fasta")
    return consensus_paths


# %%
def align_reads_to_consensuses(
    consensus_paths: dict,
    reads: Path,
    path_base: Path,
    clustering=np.array([]),
    dict_index_readIDs={},
) -> typing.Dict[int, Path]:
    alignments_paths = {}
    for cluster in consensus_paths.keys():
        alignments_paths[cluster] = path_base / f"read_to_consensus.{cluster}.bam"
        # if clustering and dict is provided, then use only selected subset of reads
        if len(clustering) > 0 and len(dict_index_readIDs) > 0:
            print("using only selected reads")
            selected_readIDs = np.array(list(dict_index_readIDs.values()))[
                clustering == cluster
            ]
            fastq_used = path_base / f"reads.{cluster}.fastq"
            select_reads(reads, selected_readIDs, fastq_used.name)
        else:
            print("using all reads")
            fastq_used = reads
        util.align_reads_with_minimap(
            reference=consensus_paths[cluster],
            reads=fastq_used,
            bamout=alignments_paths[cluster],
            tech="map-ont",
            aln_args="-O5,32",
            threads=1,
        )
    return alignments_paths


def rafs_from_alignments(
    alignments_paths: dict,
) -> typing.Dict[int, typing.List[datatypes.ReadAlignmentFragment]]:
    """computes rafs from alignments prvided in a dict of clusterID to path to bam file"""
    return {
        cluster: rafs_from_alns(list(pysam.AlignmentFile(path, "rb")))
        for cluster, path in alignments_paths.items()
    }


def scores_from_rafs(
    dict_clusters_rafs: typing.Dict[int, Path],
    dict_index_readIDs: typing.Dict[int, str],
    clustering: npt.ArrayLike,
    min_signal_size: int = 30,
) -> typing.Dict[int, typing.Dict[int, float]]:
    readID_cluster_scores = {readID: {} for readID in dict_index_readIDs.values()}
    for cluster in set(clustering):
        for raf in dict_clusters_rafs[cluster]:
            readID_cluster_scores[raf.sampleID][cluster] = score_raf_svs(
                raf, min_signal_size
            )
    return readID_cluster_scores


def scores_clustering(readID_cluster_scores: dict, clustering: npt.ArrayLike) -> dict:
    sv_scores = np.array(
        [
            readID_cluster_scores[readID][cluster]
            for readID, cluster in zip(readID_cluster_scores.keys(), clustering)
        ],
        dtype=float,
    )
    rdict = {}
    for cluster in set(clustering):
        sv_scores_cluster = sv_scores[clustering == cluster]
        rdict[cluster] = sum(sv_scores_cluster) / len(sv_scores_cluster)
        # singletons also receive 0 score
        if len(sv_scores_cluster) == 1:
            rdict[cluster] = 0
    return rdict


def reassign_clusters(clustering: npt.ArrayLike, readID_cluster_scores: dict):
    new_clustering = clustering.copy()
    for i, ((readID, scores), j) in enumerate(
        zip(readID_cluster_scores.items(), clustering)
    ):
        # if scores[j] > 0.0:
        new_clustering[i] = min(scores, key=lambda x: scores[x])
    return new_clustering


def refine_clustering(clustering, readID_cluster_scores):
    # cluster_scores = scores_clustering(readID_cluster_scores,clustering)
    new_clustering = reassign_clusters(clustering, readID_cluster_scores)
    new_cluster_scores = scores_clustering(readID_cluster_scores, new_clustering)
    return new_clustering, new_cluster_scores


def is_n_insufficient(cluster_scores) -> bool:
    return sum(cluster_scores.values()) > 1


def get_readnames_of_clustering(
    n: int, refined_readID_cluster_scores: dict, clustering: npt.ArrayLike
) -> npt.ArrayLike:
    return np.array(list(refined_readID_cluster_scores.keys()))[clustering == n]


# %%


def find_structural_haplotypes(
    path_reads: Path,
    path_base: Path,
    lamassemble_mat: Path,
    N: int,
    min_signal_size: int = 6,
    min_sv_size: int = 30,
    n_clusters: int = 0,
):
    dict_index_readIDs, _ = create_read_dicts(path_reads)
    if n_clusters == 0:
        # infer number of initial clusters automatically
        n_clusters = 2 * len(
            set(
                list(
                    map(
                        lambda x: int(x.split(".")[-1]),
                        list(dict_index_readIDs.values()),
                    )
                )
            )
        )
        log.info(f"n_clusters not provided, setting to {n_clusters}")
    clustering, _ = first_clustering(path_reads, n_clusters, min_signal_size)
    consensus_paths = make_consensuses_from_clusters(
        clustering, path_reads, dict_index_readIDs, path_base, lamassemble_mat, N
    )
    alignments_paths = align_reads_to_consensuses(
        consensus_paths, path_reads, path_base
    )
    rafs_clustered = rafs_from_alignments(alignments_paths)
    readID_cluster_scores = scores_from_rafs(
        rafs_clustered, dict_index_readIDs, clustering, min_sv_size
    )
    refined_clustering, refined_cluster_scores = refine_clustering(
        clustering, readID_cluster_scores
    )
    refined_consensus_paths = make_consensuses_from_clusters(
        refined_clustering,
        path_reads,
        dict_index_readIDs,
        path_base,
        lamassemble_mat,
        N,
    )
    refined_alignments_paths = align_reads_to_consensuses(
        refined_consensus_paths, path_reads, path_base
    )
    # scoring
    refined_rafs_clustered = rafs_from_alignments(refined_alignments_paths)
    refined_readID_cluster_scores = scores_from_rafs(
        refined_rafs_clustered, dict_index_readIDs, refined_clustering, min_sv_size
    )
    # repeat with 2*n if necessary
    if is_n_insufficient(refined_cluster_scores):
        log.info(f"increasing n_clusters from {n_clusters} to {n_clusters*2}")
        return find_structural_haplotypes(
            path_reads,
            path_base,
            lamassemble_mat,
            N,
            min_signal_size,
            min_sv_size,
            n_clusters * 2,
        )
    return (
        refined_consensus_paths,
        refined_clustering,
        refined_alignments_paths,
        refined_readID_cluster_scores,
    )


def copy_consensus_alignments_to_destination(
    path_consensus: Path,
    refined_clustering: npt.ArrayLike,
    refined_readID_cluster_scores: dict,
    refined_alignments_paths: dict,
    N: int,
) -> None:
    """copies consensus alignments to destination path and creates index"""
    with tempfile.TemporaryDirectory() as tdir:
        path_tmp_dir = Path(tdir)
        # path_tmp_dir = Path("/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/test/213")
        tmp_bamfiles = []
        for cluster in set(refined_clustering):
            readnames = get_readnames_of_clustering(
                cluster, refined_readID_cluster_scores, refined_clustering
            )
            tmp_readnames = tempfile.NamedTemporaryFile(
                mode="w", delete=True, suffix=".txt"
            )
            path_readnames = Path(tmp_readnames.name)
            path_readnames = path_tmp_dir / f"readnames.{cluster}.txt"
            with open(path_readnames, "w") as f:
                for readname in readnames:
                    print(readname, file=f)
            # copy sam header to file
            path_tmp_sam = path_tmp_dir / f"tmp.{cluster}.sam"
            path_tmp_bam = path_tmp_dir / f"tmp.{cluster}.bam"
            cmd_copy_header = f"samtools view -H {refined_alignments_paths[cluster]}"
            cmd_filter_reads = f"samtools view {refined_alignments_paths[cluster]}"
            cmd_grep = f"grep -f {str(path_readnames)}"
            cmd_view = f"samtools view -b {str(path_tmp_sam)}"
            # copy header to file
            with open(path_tmp_sam, "w") as f:
                subprocess.check_call(split(cmd_copy_header), stdout=f)
            with open(path_tmp_sam, "a") as f:
                p0 = subprocess.Popen(split(cmd_filter_reads), stdout=subprocess.PIPE)
                p1 = subprocess.Popen(split(cmd_grep), stdin=p0.stdout, stdout=f)
                p1.communicate()
            # view to bam
            with open(path_tmp_bam, "wb") as f:
                subprocess.check_call(split(cmd_view), stdout=f)
            tmp_bamfiles.append(str(path_tmp_bam))
        # merge bam files and copy to consensus path
        consensus_path = path_consensus / f"consensus.{N}.bam"
        cmd_merge = (
            f"samtools merge -f -o {str(consensus_path)} {' '.join(tmp_bamfiles)}"
        )
        with open(consensus_path, "wb") as f:
            subprocess.check_call(split(cmd_merge), stdout=f)
        # index bam file
        cmd_index = f"samtools index {str(consensus_path)}"
        subprocess.check_call(split(cmd_index))


def copy_consensus_sequences_to_destination(
    path_consensus: Path, refined_consensus_paths: dict, N: int
):
    """copies consensus sequences to destination path and creates index"""
    output_path = path_consensus / f"consensus.{N}.fasta"
    cmd_cat = (
        f"cat {' '.join([str(path) for path in refined_consensus_paths.values()])}"
    )
    cmd_index = f"samtools faidx {str(output_path)}"
    subprocess.check_call(split(cmd_cat), stdout=output_path.open(mode="w"))
    subprocess.check_call(split(cmd_index))


def infer_N_from_fastq_path(path_reads: Path) -> int:
    """infers the id number of the CR name from the fastq file name"""
    N = int(re.findall(r"\d+", path_reads.name)[-1])
    if N:
        return N
    else:
        raise ValueError(f"Could not infer N from {path_reads.name}")


def consensus(
    path_reads: Path,
    lamassemble_mat: Path,
    path_consensus: Path,
    n_samples: int,
    N: int = -1,
    min_signal_size: int = 6,
    min_sv_size: int = 30,
) -> None:
    n_clusters = int(ceil(n_samples * 1.5))
    if N < 0:
        N = infer_N_from_fastq_path(path_reads)
    with tempfile.TemporaryDirectory() as tempdir:
        path_base = Path(tempdir)
        (
            refined_consensus_paths,
            refined_clustering,
            refined_alignments_paths,
            refined_readID_cluster_scores,
        ) = find_structural_haplotypes(
            path_reads,
            path_base,
            lamassemble_mat,
            N,
            min_signal_size,
            min_sv_size,
            n_clusters=n_clusters,
        )
        copy_consensus_alignments_to_destination(
            path_consensus,
            refined_clustering,
            refined_readID_cluster_scores,
            refined_alignments_paths,
            N,
        )
        copy_consensus_sequences_to_destination(
            path_consensus, refined_consensus_paths, N
        )


# %%

# path_reference = Path("/home/vinzenz/development/LRSV-detection/development/test/hs37d5.fa")

path_reads = Path(
    "/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/fastq/213.fastq"
)
lamassemble_mat = Path("/home/vinzenz/development/lamassemble/train/promethion.mat")
path_consensus = Path(
    "/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/consensus"
)

# %%

consensus(path_reads, lamassemble_mat, path_consensus, 2)
# %%
