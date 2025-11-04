# %%
# =============================================================================
# this function receives a cmd (list of args) and returns the result as a
# pandas dataframe

import csv
import gzip
import json
import pickle
import shlex
import shutil
import sqlite3
import subprocess
import tempfile
import typing
from collections import Counter
from enum import Enum
from io import StringIO
from math import floor
from pathlib import Path

import cattrs
import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval
from logzero import logger as log
from sklearn import preprocessing
from stopit import TimeoutException
from stopit import threading_timeoutable as timeoutable
from tqdm import tqdm

# %%
from . import alignments_to_rafs, consensus_class, datatypes

# %%


def yield_last_column(input: Path, dtype: type) -> typing.Iterable:
    csv.field_size_limit(2147483647)
    # open gzipped file if compressed, else open as text file
    try:
        with gzip.open(input, "rt") as f:
            reader = csv.reader(f, delimiter="\t", quotechar='"')
            for line in tqdm(reader):
                yield cattrs.structure(json.loads(line[-1]), dtype)
    except:
        try:
            with open(input, "r") as f:
                reader = csv.reader(f, delimiter="\t", quotechar='"')
                for line in tqdm(reader):
                    yield cattrs.structure(json.loads(line[-1]), dtype)
        except:
            raise FileNotFoundError(
                f"Could not open {input}. Make sure the file exists and is not corrupted."
            )


def yield_from_raf(input: Path) -> typing.Iterable[datatypes.ReadAlignmentFragment]:
    yield from yield_last_column(input=input, dtype=datatypes.ReadAlignmentFragment)


def yield_from_extendedSVsignal(
    input: Path,
) -> typing.Iterable[datatypes.ExtendedSVsignal]:
    yield from yield_last_column(input=input, dtype=datatypes.ExtendedSVsignal)


def yield_consensus_objects(
    path_db: Path, consensusIDs: set[str] | None = None, silent: bool = True
) -> typing.Generator[consensus_class.Consensus, None, None]:
    """produces a dict consensusID:consensus_sequence"""
    # iterate all consensus objects in database and construct a dict consensusID:consensusObject
    if not silent:
        log.info(f"loading consensus sequences from {path_db}...")
    try:
        conn = sqlite3.connect("file:" + str(path_db) + "?mode=ro", uri=True)
        c = conn.cursor()
    except sqlite3.OperationalError as e:
        log.error(
            f"Could not open database {path_db}. Make sure the file exists and is not corrupted."
        )
        raise e
    if not consensusIDs:
        c.execute("SELECT id,consensus FROM consensuses")
        for row in c:
            consensus_object: consensus_class.Consensus = cattrs.structure(
                pickle.loads(row[1]), consensus_class.Consensus
            )
            yield consensus_object
    else:
        placeholders = ",".join(
            ["?"] * len(consensusIDs)
        )  # Correctly format placeholders
        query = f"SELECT id, consensus FROM consensuses WHERE id IN ({placeholders})"
        c.execute(query, tuple(consensusIDs))  # Pass parameters correctly
        for row in c:
            consensus_object: consensus_class.Consensus = cattrs.structure(
                pickle.loads(row[1]), consensus_class.Consensus
            )
            yield consensus_object
    c.close()
    conn.close()


def load_consensus_sequences_from_db(
    path_db: Path, silent: bool = False
) -> dict[str, str]:
    """produces a dict consensusID:consensus_sequence"""
    # iterate all consensus objects in database and construct a dict consensusID:consensusObject
    if not silent:
        log.info(f"loading consensus sequences from {path_db}...")
    conn = sqlite3.connect("file:" + str(path_db) + "?mode=ro", uri=True)
    c = conn.cursor()
    c.execute("SELECT id,consensus FROM consensuses")
    results = {}
    for row in c:
        try:
            consensus_object: consensus_class.Consensus = cattrs.structure(
                pickle.loads(row[1]), consensus_class.Consensus
            )
            results[row[0]] = consensus_object.consensus_sequence
        except:
            import pdb

            pdb.set_trace()
    c.close()
    conn.close()
    return results


def load_alignments(path: Path) -> typing.List[pysam.AlignedSegment]:
    return list(pysam.AlignmentFile(path, "rb"))


def load_crs_connections(input: Path) -> dict:
    reader = csv.reader(input.open("r"), delimiter="\t")
    return {tuple(json.loads(l[0])): [set(k) for k in json.loads(l[1])] for l in reader}


def execute_to_df(cmd: list, sep="\t", header=None, index_col=None) -> pd.DataFrame:
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    csv = StringIO(process.stdout.read().decode())
    data = pd.read_csv(csv, header=header, index_col=index_col, sep=sep)
    csv.close()
    return data


def execute_to_np(cmd: list, sep="\t"):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    csv = StringIO(process.stdout.read().decode())
    data = np.loadtxt(csv, delimiter=sep)
    return data


def concat_txt_files(input_files: typing.List[Path], output: Path):
    with open(output, "w") as o:
        for f in input_files:
            with open(f, "r") as i:
                o.write(i.read())


def concat_binary_files(input_files: typing.List[Path], output: Path):
    with open(output, "wb") as o:
        for f in input_files:
            with open(f, "rb") as i:
                o.write(i.read())


def parse_alignment(
    aln: pysam.AlignedSegment,
    samplename: int | None = None,
    parse_sequence: bool = True,
    parse_qualities: bool = True,
) -> datatypes.Alignment:
    """Parse a pysam AlignedSegment into a datatypes.Alignment object."""
    # if sampleID is give, then just take the read name as it is. Else, remove last . separated part of the readname
    # if sampleID is None:
    #     readname = '.'.join(aln.query_name.split('.')[:-1]) if '.' in aln.query_name else aln.query_name
    #     sampleID = int(aln.query_name.split('.')[-1]) if '.' in aln.query_name else 0
    # else:
    #     readname = aln.query_name
    result = datatypes.Alignment(
        readname=aln.query_name,
        samplename=samplename,
        reference_name=aln.reference_name,
        reference_ID=aln.reference_id,
        reference_start=aln.reference_start,
        reference_end=aln.reference_end,
        samdict=aln.to_dict(),
        headerdict=aln.header.to_dict(),
    )
    if not parse_sequence:
        result.samdict["seq"] = "*"
    if not parse_qualities:
        # remove 'qual' key from samdict
        result.samdict["qual"] = "*"
    return result


@timeoutable()
def align_reads_with_minimap(
    reference: Path | str,
    reads: Path | str,
    bamout: Path | str,
    tech: str = "map-ont",
    threads: int = 3,
    aln_args: str = "",
    logfile: Path | None = None,
) -> None:
    try:
        cmd_align = shlex.split(
            f"minimap2 -a -x {tech} -t {threads} {aln_args} {reference} {reads}"
        )
        cmd_sort = shlex.split(f"samtools sort -O BAM -o {bamout}")
        cmd_index = shlex.split(f"samtools index {bamout}")
        if logfile:
            log_handle = open(logfile, "wt")
        log.info(" ".join(cmd_align))
        p_align = subprocess.Popen(
            cmd_align, stdout=subprocess.PIPE, stderr=log_handle if logfile else None
        )
        p_sort = subprocess.Popen(
            cmd_sort, stdin=p_align.stdout, stderr=log_handle if logfile else None
        )
        p_sort.communicate()
        log.info(f"indexing {bamout}")
        subprocess.check_call(cmd_index, stderr=log_handle if logfile else None)
        if logfile:
            log_handle.close()
    except TimeoutException:
        raise TimeoutException(
            "Alignment with minimap2 timed out after the specified time limit."
        )


# def align_reads_with_minimap(
#         reference:Path,
#         reads:str,
#         bamout:str,
#         tech:str="map-ont",
#         threads:int=3,
#         aln_args:str="",
#         logfile:Path|None=None) -> None:
#     cmd_align_ = f"minimap2 -a -x {tech} -t {threads} {aln_args} {reference} {reads}"
#     print(f"+ {cmd_align_}", file=sys.stderr)
#     cmd_align = shlex.split(cmd_align_)
#     cmd_compress = shlex.split("samtools view -b")
#     cmd_sort = shlex.split(f"samtools sort -o {bamout}")
#     cmd_index = shlex.split(f"samtools index {bamout}")
#     if logfile is None:
#         logfile = Path(f"{bamout}.log")
#     with open(logfile, "wt") as logf:
#         #logf = open(logfile, "wt")
#         p_align = subprocess.Popen(cmd_align, stdout=subprocess.PIPE, stderr=logf)
#         p_compress= subprocess.Popen(cmd_compress,stdin=p_align.stdout, stdout=subprocess.PIPE, stderr=logf)
#         p_align.stdout.close()
#         p_sort= subprocess.Popen(cmd_sort,stdin=p_compress.stdout, stderr=logf)
#         p_compress.stdout.close()
#         p_sort.communicate()
#         subprocess.call(cmd_index)
#     #logf.close()


def align_reads_with_ngmlr(
    reference: Path,
    reads: str,
    bamout: Path,
    tech: str = "ont",
    threads: int = 3,
    aln_args: str = "",
):
    # first cat all consensus sequences into one fasta file
    cmd_cat = shlex.split("cat " + reads)
    # write to temp fasta
    tmp_fasta = tempfile.NamedTemporaryFile(mode="w", delete=True, suffix=".fasta")
    with open(tmp_fasta.name, "w") as f:
        subprocess.check_call(cmd_cat, stdout=f)
        reads = tmp_fasta.name
    cmd_align = shlex.split(
        f"ngmlr -t {threads} -x {tech} {aln_args} -r {reference} -q {reads}"
    )
    log.info(" ".join(cmd_align))
    cmd_compress = shlex.split("samtools view -b")
    cmd_sort = shlex.split(f"samtools sort -o {bamout}")
    cmd_index = shlex.split(f"samtools index {bamout}")
    p_align = subprocess.Popen(cmd_align, stdout=subprocess.PIPE)
    p_compress = subprocess.Popen(
        cmd_compress, stdin=p_align.stdout, stdout=subprocess.PIPE
    )
    p_align.stdout.close()
    p_sort = subprocess.Popen(cmd_sort, stdin=p_compress.stdout)
    p_compress.stdout.close()
    p_sort.communicate()
    subprocess.call(cmd_index)


def align_reads_with_winnowmap(
    reference: Path, reads: Path, bamout: Path, tech: str, rep_k15: Path = Path("")
):
    if tech not in ["map-pb", "map-ont", "asm5", "asm10", "asm20"]:
        raise ValueError(f"tech must be either 'map-pb' or 'map-ont', not {tech}")
    # if no rep_k15 is provided, then it must be computed beforehand.
    # create temporary merylDB file
    tmp_merylDB = tempfile.NamedTemporaryFile(suffix=".merylDB", delete=False)
    rep_k15_file = tempfile.NamedTemporaryFile(suffix=".rep15.txt", delete=False)
    if rep_k15 == Path(""):
        rep_k15 = Path(rep_k15_file.name)
        cmd_k15_count = shlex.split(
            f"meryl count k=15 output {tmp_merylDB.name} {reference}"
        )
        cmd_k15_print = shlex.split(
            f"meryl print greater-than distinct=0.9998 {tmp_merylDB}"
        )
        # subprocess count
        subprocess.check_call(cmd_k15_count)
        # subprocess print, write output to rep_k15
        with open(rep_k15, "w") as f:
            subprocess.check_call(cmd_k15_print, stdout=f)
    tmp_merylDB.close()
    # reads need to be gzipped?
    # tech: map-pb, map-ont.
    cmd_align = shlex.split(
        f"winnowmap -Y -W {rep_k15} --MD --eqx -ax {tech} {reference} {reads}"
    )
    cmd_compress = shlex.split("samtools view -b")
    cmd_sort = shlex.split(f"samtools sort -o {bamout}")
    cmd_index = shlex.split(f"samtools index {bamout}")
    p_align = subprocess.Popen(cmd_align, stdout=subprocess.PIPE)
    p_compress = subprocess.Popen(
        cmd_compress, stdin=p_align.stdout, stdout=subprocess.PIPE
    )
    p_align.stdout.close()
    p_sort = subprocess.Popen(cmd_sort, stdin=p_compress.stdout)
    p_compress.stdout.close()
    p_sort.communicate()
    subprocess.call(cmd_index)
    # remove temporary files
    rep_k15_file.close()


# --- data loading --- #


def load_crs(path_crs: Path) -> typing.List[datatypes.CandidateRegion]:
    csv.field_size_limit(2147483647)
    candidate_regions = []
    with open(path_crs, "r") as f:
        reader = csv.reader(f, delimiter="\t", quotechar='"')
        for line in reader:
            candidate_region = cattrs.structure(
                json.loads(line[-1]), datatypes.CandidateRegion
            )
            candidate_regions.append(candidate_region)
    return candidate_regions


def create_ref_dict(reference: Path) -> typing.Dict[int, str]:
    """returns a dictionary that maps reference IDs to chromosome names"""
    return {
        i: str(line.rstrip().split("\t")[0])
        for i, line in enumerate(open(Path(str(reference) + ".fai"), "r"))
    }


def bed_chr_to_chrID(
    input: Path, reference: Path, output: Path = Path(""), add_id: bool = True
):
    """loads a bed file (input) and converts the chromosome names to chromosome IDs.
    Adds an ID to the end of each line."""
    ref_dict = create_ref_dict(reference)
    ref_dict = {v: k for k, v in ref_dict.items()}
    if output != Path(""):
        f = open(output, "w")
    for i, line in enumerate(open(input, "r")):
        line = line.rstrip().split("\t")
        line[0] = ref_dict[str(line[0])]
        if add_id:
            line.append(i)
        s = "\t".join(list(map(str, line)))
        if output != Path(""):
            print(s, file=f)
        else:
            print(s)


def bed_chrID_to_chr(input: Path, reference: Path, output: Path = Path("")):
    """loads a bed file (input) and converts the chromosome IDs to chromosome names."""
    ref_dict = create_ref_dict(reference)
    if output != Path(""):
        f = open(output, "w")
    for line in open(input, "r"):
        line = line.rstrip().split("\t")
        line[0] = ref_dict[int(line[0])]
        s = "\t".join(list(map(str, line)))
        if output != Path(""):
            print(s, file=f)
        else:
            print(s)


def create_fai_if_not_exists(reference: Path) -> Path:
    if not Path(str(reference) + ".fai").exists():
        log.warning(f"{str(reference)+'.fai'} not found. Trying to create one..")
        try:
            cmd = shlex.split(f"samtools faidx {reference}")
            subprocess.check_call(cmd)
        except:
            raise FileNotFoundError(
                f"{str(reference)+'.fai'} not found. And could not be created. Make sure to provide a path to a .fasta file for which a .fai file exists in the same directory."
            )
    return Path(str(reference) + ".fai")


def fai_to_bed(reference: Path, output: Path) -> None:
    create_fai_if_not_exists(reference=reference)
    with open(output, "w") as f:
        for line in open(Path(str(reference) + ".fai"), "r"):
            line = line.rstrip().split("\t")
            print(line[0], 0, line[1], sep="\t", file=f)


def genome_file_for_bedtools(
    reference: Path, output: Path | None = None
) -> dict[str, int]:
    create_fai_if_not_exists(reference=reference)
    rdict = dict()
    if output:
        with open(output, "w") as f:
            for line in open(Path(str(reference) + ".fai"), "r"):
                line = line.rstrip().split("\t")
                print(f"{line[0]}\t{line[1]}", file=f)
    else:
        for line in open(Path(str(reference) + ".fai"), "r"):
            line = line.rstrip().split("\t")
            rdict[line[0]] = int(line[1])
    return rdict


# =============================================================================
#  complexity
# =============================================================================

import numpy.typing as npt

# def hash_dna(seq,alphabet_hash_dict:dict={dna:np.uint8(i) for i,dna in enumerate('ACGTN')}):
#     for s in seq:
#         yield alphabet_hash_dict[s.upper()]

# def hash_dna4(seq:npt.NDArray[np.float16]) -> np.uint16:
#     """Hash DNA sequence with 4-letter alphabet (ACGT) using 2-bit encoding.
#     Values should be in range 0-3."""
#     value:int = 0
#     for s in seq:
#         value = (value << 2)+int(s)
#     return np.uint16(value)


def hash_dna5(seq: npt.NDArray[np.float16]) -> np.uint16:
    """Hash DNA sequence with 5-letter alphabet (ACGT + N/unknown) using 3-bit encoding.
    Values should be in range 0-4."""
    value: int = 0
    for s in seq:
        value = (value << 3) + int(s)  # 3 bits per base (supports 0-7)
    return np.uint16(value)


# def kmers_on_dna4(seq:npt.NDArray[np.float16],k:int) -> npt.NDArray[np.float16]:
#     n = len(seq)-k+1
#     return np.array([hash_dna4(seq[i:i+k]) for i in range(0,n)],dtype=np.uint16)


def kmers_on_dna5(seq: npt.NDArray[np.float16], k: int) -> npt.NDArray[np.float16]:
    """Compute k-mer hashes for 5-letter alphabet (ACGT + N/unknown)."""
    n = len(seq) - k + 1
    return np.array([hash_dna5(seq[i : i + k]) for i in range(0, n)], dtype=np.uint16)


def complexity_local_track(
    dna_iter: typing.Iterator,
    w: int = 11,
    K: list[int] = [1, 2, 3, 4, 5],
    padding: bool = False,
) -> npt.NDArray[np.float16]:
    """
    Calculate local sequence complexity along a DNA sequence.

    Args:
        dna_iter: Iterator over DNA sequence characters
        w: Window size for complexity calculation
        K: List of k-mer sizes to use for complexity calculation
        padding: If True, pad output to match input sequence length

    Returns:
        Array of complexity values (0-1 range)
    """
    # Set up alphabet based on use_5letter flag
    # 5-letter alphabet: A=0, C=1, G=2, T=3, N=4, unknown=4
    alphabet_hash_dict = {
        "A": np.uint8(0),
        "a": np.uint8(0),
        "C": np.uint8(1),
        "c": np.uint8(1),
        "G": np.uint8(2),
        "g": np.uint8(2),
        "T": np.uint8(3),
        "t": np.uint8(3),
        "N": np.uint8(4),
        "n": np.uint8(4),
    }

    default_value = 4  # Default for unknown letters

    kmers = [np.zeros(w - k + 1, dtype=np.uint8) for k in K]
    track: list[float] = []
    # init
    MAXCOUNTS = [min(5**k, w - k + 1) for k in K]
    # init iteration
    try:
        dna_window = np.array(
            list(
                map(
                    lambda x: alphabet_hash_dict.get(x, default_value),
                    [next(dna_iter) for _ in range(w)],
                )
            ),
            dtype=np.uint8,
        )
    except StopIteration:
        log.warning(
            f"complexity_local_track: sequence is shorter than window size {w}. Returning empty array."
        )
        return np.array([], dtype=np.float16)

    # calc kmer codes
    for i, k in enumerate(K):
        kmers[i] = kmers_on_dna5(seq=dna_window, k=k)

    track.append(np.prod([len(set(kmers[i])) / MAXCOUNTS[i] for i in range(len(K))]))
    # padding
    if padding:
        # Prepend padding values instead of replacing
        pad_count = int(floor(w / 2))
        track = [track[0]] * pad_count + track
    # loop iteration
    j = w
    while True:
        try:
            letter = next(dna_iter)
            dna_window[:-1] = dna_window[1:]
            dna_window[-1] = alphabet_hash_dict.get(letter, default_value)
            # update kmers
            for i, k in enumerate(K):
                if j % k == 0:
                    kmers[i][:-1] = kmers[i][1:]
                    kmers[i][-1] = hash_dna5(dna_window[-k:])
            track.append(
                np.prod([len(set(kmers[i])) / MAXCOUNTS[i] for i in range(len(K))])
            )
            j += 1
        except StopIteration:
            break

    # stop iteration
    # padding
    if padding:
        # Append padding values at the end
        pad_count = int(floor(w / 2))
        for i in range(pad_count):
            track.append(track[-1])
    return np.array(track, dtype=np.float16)


# def kmer_hashes_from_sequence(K:int,seq:str) -> npt.NDArray[np.uint16]:
#     return kmers_on_dna4(np.array(hash_dna(seq),dtype=np.uint32),K)

# def kmer_counter_dna4(seq:str,K:int) -> dict[int,int]:
#     # make sure all letters are in the alphabet "ACGT"
#     assert all([s in 'ACGT' for s in seq])
#     kmers = kmer_hashes_from_sequence(K,seq)
#     # cast to int
#     kmers = [int(k) for k in kmers]
#     return dict(Counter(kmers))

# try a different approach of computing complexity
# the idea is to iterate a sequence of length N
# with a range of windows of size s0, s1, ..., sm
# e.g. S = [1,2,3]
# then iterate the sequence and check if the current window equals the last non-overlapping window
# def sequence_complexity(seq:str,window_sizes:list[int]=[1,2,3,4,5]) -> np.ndarray:
#     result = np.zeros(len(seq),dtype=np.float16)
#     for w in window_sizes:
#         for i in range(len(seq)-w):
#             if i-w < 0:
#                 result[i] += 1.0
#             else:
#                 if seq[i-w:i] == seq[i:i+w]:
#                     result[i] += 0.0
#                 else:
#                     result[i] += 1.0
#         result[-w:] += 1.0
#     result = result / len(window_sizes)
#     return result


# =============================================================================
#  read intervals
# =============================================================================


def get_interval_on_read_in_region(
    a: pysam.libcalignedsegment.AlignedSegment,
    start: int,
    end: int,
    buffer_clipped_length: int = 0,
) -> typing.Tuple[int, int]:
    # find start and end position on ref on read
    start, end = (end, start) if start > end else (start, end)
    istart = get_ref_position_on_read(
        alignment=a,
        position=start,
        direction=Direction.LEFT,
        buffer_clipped_length=buffer_clipped_length,
    )
    iend = get_ref_position_on_read(
        alignment=a,
        position=end,
        direction=Direction.RIGHT,
        buffer_clipped_length=buffer_clipped_length,
    )
    return (iend, istart) if a.is_reverse else (istart, iend)


def get_interval_on_ref_in_region(
    a: pysam.libcalignedsegment.AlignedSegment, start: int, end: int
) -> typing.Tuple[int, int]:
    # find start and end position on ref
    start, end = (end, start) if start > end else (start, end)
    istart = get_read_position_on_ref(
        alignment=a, position=start, direction=Direction.LEFT
    )
    iend = get_read_position_on_ref(
        alignment=a, position=end, direction=Direction.RIGHT
    )
    return istart, iend


# %%
# =============================================================================
#  used only in tests
# =============================================================================


def generate_sequence(size: int, seed: int = 0) -> typing.List[str]:
    np.random.seed(seed)
    return [np.random.choice(list("ACGT"), p=[0.3, 0.2, 0.2, 0.3]) for _ in range(size)]


def delete_interval(seq: list, a: int, b: int):
    return [*seq[:a], *seq[b:]]


def insertion(
    seq: list, pos: int, size: int = 0, sequence: list = [], seed: int = 0
) -> typing.List[str]:
    if sequence == [] and size == 0:
        return seq
    if sequence == []:
        sequence = generate_sequence(size=size, seed=seed)
    r = [*seq[:pos], *sequence, *seq[pos:]]
    return [*seq[:pos], *sequence, *seq[pos:]]


def duplication(
    seq: list, a: int = 0, b: int = 0, copynumber: int = 2, insertion_site: int = None
) -> typing.List[str]:
    if b == 0:
        b = len(seq)
    if insertion_site:
        return [
            *seq[:insertion_site],
            *(seq[a:b] * (copynumber - 1)),
            *seq[insertion_site:],
        ]
    return [*seq[:a], *(seq[a:b] * copynumber), *seq[b:]]


def complement(seq: list) -> typing.List[str]:
    d = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return list(map(lambda s: d[s], seq))


def reverse_complement(seq: list, _=None, __=None) -> typing.List[str]:
    return complement(seq[::-1])


def inversion(seq: list, a: int = 0, b: int = 0) -> typing.List[str]:
    if b == 0:
        b = len(seq)
    return [*seq[:a], *complement(seq[a:b])[::-1], *seq[b:]]


def translocation(
    seq: list, gapStart: int, gapLength: int, fillPos: int
) -> typing.List[str]:
    seqcache = seq[gapStart : (gapStart + gapLength)]
    seq = [*seq[:gapStart], *seq[(gapStart + gapLength) :]]
    return [*seq[:fillPos], *seqcache, *seq[fillPos:]]


def add_noise(seq: list, noise_level: float = 0.1, seed: int = 0) -> typing.List[str]:
    np.random.seed(seed)
    return [
        (
            s
            if np.random.rand() > noise_level
            else np.random.choice(list(set("ACGT") - set(s)))
        )
        for s in seq
    ]


# make indexed alignments with indexed reference in provided dir
# e.g.
# seq=generate_sequence(100)
# rules = [(duplication,20,50,2),
#          (inversion,80,110),
#          (delete_interval,110,120)]
# apply_seq_synth_rules(seq,rules)
def apply_seq_synth_rules(seq: list, rules: typing.List[tuple]) -> typing.List[str]:
    r = seq.copy()
    for t in rules[::-1]:
        r = t[0](r, *t[1:])
    return r


# a rule looks like this:
# (function, *args)
# where function is one of the functions defined above and args are the arguments,
# e.g. (duplication,20,50,2)
def generate_read_from_rules(
    refStart: int, refEnd: int, rules: typing.List[typing.Tuple], reference: list
) -> typing.List[str]:
    # calc final extends of read from rules
    # rules are applied to reference
    ref = apply_seq_synth_rules(reference, rules)
    # read is generated from applied reference
    return ref[refStart:refEnd]


def write_sequences_to_fasta(
    seqs: typing.List[list],
    path: Path,
    chrnames: bool = False,
    prefix: str = "seq",
    suffix: str = "",
    dir: Path | None = None,
):
    if chrnames:
        seqnames = [prefix + str(i) + suffix for i in range(len(seqs))]
    else:
        seqnames = [
            prefix + str(hash("".join(["".join(s) for s in seqs])))[:5] + suffix
        ] * len(seqs)
    records = [
        SeqRecord(
            Seq("".join(seqs[i])),
            id=seqnames[i],
            name=seqnames[i],
            description=f"synthetic sequence of length {str(len(seqs[i]))}",
        )
        for i in range(len(seqs))
    ]
    tmp_buffer = tempfile.NamedTemporaryFile(
        dir=dir, mode="w", suffix=".fasta", delete=True
    )
    SeqIO.write(records, tmp_buffer.name, "fasta")
    # then format with seqtk seq -l 60

    cmd_format = shlex.split(f"seqtk seq -l 60 {tmp_buffer.name}")
    with open(path, "w") as f:
        subprocess.check_call(cmd_format, stdout=f)


def index_reference(path: Path):
    cmd = shlex.split(f"samtools faidx {path}")
    subprocess.check_call(cmd)


# re-implement, but so that the path to the alignments is passed as a parameter
def create_alignments_to_reference(
    reads: list[list[str]],
    reference: list[str],
    alignments: Path | None = None,
    aln_args: str | None = None,
) -> list[pysam.AlignedSegment]:
    tmp_dir = tempfile.TemporaryDirectory(suffix="tmp_alignment_creation.")
    # write reference to tmp fasta
    tmp_reference = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=True, dir=tmp_dir.name
    )
    ref_record = SeqRecord(
        Seq("".join(reference)),
        id="ref",
        name="ref",
        description=f"synthetic reference of length {str(len(reference))}",
    )
    with open(tmp_reference.name, "w") as f:
        SeqIO.write(ref_record, f, "fasta")
    # index reference
    index_reference(tmp_reference.name)
    # write reads to tmp fasta
    tmp_reads = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=True, dir=tmp_dir.name
    )
    reads_records = [
        SeqRecord(
            Seq("".join(r)),
            id=f"read-{i}",
            name=f"read-{i}",
            description=f"synthetic read of length {str(len(r))}",
        )
        for i, r in enumerate(reads)
    ]
    with open(tmp_reads.name, "w") as f:
        SeqIO.write(reads_records, f, "fasta")
    # align reads to reference
    if not alignments:
        tmp_alignments = tempfile.NamedTemporaryFile(
            mode="w", suffix=".bam", delete=True, dir=tmp_dir.name
        )
        alignments = Path(tmp_alignments.name)
    align_reads_with_minimap(
        aln_args=aln_args if aln_args else "",
        bamout=alignments,
        reference=Path(tmp_reference.name),
        reads=tmp_reads.name,
        tech="map-ont",
        threads=1,
    )
    # parse alignments to list of AlignedSegments
    alns = list(pysam.AlignmentFile(alignments, "rb"))
    # return list of AlignedSegments
    # close tmp_dir
    tmp_dir.cleanup()
    return alns


def get_alignments_to_reference(
    reads: list[SeqRecord], reference: SeqRecord, aln_args: str = ""
) -> list[pysam.AlignedSegment]:
    alns: list[pysam.AlignedSegment] = []
    with tempfile.TemporaryDirectory(suffix="tmp_get_alignments") as tdir:
        path_reads = Path(tdir) / "reads.fasta"
        path_reference = Path(tdir) / "reference.fasta"
        path_alignments = Path(tdir) / "alignments.bam"
        # write reference to tmp fasta
        with open(path_reference, "w") as f:
            SeqIO.write(reference, f, "fasta")
        # index reference
        index_reference(path_reference)
        # write reads to tmp fasta
        with open(path_reads, "w") as f:
            SeqIO.write(reads, f, "fasta")
        # align reads to reference
        align_reads_with_minimap(
            aln_args=aln_args,
            bamout=path_alignments,
            reference=path_reference,
            reads=path_reads,
            tech="map-ont",
            threads=1,
        )
        # parse alignments to list of AlignedSegments
        alns = list(pysam.AlignmentFile(path_alignments, "rb"))
    # return list of AlignedSegments
    return alns


def get_samplenames_from_samfile(samfile: pysam.Samfile) -> typing.List[str]:
    """returns a list of unique sample names from the RG tags in the samfile header."""
    try:
        return list(set([item["SM"] for item in samfile.header.to_dict()["RG"]]))
    except KeyError:
        return []


def create_sample_dict_from_alignments(
    alignments: typing.List[Path],
) -> typing.Tuple[typing.Dict[str, Path], typing.Dict[str, Path]]:
    """creates two dictionaries. The first one maps sample names to the path of the alignments file.
    The sample name is either taken from the RG tags in the alignment file's headers or
    if no RG tag is present, then the sample name is the file name stem."""
    sample_dict = {}
    for a in alignments:
        samfile = pysam.AlignmentFile(a.absolute(), "rb")
        if "RG" not in samfile.header.to_dict().keys():
            sample_dict[a.stem] = str(a.absolute())
        else:
            sns = get_samplenames_from_samfile(samfile)
            if len(sns) == 1:
                sample_dict[sns[0]] = str(a.absolute())
            if len(sns) > 1:
                raise ValueError(f"more than one sample name in {a.absolute()}")
    return sample_dict


def generate_sampledicts(alignments: typing.List[Path], output: Path):
    """creates three dictionaries. The first one maps sample names to the path of the alignments file.
    The second one maps sample names to sample IDs. The third maps sample IDs to the path of the alignments file.
    """
    dict_samplename_to_path = create_sample_dict_from_alignments(alignments)
    dict_samplename_to_id = {
        key: i for i, key in enumerate(dict_samplename_to_path.keys())
    }
    dict_id_to_bam = {
        v: dict_samplename_to_path[k] for k, v in dict_samplename_to_id.items()
    }
    with open(output, "w") as f:
        json.dump([dict_samplename_to_path, dict_samplename_to_id, dict_id_to_bam], f)


def load_sampledicts(input: Path) -> typing.List[dict]:
    """loads three dictionaries. The first one maps sample names to the path of the alignments file.
    The second one maps sample names to sample IDs. The third maps sample IDs to the path of the alignments file.
    Returns:
        typing.Tuple[typing.Dict[str,Path],typing.Dict[str,int]]
    """
    with open(input, "r") as f:
        dicts = json.load(f)
    # correct sampleIDs in dicts[2]: change type of keys to ints
    dicts[2] = {int(k): v for k, v in dicts[2].items()}
    return dicts


# %%


# plot to terminal a quick visualization that works like imshow with a provided list
# of ascii letters as shades from light to dark.
def plot_to_terminal(data: list, shades: list = [".", "-", "+", "#"]):
    ndata = preprocessing.normalize(data)
    for row in ndata:
        for col in row:
            print(shades[floor(col * (len(shades) - 1))], end="")
        print()


def overlapping_interval(
    A: typing.Tuple[int, int, int], B: typing.Tuple[int, int, int]
) -> int:
    """A and B must be triples of format
    (referenceID,start,end) and of type int. returns size of overlap."""
    if A[0] != B[0]:
        return 0
    A, B = (B, A) if B[1] < A[1] else (A, B)
    return max(0, min(A[2], B[2]) - max(A[1], B[1]))


def overlap(a: typing.Tuple[int, int], b: typing.Tuple[int, int]) -> int:
    """returns the size of overlap. returns the the size of the overlap [1 to n] or the distance between the intervals [-1 to -n]"""
    # test if all values are positive and in order
    if a[0] < 0 or a[1] < 0 or b[0] < 0 or b[1] < 0:
        raise ValueError(f"all values must be positive. given values: {a} and {b}")
    if a[0] > a[1]:
        a = (a[1], a[0])
    if b[0] > b[1]:
        b = (b[1], b[0])
    a, b = (a, b) if a[0] <= b[0] else (b, a)
    if a[1] - a[0] == 0 and b[0] <= a[0] <= b[1]:
        return 1
    if b[1] - b[0] == 0 and a[0] <= b[0] <= a[1]:
        return 1
    if a[1] == b[0]:
        return 1
    return min(a[1], b[1]) - max(a[0], b[0])


# %%
def get_reference_sequence_from_read_alignment(
    a: pysam.AlignedSegment, reference: Path
) -> Seq:
    """returns the reference sequence that is covered by the read alignment."""
    # get start and end of alignment
    start = a.reference_start
    end = a.reference_end
    # get reference sequence
    cmd = shlex.split(f"samtools faidx {reference} {a.reference_name}:{start}-{end}")
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # parse output to string
    out_str = process.stdout.read().decode()
    _, seq = out_str.split("\n", 1)
    process.stdout.close()
    return Seq(seq)


def get_unaligned_intervals(
    total_length: int, covered_intervals: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    # sort all intervals by start and end
    # assert: no overlaps allowed
    # iterate all intervals and save the gaps bewteen a[1],b[0] with a <= b
    # in the range 0,total_lenght
    # 0 to a and end to total_length are handled separately
    sorted_intervals = sorted(covered_intervals, key=lambda x: (x[0], x[1]))
    assert all(
        [
            sorted_intervals[i][1] <= sorted_intervals[i + 1][0]
            for i in range(len(sorted_intervals) - 1)
        ]
    )

    unaligned_intervals = []
    first_interval = (0, sorted_intervals[0][0])
    if first_interval[1] > 0:
        unaligned_intervals.append(first_interval)
    for i in range(len(sorted_intervals) - 1):
        if sorted_intervals[i][1] < sorted_intervals[i + 1][0]:
            unaligned_intervals.append(
                (sorted_intervals[i][1], sorted_intervals[i + 1][0])
            )
    last_interrval = (sorted_intervals[-1][1], total_length)
    if last_interrval[0] < total_length:
        unaligned_intervals.append(last_interrval)
    return unaligned_intervals


def add_hard_clips_to_aligned_tails(input: Path, output: Path) -> None:
    # iterate alignments in input
    # each alignment has a query_name of the form: name _ readstart _ readend
    # all aligned tails need to get the is_supplementary flag set to True
    # An alignment has the componenets:
    # hard / soft clipped start, matched segment, hard/ soft clipped end.
    # All soft clips need to be converted to hard clips (trim the query_sequence and qualities)
    # then all hard clips are merged with the hard clips taken from the query_name
    # the query_name is reduced to the first element after splitting by '_'
    # SA tags of all alignments need updates:
    # add their ref,pos,strand,cigarstring,mapq,?;
    pass


def add_clipped_tails_to_alignments(
    alignments: Path | str,
    reference: Path | str,
    output: Path | str,
    tech: str = "asm5",
    aln_args: str = "",
    threads: int = 0,
    min_clipped_length: int = 300,
    verbose: bool = False,
) -> None:
    log.warning(
        "%s is deprecated and will be removed in the future. It has no use anymore since aligned virtual fragments are now padded with read sequence."
    )
    """adds clipped tails to alignments. The clipped tails are aligned to the reference and added to the original alignment. Known issue: hard clips are missing on the added tails."""
    read_covered_intervals: dict[str, list[tuple]] = dict()
    readlengths: dict[str, int] = dict()
    # selected_readalignments = dict()
    for aln in pysam.AlignmentFile(
        alignments, "rb"
    ):  # secondary and supplementary alignments are not considered
        if not aln.is_mapped or aln.is_secondary:
            continue
        rl = aln.infer_read_length()  # total read length with clipped ends
        readlengths[aln.query_name] = rl
        s, e = alignments_to_rafs.get_start_end(aln)[
            2:
        ]  # get start and end of read interval in unclipped read coords
        if aln.query_name not in read_covered_intervals:
            read_covered_intervals[aln.query_name] = []
        read_covered_intervals[aln.query_name].append((s, e))
        # selected_readalignments[a.query_name] = a
    if verbose:
        log.info(f"found {len(read_covered_intervals)} reads with clipped tails")
    if len(read_covered_intervals) == 0:
        log.info("no reads with clipped tails found. Copying input to output.")
        shutil.copy(alignments, output)
        return

    # open temporary dir to store fasta files and save alignments to.
    tmp_dir = tempfile.TemporaryDirectory()
    # create temporary fasta file in tmp_dir
    tmp_fasta_unaligned_tails = tempfile.NamedTemporaryFile(
        suffix=".fasta", dir=tmp_dir.name, delete=False
    )
    # sort each read's intervals and process
    unaligned_intervals: dict[str, list[tuple]] = dict()
    with open(tmp_fasta_unaligned_tails.name, "w") as f:
        for aln in pysam.AlignmentFile(alignments, "rb"):
            # check if the alignment is primary (has no hard clippings)
            if aln.is_secondary or aln.is_supplementary:  # filter all non-primary
                continue

            unaligned_intervals[aln.query_name] = get_unaligned_intervals(
                total_length=readlengths[aln.query_name],
                covered_intervals=read_covered_intervals[aln.query_name],
            )
            # save each unaligned interval whose size exceeds min_clipped_length
            # to the temporary fasta file. add to the name its interval
            # the uncovered sequences hav to be revcomp if the alignment is revcomp to ensure
            # direction consistency
            seq_full = (
                Seq(aln.query_sequence).reverse_complement()
                if aln.is_reverse
                else Seq(aln.query_sequence)
            )
            for s, e in unaligned_intervals[aln.query_name]:
                if e - s >= min_clipped_length:
                    seq = seq_full[s:e]
                    print(f">{aln.query_name}_{s}_{e}", file=f)
                    print(seq, file=f)

    log.info("aligning unaligned tails to reference")
    tmp_aligned_tails = tmp_dir.name + "/alnigned_tails.bam"
    align_reads_with_minimap(
        reference=reference,
        reads=tmp_fasta_unaligned_tails.name,
        bamout=tmp_aligned_tails,
        tech=tech,
        aln_args=aln_args + " --secondary=no",
        threads=threads,
    )
    # BUG: New alignments refer to shortened reads. Their on-read position is cut off.
    # fix: check the unaligned_intervals for each new alignment and annotate a hard clip to it.
    # -> add_hard_clips_to_aligned_tails

    # merge new and old alignments
    log.info("merging alignments")
    cmd_merge = shlex.split(
        f"samtools merge -f {tmp_dir.name}/merged.bam {alignments} {tmp_dir.name}/aln.bam"
    )
    subprocess.check_call(cmd_merge)
    # output
    log.info("sorting and indexing merged alignments")
    cmd_out = shlex.split(f"samtools sort -o {output} {tmp_dir.name}/merged.bam")
    subprocess.check_call(cmd_out)
    # index
    cmd_index = shlex.split(f"samtools index {output}")
    subprocess.check_call(cmd_index)
    # remove temporary dir
    tmp_dir.cleanup()


# %%


def scale_coords(
    coords: typing.List[int], maxlen: int, line_width: int
) -> typing.List[int]:
    return [int(floor((coord / maxlen) * line_width)) for coord in coords]


def display_ascii_alignments(
    alignments: typing.List[pysam.AlignedSegment],
    terminal_width: int = 100,
    filter_density_min_bp: int = 30,
    filter_density_radius: int = 1000,
    min_bnd_size: int = 300,
    min_signal_size: int = 8,
    min_aln_length: int = 500,
    max_ref_len: int = 0,
):
    if len(alignments) == 0:
        return
    dict_alns = {a.reference_name: [] for a in alignments}
    for refname in dict_alns.keys():
        dict_alns[refname] = [a for a in alignments if a.reference_name == refname]
        alns = dict_alns[refname]
        if max_ref_len < 1:
            try:
                max_ref_len = max(
                    [
                        a.reference_start + a.reference_length
                        for a in alns
                        if a.reference_length
                    ]
                )
            except:
                continue
        print(f"{refname}: {len(alns)} aligned segments, max_ref_len: {max_ref_len}")
        for a in alns:
            if a.reference_length >= min_aln_length:
                raf = alignments_to_rafs.parse_ReadAlignmentFragment_from_alignment(
                    alignment=a,
                    samplename="no samplename",
                    filter_density_min_bp=filter_density_min_bp,
                    filter_density_radius=filter_density_radius,
                    min_bnd_size=min_bnd_size,
                    min_signal_size=min_signal_size,
                )
                # construct the string.
                rstart, rend = scale_coords(
                    [raf.reference_alignment_start, raf.reference_alignment_end],
                    max_ref_len,
                    terminal_width,
                )
                rlen = rend - rstart
                line = (
                    ([" "] * rstart)
                    + (["<" if a.is_reverse else ">"] * rlen)
                    + ([" "] * (terminal_width - rend))
                )
                # insert insertions and deletions
                for sv in raf.SV_signals:
                    if sv.sv_type == 1:  # DEL
                        start, end = scale_coords(
                            [sv.ref_start, sv.ref_end], max_ref_len, terminal_width
                        )
                        # modify the interval on line
                        line[start:end] = ["-"] * (end - start)
                    if sv.sv_type == 0:  # INS
                        start = scale_coords(
                            [sv.ref_start], max_ref_len, terminal_width
                        )[0]
                        line[start] = "V"
                    if sv.sv_type == 4:  # BNDR
                        start = scale_coords(
                            [sv.ref_start], max_ref_len, terminal_width
                        )[0]
                        # check if bndr is hard or soft clipped
                        line[start - 1] = "R" if a.cigartuples[-1][0] == 5 else "r"
                    if sv.sv_type == 3:  # BNDL
                        start = scale_coords([sv.ref_end], max_ref_len, terminal_width)[
                            0
                        ]
                        line[start] = (
                            ("L" if a.cigartuples[0][0] == 5 else "l")
                            if line[start] != "R"
                            else "X"
                        )
                print("".join(line) + " " + a.query_name)


# def display_ascii_rafs(
#         rafs:typing.List[datatypes.ReadAlignmentFragment],
#         terminal_width:int=100,
#         min_indel_size:int=30,
#         min_cliped_size:int=300,
#         min_aln_length:int=500) -> None:
#     if len(rafs) == 0:
#         return
#     dict_rafs = {a.reference_name:[] for a in rafs}
#     for refname in dict_rafs.keys():
#         dict_rafs[refname] = [a for a in rafs if a.reference_name == refname]
#         r = dict_rafs[refname]
#         try:
#             max_ref_len = max([a.referenceAlignmentEnd for a in r])
#         except:
#             continue
#         print(f"{refname}: {len(r)} aligned segments, max_ref_len: {max_ref_len}")
#         for raf in r:
#             if raf.referenceAlignmentEnd - raf.referenceAlignmentStart >= min_aln_length:
#                 # construct the string.
#                 rstart,rend = scale_coords([raf.referenceAlignmentStart,raf.referenceAlignmentEnd],max_ref_len,terminal_width)
#                 rlen = rend-rstart
#                 line = ([' ']*rstart)+(['=']*rlen)+([' ']*(terminal_width-rend))
#                 # insert insertions and deletions
#                 for deletion in raf.deletions:
#                     start,end = scale_coords([deletion.ref_pos,deletion.ref_pos+deletion.ref_length],max_ref_len,terminal_width)
#                     # modify the interval on line
#                     line[start:end] = ['-']*(end-start)
#                 for insertion in raf.insertions:
#                     start = scale_coords([insertion.ref_pos],max_ref_len,terminal_width)[0]
#                     line[start] = 'V'
#                 if raf.breakendRight:
#                     start = scale_coords([raf.breakendRight.referencePosition],max_ref_len,terminal_width)[0]
#                     line[start-1] = '>'
#                 if raf.breakendLeft:
#                     start = scale_coords([raf.breakendLeft.referencePosition],max_ref_len,terminal_width)[0]
#                     clen = raf.breakendLeft.clippedLength
#                     line[start] = '<' if line[start] != '>' else 'X'
#                 print(''.join(line)+raf.read_name)

# =============================================================================


class Direction(Enum):
    NONE = 0
    LEFT = 1
    RIGHT = 2
    BOTH = 3


def get_starts_ends(
    t: np.ndarray, x: np.ndarray, reference_start: int, is_reverse: bool
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # compute a list of tuples of start and end positions for each cigar block
    read_start = x[0] if t[0] in (4, 5) else 0
    length_read = sum(
        x * ((t == 0) | (t == 1) | (t == 4) | (t == 5) | (t == 7) | (t == 8))
    )
    x_read = x * ((t == 0) | (t == 1) | (t == 7) | (t == 8))
    x_read = np.cumsum(x_read) + read_start
    x_read_starts = np.array([read_start, *(x_read)[:-1]])
    x_read_ends = x_read
    #
    x_ref = x * ((t == 0) | (t == 2) | (t == 3) | (t == 7) | (t == 8))
    x_ref = np.cumsum(x_ref) + reference_start
    x_ref_starts = np.array([reference_start, *(x_ref)[:-1]])
    x_ref_ends = x_ref
    if is_reverse:
        return (
            length_read - x_read_ends,
            length_read - x_read_starts,
            x_ref_starts,
            x_ref_ends,
        )
    else:
        return x_read_starts, x_read_ends, x_ref_starts, x_ref_ends


# for forward and backward alignments
def get_ref_pitx_on_read(
    alignment: pysam.AlignedSegment,
    position: int,
    direction: Direction,
    buffer_clipped_length: int = 0,
) -> typing.Tuple[int, int, int]:
    """Get the position on the query sequence that corresponds to a position on the reference. \
Forward or backward orientation of the aligned query sequence is respected.

    Args:
        alignment (pysam.AlignedSegment): the pysam alignment object
        position (int): the position on the reference
        direction (Direction): pick from util.Direction.LEFT, util.Direction.RIGHT, util.Direction.NONE

    Returns:
        typing.Tuple[int,int,int]: the position on the read, the block in which the position is located, the t,x arrays
    """
    is_reverse = alignment.is_reverse
    ref_start = alignment.reference_start
    ref_end = alignment.reference_end
    t, x = zip(*alignment.cigartuples)
    x = np.array(x)
    t = np.array(t)
    aln_length = int(
        sum(x * ((t == 0) | (t == 1) | (t == 4) | (t == 5) | (t == 7) | (t == 8)))
    )
    left_clipped = int(x[0]) if t[0] in (4, 5) else 0
    right_clipped = int(x[-1]) if t[-1] in (4, 5) else 0
    # if the position is left of the read, return alignment.query_alignment_start.
    if position <= ref_start:
        buffer = max(0, left_clipped - buffer_clipped_length)
        if is_reverse:
            return aln_length - buffer, -1, t, x
        else:
            return buffer, 0, t, x
    # if the position is right of the read, return alignment.query_alignment_end.
    if position >= ref_end:
        last_matching_block = np.where(
            (t == 0) | (t == 2) | (t == 3) | (t == 7) | (t == 8)
        )[0][-1]
        buffer = max(0, right_clipped - buffer_clipped_length)
        if is_reverse:
            return buffer, last_matching_block, t, x
        else:
            return aln_length - buffer, last_matching_block, t, x
    x_read_starts, x_read_ends, x_ref_starts, x_ref_ends = get_starts_ends(
        t=t, x=x, reference_start=ref_start, is_reverse=is_reverse
    )
    # find the block in which the position is located. It can be multiple blocks if the position is in an insertion
    # in that case, select fromt he blocks that are not insertions
    block = np.where((x_ref_starts <= position) & (position <= x_ref_ends))[0][0]
    # block = np.where((x_ref_starts <= position) & (position <= x_ref_ends))[0][0]
    # if the block is a match, return the position on the read
    final_pos = -1
    if t[block] in (0, 7, 8):
        if is_reverse:
            final_pos = x_read_starts[block] + (x_ref_ends[block] - position)
        else:
            final_pos = x_read_starts[block] + (position - x_ref_starts[block])
        # final_pos = x_read_starts[block] + \
        #     ((x_ref_ends[block] - position) if is_reverse \
        #         else (position - x_ref_starts[block]))
    # if position is on a deletion, find the next lock that is a match.
    # if direction is left: pick the next left index, where t[index] is a match t[index]==0
    # if direction is right: pick the next right index, where t[index] is a match t[index]==0
    if t[block] == 1:
        # if the position is exactly on an insertion, then go 1 to the right, until 'block' is on a match block,
        # if direction is RIGHT
        # else go 1 to the left.
        if direction == Direction.RIGHT:
            while t[block] not in (0, 7, 8) and block < len(t) - 1:
                block += 1
            final_pos = x_read_starts[block]
        else:
            while t[block] not in (0, 7, 8) and block > 0:
                block -= 1
            final_pos = x_read_ends[block]
    if t[block] in (2, 3):
        if direction == Direction.LEFT:
            # iterate index to the left (-1) until t[index] is a match
            while t[block] not in (0, 7, 8) and block > 0:
                block -= 1
            if is_reverse:
                final_pos = x_read_starts[block]
            else:
                final_pos = x_read_ends[block]
        if direction == Direction.RIGHT:
            # iterate index to the right (+1) until t[index] is a match
            while t[block] not in (0, 7, 8) and block < len(t) - 1:
                block += 1
            if is_reverse:
                final_pos = x_read_ends[block]
            else:
                final_pos = x_read_starts[block]
        if direction == Direction.NONE or direction == Direction.BOTH:
            # find the closest position that is aligned to the left or right of the position.
            block_l, block_r = block, block
            while t[block_l] not in (0, 7, 8) and block_l > 0:
                block_l -= 1
            while t[block_r] not in (0, 7, 8) and block_r < len(t) - 1:
                block_r += 1
            # the final position of block_l is the end of the chosen block
            # the final position of block_r is the start of the chosen block
            pos_l = x_ref_ends[block_l]
            pos_r = x_ref_starts[block_r]
            if position - pos_l < pos_r - position:
                if is_reverse:
                    final_pos = x_read_starts[block_l]
                else:
                    final_pos = x_read_ends[block_l]
                block = block_l
            else:
                if is_reverse:
                    final_pos = x_read_ends[block_r]
                else:
                    final_pos = x_read_starts[block_r]
                block = block_r
    return int(final_pos), int(block), t, x


def get_ref_position_on_read(
    alignment: pysam.AlignedSegment,
    position: int,
    direction: Direction,
    buffer_clipped_length: int = 0,
) -> int:
    return get_ref_pitx_on_read(alignment, position, direction, buffer_clipped_length)[
        0
    ]


# -----------------------------------------------------------------------------


def get_read_pitx_on_ref(
    alignment: pysam.AlignedSegment, position: int, direction: Direction
) -> typing.Tuple[int, int, int]:
    # compute a list of tuples of start and end positions for each cigar block
    is_reverse = alignment.is_reverse
    ref_start = alignment.reference_start
    ref_end = alignment.reference_end
    t, x = zip(*alignment.cigartuples)
    x = np.array(x)
    t = np.array(t)
    read_start = 0
    if is_reverse:
        if t[-1] in (4, 5):
            read_start = x[-1]
    else:
        if t[0] in (4, 5):
            read_start = x[0]
    read_end = read_start + sum(x * ((t == 0) | (t == 1) | (t == 7) | (t == 8)))
    # if the position is on the read, but before the query alignment starts, return alignment.reference_start.
    if position <= read_start:
        if is_reverse:
            last_matching_block = np.where((t == 0) | (t == 1) | (t == 7) | (t == 8))[
                0
            ][
                -1
            ]  # maybe not with 8 (sequence mismatch)
            return int(ref_end), last_matching_block, t, x
        else:
            last_matching_block = 0
            return int(ref_start), last_matching_block, t, x
    # if the position is on the read, but after the query alignment ends, return alignment.reference_end.
    if position >= read_end:
        if is_reverse:
            # the first matching block
            last_matching_block = 0
            return int(ref_start), last_matching_block, t, x
        else:
            last_matching_block = np.where((t == 0) | (t == 1) | (t == 7) | (t == 8))[
                0
            ][
                -1
            ]  # maybe not with 8 (sequence mismatch
            return int(ref_end), last_matching_block, t, x
    x_read_starts, x_read_ends, x_ref_starts, x_ref_ends = get_starts_ends(
        t=t, x=x, is_reverse=is_reverse, reference_start=ref_start
    )
    # print(f"read_starts: {x_read_starts}, read_ends: {x_read_ends}")
    # print(f"ref_starts: {x_ref_starts}, ref_ends: {x_ref_ends}")
    # find the block in which the position is located. It can be multiple blocks if the position is in an insertion
    # in that case, select fromt he blocks that are not insertions
    # TODO: can be sped up with binary search
    block = np.where((x_read_starts <= position) & (position <= x_read_ends))[0][0]
    # block = np.where((x_ref_starts <= position) & (position <= x_ref_ends))[0][0]
    # if the block is a match, return the position on the read
    final_pos = -1
    if t[block] in (0, 7, 8):
        if is_reverse:
            # print(f"final_pos: {x_ref_ends[block]} - {position} + {x_read_ends[block]} = {x_ref_ends[block] - (position - x_read_ends[block])}")
            final_pos = x_ref_ends[block] - (position - x_read_starts[block])
        else:
            # print(f"final_pos: {x_ref_starts[block]} + {position} - {x_read_starts[block]} = {x_ref_starts[block] + position - x_read_starts[block]}")
            final_pos = x_ref_starts[block] + position - x_read_starts[block]
    # if position is on a deletion, find the next lock that is a match.
    # if direction is left: pick the next left index, where t[index] is a match t[index]==0
    # if direction is right: pick the next right index, where t[index] is a match t[index]==0
    if t[block] == 2:
        # if the position is exactly on an insertion, then go 1 to the right, until 'block' is on a match block,
        # if direction is RIGHT
        # else go 1 to the left.
        if direction == Direction.RIGHT:
            while t[block] not in (0, 7, 8) and block < len(t) - 1:
                block += 1
            final_pos = x_ref_starts[block]
        else:
            while t[block] not in (0, 7, 8) and block > 0:
                block -= 1
            final_pos = x_ref_ends[block]
    if t[block] == 1:
        if direction == Direction.LEFT:
            # iterate index to the left (-1) until t[index] is a match
            while t[block] not in (0, 7, 8) and block > 0:
                block -= 1
                if is_reverse:
                    final_pos = x_ref_starts[block]
                else:
                    final_pos = x_ref_ends[block]
        if direction == Direction.RIGHT:
            # iterate index to the right (+1) until t[index] is a match
            while t[block] not in (0, 7, 8) and block < len(t) - 1:
                block += 1
            if is_reverse:
                final_pos = x_ref_ends[block]
            else:
                final_pos = x_ref_starts[block]
        if direction == Direction.NONE:
            # find the closest position that is aligned to the left or right of the position.
            block_l, block_r = block, block
            while t[block_l] not in (0, 7, 8) and block_l > 0:
                block_l -= 1
            while t[block_r] not in (0, 7, 8) and block_r < len(t) - 1:
                block_r += 1
            # the final position of block_l is the end of the chosen block
            # the final position of block_r is the start of the chosen block
            pos_l = x_read_ends[block_l]
            pos_r = x_read_starts[block_r]
            if position - pos_l < pos_r - position:
                if is_reverse:
                    final_pos = x_ref_starts[block_l]
                else:
                    final_pos = x_ref_ends[block_l]
                block = block_l
            else:
                if is_reverse:
                    final_pos = x_ref_ends[block_r]
                else:
                    final_pos = x_ref_starts[block_r]
                block = block_r
    return int(final_pos), int(block), t, x


# start with a position on the read and find the corresponding position on the reference
def get_read_position_on_ref(
    alignment: pysam.AlignedSegment, position: int, direction: Direction
) -> int:
    "start with a position on the read and find the corresponding position on the reference"
    return get_read_pitx_on_ref(alignment, position, direction)[0]


def get_interval_on_consensus_in_region(
    a: pysam.libcalignedsegment.AlignedSegment, start: int, end: int
) -> typing.Tuple[int, int]:
    istart = get_ref_position_on_read(
        alignment=a, position=start, direction=Direction.LEFT
    )
    iend = get_ref_position_on_read(
        alignment=a, position=end, direction=Direction.RIGHT
    )
    return istart, iend


def query_start_end_on_read(aln: pysam.AlignedSegment) -> typing.Tuple[int, int]:
    """This function returns the position of the start and end of the query alignment on the read, following the alignment orientation."""
    t, x = list(zip(*aln.cigartuples))
    t = np.array(t)
    x = np.array(x)
    s = int(x[0]) if t[0] in (4, 5) else int(0)
    body = int(sum(x * ((t == 0) | (t == 1) | (t == 7) | (t == 8))))
    e = int(x[-1]) if t[-1] in (4, 5) else int(0)
    if aln.is_reverse:
        return e, e + body
    else:
        return s, s + body


def query_total_length(aln: pysam.AlignedSegment) -> int:
    """This function returns the total length of the query alignment on the read."""
    t, x = list(zip(*aln.cigartuples))
    t = np.array(t)
    x = np.array(x)
    return int(
        sum(x * ((t == 0) | (t == 1) | (t == 4) | (t == 5) | (t == 7) | (t == 8)))
    )


# %%
# =============================================================================
# function to compute all k-mers from a provided DNA sequence
# given is the alphabet size and alphabet letters.
# the function returns a list of ids of all k-mers in the sequence.
# a k-mer ID is the position in the sorted list of all possible k-mers of length k.


# first, compute a dict to translate any letter from the alphabet to its hash code
def compute_letter_dict(alphabet) -> dict[str, int]:
    if len(alphabet) == 0:
        raise ValueError("alphabet must not be empty")
    if type(alphabet) != str:
        raise ValueError("alphabet must be a string")
    d = {c: i for i, c in enumerate(sorted(list(set(map(str.upper, alphabet)))))}
    d.setdefault("N", len(d))
    return {c: i for i, c in enumerate(sorted(list(set(map(str.upper, alphabet)))))}


def get_hash_of_kmer(kmer: str, letter_dict: dict) -> int:
    if not kmer or len(kmer) == 0:
        raise ValueError("kmer must not be empty")
    n_letters = len(letter_dict)
    k = len(kmer)
    hash_a = sum(
        [
            letter_dict.get(kmer[i], len(letter_dict)) * n_letters ** (k - i - 1)
            for i in range(k)
        ]
    )
    # kmer_rc = str(Seq(kmer).reverse_complement())
    # hash_b = sum([letter_dict.get(kmer_rc[i],len(letter_dict)) * n_letters**(k-i-1) for i in range(k)])
    return hash_a  # min(hash_a,hash_b)


## BEWARE: only computes positive kmers, no rev comp respected
def rolling_hashes(string: str, letter_dict: dict, k: int) -> list[int]:
    if len(string) < k:
        return []
    hashes = [0] * (len(string) - k + 1)
    hashes[0] = get_hash_of_kmer(string[:k], letter_dict)
    first_hash = letter_dict.get(string[0], len(letter_dict))
    N = len(letter_dict)
    # now roll over all letters and always subtract the hash of the first letter (first_hash)
    # and add the hash of the next letter (letter_dict[string[i+k]])
    # to get the new hash value for the next position
    for i in range(1, len(string) - k + 1):
        hashes[i] = int(
            (hashes[i - 1] - first_hash * N ** (k - 1)) * N
            + letter_dict.get(string[i + k - 1], len(letter_dict))
        )
        first_hash = letter_dict.get(string[i], len(letter_dict))
    return hashes


def kmer_counter_from_string(string: str, letter_dict: dict, k: int) -> Counter:
    hashes_a = rolling_hashes(string, letter_dict, k)
    hashes_b = rolling_hashes(str(Seq(string).complement()), letter_dict, k)
    hashes = [min(a, b) for a, b in zip(hashes_a, hashes_b)]
    return Counter(hashes)


# def kmer_set_from_string(string:str,letter_dict:dict,k:int) -> set:
#     hashes_a = rolling_hashes(string,letter_dict,k)
#     hashes_b = rolling_hashes(str(Seq(string).reverse_complement()),letter_dict,k)
#     return set([min(a,b) for a,b in zip(hashes_a,hashes_b)])


def kmer_sketch_from_strings(
    strings: list[str],
    letter_dict: dict,
    k: int,
    min_abs: int = 3,
    min_rel: float = 0.02,
) -> set[int]:
    # get the counter of all k-mers
    counters: list[dict[int, int]] = [
        dict(kmer_counter_from_string(s, letter_dict, k)) for s in strings
    ]
    # sum the values of dicts in counters
    summed_counter_dict: dict[int, int] = dict(
        Counter({k: v for counter in counters for k, v in counter.items()})
    )
    # filter the merged counter dict
    sum_letters = sum(map(len, strings))
    return set(
        [
            k
            for k, v in summed_counter_dict.items()
            if v >= min_abs and v >= min_rel * sum_letters
        ]
    )


def kmer_similarity(
    string_a: str,
    string_b: str,
    letter_dict: dict = {"A": 0, "C": 1, "G": 2, "T": 3},
    k: int = 9,
) -> float:
    a: dict[int, int] = dict(
        kmer_counter_from_string(string=string_a.upper(), letter_dict=letter_dict, k=k)
    )
    b: dict[int, int] = dict(
        kmer_counter_from_string(string=string_b.upper(), letter_dict=letter_dict, k=k)
    )
    if len(a) == 0 or len(b) == 0:
        return 0.0
    # count the number of hits of the smaller string's k-mers in the larger string's k-mers
    smaller, larger = (a, b) if len(a) <= len(b) else (b, a)
    min_hits = sum(
        min(larger.get(kmer, 0), smaller.get(kmer, 0)) for kmer in smaller.keys()
    )
    # calculate the fraction of hits in min_hits over the number of k-mers in the smaller string
    sum_smaller = sum(smaller.values())
    if sum_smaller == 0:
        return 0.0
    return min_hits / sum_smaller


def kmer_similarity_of_groups(
    group_a: list[str],
    group_b: list[str],
    letter_dict: dict = {"A": 0, "C": 1, "G": 2, "T": 3},
    k: int = 9,
) -> float:
    # construct kmer counter dicts for each group and combine them. Combining them means is summing their values
    combined_counter_a = Counter()
    for seq in group_a:
        assert isinstance(seq, str), f"seq {seq} is not a string"
        combined_counter_a += kmer_counter_from_string(
            string=seq, letter_dict=letter_dict, k=k
        )
    a = dict(combined_counter_a)
    combined_counter_b = Counter()
    for seq in group_b:
        assert isinstance(seq, str), f"seq {seq} is not a string"
        combined_counter_b += kmer_counter_from_string(
            string=seq, letter_dict=letter_dict, k=k
        )
    b = dict(combined_counter_b)
    # filter the combined counter by min_abs and min_rel
    smaller, larger = (a, b) if len(a) <= len(b) else (b, a)
    min_hits = sum(
        min(larger.get(kmer, 0), smaller.get(kmer, 0)) for kmer in smaller.keys()
    )
    # calculate the fraction of hits in min_hits over the number of k-mers in the smaller string
    sum_smaller = sum(smaller.values())
    if sum_smaller == 0:
        return 0.0
    return min_hits / sum_smaller


# def jaccard_index(set1:set,set2:set) -> float:
#     return len(set1.intersection(set2)) / len(set1.union(set2))


def intersection_ratio_of_smaller_set(set1: set, set2: set) -> float:
    if len(set1) == 0 or len(set2) == 0:
        return 0.0
    a, b = (set1, set2) if len(set1) <= len(set2) else (set2, set1)
    return len(a.intersection(b)) / len(a)


def similar_regions_on_string(
    start: int, end: int, string: str, max_mismatches: int = 5
) -> list[list[int, int]]:
    """Finds regions that contain substrings similar to the substring at the given position."""
    max_mismatches = max(1, max_mismatches)

    # Input validation
    if not string:
        return []

    if start < 0:
        start = 0
    if end > len(string):
        end = len(string)
    if start > end:
        end, start = start, end  # Ensure start is less than or equal to end
    if start == end:  # No window to check
        return []

    window_size = end - start
    if window_size <= 0:
        return []

    # Additional safety checks
    if start >= len(string) or end > len(string):
        return []

    target = string[start:end]
    if not target:  # Empty target string
        return []

    results = []

    # Slide window across the entire string with bounds checking
    max_start_pos = len(string) - window_size
    if max_start_pos < 0:  # String is shorter than window
        return []

    for i in range(max_start_pos + 1):
        if abs(i - start) < window_size:  # Skip overlapping regions
            continue

        # Bounds check before accessing string indices
        if i + window_size > len(string):
            break

        # Count mismatches with early termination and bounds checking
        mismatches = 0
        for j in range(window_size):
            # Double check bounds (should be safe due to loop condition, but extra safety)
            if i + j >= len(string) or j >= len(target):
                mismatches = max_mismatches + 1  # Force break
                break
            if string[i + j] != target[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break

        if mismatches <= max_mismatches:
            if not results:
                results.append([i, i + window_size])
            elif Interval(i, i + window_size).overlaps(
                Interval(results[-1][0], results[-1][1])
            ):
                # extend last interval
                results[-1][0] = min(results[-1][0], i)
                results[-1][1] = max(results[-1][1], i + window_size)
            else:
                results.append([i, i + window_size])

    return results


# def jaccard_index_of_two_strings(s0:str,s1:str,letter_dict:dict={'A':0,'C':1,'G':2,'T':3},k:int=8) -> float:
#     return jaccard_index(kmer_set_from_string(s0,letter_dict,k),kmer_set_from_string(s1,letter_dict,k))
# %%


def get_read_alignment_intervals_in_region(
    region_start: int,
    regions_end: int,
    alignments: list[pysam.AlignedSegment],
    buffer_clipped_length: int,
) -> dict[str, tuple[int, int, int, int]]:
    # check if all elements in alignments are of type pysam.AlignedSegment
    for aln in alignments:
        assert isinstance(
            aln, pysam.AlignedSegment
        ), f"aln {aln} is not a pysam.AlignedSegment"
    dict_intervals = {}
    for aln in alignments:
        cr_start, cr_end = region_start, regions_end
        read_start, read_end = get_interval_on_read_in_region(
            a=aln,
            start=cr_start,
            end=cr_end,
            buffer_clipped_length=buffer_clipped_length,
        )
        ref_start, ref_end = get_interval_on_ref_in_region(
            a=aln, start=read_start, end=read_end
        )
        # check read_start and read_end if they are numerical and end-start is greater than 0
        if read_start is not None and read_end is not None and read_end > read_start:
            if aln.query_name not in dict_intervals:
                dict_intervals[aln.query_name] = []
            dict_intervals[aln.query_name].append(
                (read_start, read_end, ref_start, ref_end)
            )
            # dict_intervals[read_aug_name].append((min(read_start,read_end),max(read_start,read_end),min(ref_start,ref_end),max(ref_start,ref_end)))
    return dict_intervals


def cut_pysam_alignment_to_region(
    aln: pysam.AlignedSegment, start: int, end: int
) -> pysam.AlignedSegment:
    """cut a pysam AlignedSegment to a region defined by start and end. Returns a new, shortened AlignedSegment."""
    # at first make sure that the alignment overlaps with the region
    if aln.reference_start > end or aln.reference_end < start:
        raise ValueError("alignment does not overlap with region")
    # shorten the cigar string


# cut alignments given a single region
def cut_alignments(
    alignments: list[pysam.AlignedSegment], region: tuple[str, int, int]
) -> list[pysam.AlignedSegment]:
    """cuts alignments such that they seem to come from a read with the max extents of the region"""
    """As a result, all alignments start and end within the region."""
    chr, start, end = region
    # verify that all alignments are on the same chromosome
    assert all(
        [a.reference_name == chr for a in alignments]
    ), "all alignments must be on the same chromosome"
    cut_alns: list[pysam.AlignedSegment] = []
    for alignment in alignments:
        # find the positions on the ref on the read
        start_final_pos, start_block, start_t, start_x = get_ref_pitx_on_read(
            alignment=alignment,
            position=region[1],
            direction=Direction.BOTH,
            buffer_clipped_length=0,
        )
        end_final_pos, end_block, end_t, end_x = get_ref_pitx_on_read(
            alignment=alignment,
            position=region[2],
            direction=Direction.BOTH,
            buffer_clipped_length=0,
        )
        # find the positions of the read on the reference
        ref_start, ref_start_block, ref_start_t, ref_start_x = get_read_pitx_on_ref(
            alignment=alignment, position=start_final_pos, direction=Direction.BOTH
        )
        ref_end, ref_end_block, ref_end_t, ref_end_x = get_read_pitx_on_ref(
            alignment=alignment, position=end_final_pos, direction=Direction.BOTH
        )
        # cut the cigar string. it should be cut in the ref_start_block and ref_end_block
        # but what is the exact offset? i.e. the match block would need to be changed.


def cigartuples_to_cigarstring(cigartuples: list[tuple[int, int]]) -> str:
    code = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X"}
    return "".join([f"{x[1]}{code[x[0]]}" for x in cigartuples])


def cut_alignment(
    alignment: pysam.AlignedSegment, region: tuple[str, int, int]
) -> pysam.AlignedSegment | None:
    # find the positions on the ref on the read
    start_final_pos, start_block, start_t, start_x = get_ref_pitx_on_read(
        alignment=alignment,
        position=region[1],
        direction=Direction.BOTH,
        buffer_clipped_length=0,
    )
    end_final_pos, end_block, end_t, end_x = get_ref_pitx_on_read(
        alignment=alignment,
        position=region[2],
        direction=Direction.BOTH,
        buffer_clipped_length=0,
    )
    # find the positions of the read on the reference
    ref_start, ref_start_block, ref_start_t, ref_start_x = get_read_pitx_on_ref(
        alignment=alignment, position=start_final_pos, direction=Direction.BOTH
    )
    ref_end, ref_end_block, ref_end_t, ref_end_x = get_read_pitx_on_ref(
        alignment=alignment, position=end_final_pos, direction=Direction.BOTH
    )

    x_read_starts, x_read_ends, x_ref_starts, x_ref_ends = get_starts_ends(
        t=start_t,
        x=start_x,
        reference_start=alignment.reference_start,
        is_reverse=alignment.is_reverse,
    )

    if alignment.is_reverse:
        new_start_number = start_x[start_block] - (
            x_read_ends[start_block] - start_final_pos
        )
        new_end_number = x_read_ends[end_block] - end_final_pos
        new_cigar_first = (start_t[start_block], new_start_number)
        new_cigar_last = (end_t[end_block], new_end_number)
        new_cigar = alignment.cigartuples[start_block : end_block + 1]
        if len(new_cigar) == 0:
            return None
        new_cigar[0] = new_cigar_first
        new_cigar[-1] = new_cigar_last
    else:
        new_start_number = start_x[start_block] - (
            start_final_pos - x_read_starts[start_block]
        )
        new_end_number = end_final_pos - x_read_starts[end_block]
        new_cigar_first = (start_t[start_block], new_start_number)
        new_cigar_last = (end_t[end_block], new_end_number)
        new_cigar = alignment.cigartuples[start_block : end_block + 1]
        new_cigar[0] = new_cigar_first
        new_cigar[-1] = new_cigar_last

    new_aln = pysam.AlignedSegment(header=alignment.header)
    new_aln.query_name = alignment.query_name
    new_aln.reference_id = alignment.reference_id
    new_aln.reference_start = ref_start
    # new_aln.reference_end = ref_end
    # new_aln.query_length = abs(end_final_pos - start_final_pos)
    # new_aln.query_alignment_start = 0
    # new_aln.query_alignment_end = abs(end_final_pos - start_final_pos)
    new_aln.cigartuples = new_cigar

    new_aln.cigarstring = cigartuples_to_cigarstring(new_cigar)
    new_aln.flag = alignment.flag
    new_aln.mapping_quality = alignment.mapping_quality

    if alignment.query_sequence:
        if alignment.is_reverse:
            new_aln.query_sequence = alignment.query_sequence[
                end_final_pos : start_final_pos + 1
            ]
        else:
            new_aln.query_sequence = str(
                Seq(
                    alignment.query_sequence[start_final_pos : end_final_pos + 1]
                ).reverse_complement()
            )

    return new_aln


def cut_alignments(
    alignments_file: Path, region: tuple[str, int, int]
) -> list[pysam.AlignedSegment]:
    samfile = pysam.AlignmentFile(alignments_file, "rb")
    return [
        value
        for value in [
            cut_alignment(alignment, region)
            for alignment in samfile.fetch(region[0], region[1], region[2])
        ]
        if value
    ]


def cut_alignments_in_regions(
    samfile: Path, regions: list[tuple[str, int, int]], output: Path
) -> None:
    with pysam.AlignmentFile(
        output, "wb", template=pysam.AlignmentFile(samfile, "rb")
    ) as f:
        for region in tqdm(regions):
            print(region)
            for i, a in enumerate(cut_alignments(samfile, region)):
                # test query lengths
                if a.query_length != a.infer_query_length():
                    print(
                        f"query length mismatch: {a.query_length} vs {a.infer_query_length()}"
                    )
                if abs(a.query_length - a.infer_query_length()) > 1:
                    print(i, a.query_name, a.query_length, a.infer_query_length())
                f.write(a)
    return None



def insert_table_metadata(path_database: str, table_name: str, tag: str):
    """
    Inserts metadata for a table into a hardcoded metadata table in the SQLite database.

    Parameters:
    - path_database (str): Path to the SQLite database file.
    - tag (str): A short label or tag, e.g. 'original_data'.
    """
    metadata_table = "table_metadata"

    with sqlite3.connect(path_database) as conn:
        cursor = conn.cursor()

        # Create metadata table if it doesn't exist
        cursor.execute(
            f"""
            CREATE TABLE IF NOT EXISTS {metadata_table} (
                table_name TEXT PRIMARY KEY,
                tag TEXT );"""
        )

        # Insert or replace metadata for the given table
        cursor.execute(
            f"""
            INSERT OR REPLACE INTO {metadata_table} (table_name, tag)
            VALUES (?, ?);""",
            (table_name, tag),
        )


def get_all_table_metadata(path_database: str) -> dict[str, str]:
    """
    Retrieves all metadata entries from the hardcoded metadata table in the SQLite database.

    Parameters:
    - path_database (str): Path to the SQLite database file.

    Returns:
    - List of tuples containing (table_name, tag)
    """
    metadata_table = "table_metadata"
    with sqlite3.connect(path_database) as conn:
        cursor = conn.cursor()

        # Check if metadata table exists
        cursor.execute(
            """
            SELECT name 
            FROM sqlite_master 
            WHERE type='table' AND name=?;
        """,
            (metadata_table,),
        )

        if not cursor.fetchone():
            log.warning(
                f"Metadata table '{metadata_table}' does not exist in the database."
            )
            return []  # Metadata table doesn't exist yet

        cursor.execute(
            f"""
            SELECT table_name, tag
            FROM {metadata_table}
        """
        )
        return {k: v for k, v in cursor.fetchall()}


# =============================================================================
# Weight functions for distance-based weighting
# =============================================================================


def inverse_power_weight(distance, scale=100.0, falloff=1.0):
    """
    Inverse power law decay function.

    Args:
        distance: Distance from center point
        scale: Distance scaling factor (controls width)
        falloff: Controls falloff shape (1.0=~linear, >1=sub-linear, <1=supra-linear)

    Returns:
        Weight between 0 and 1
    """
    return 1.0 / (1.0 + (distance / scale) ** falloff)


def gaussian_weight(distance, scale=100.0, falloff=2.0):
    """
    Gaussian-based weight function.

    Args:
        distance: Distance from center point
        scale: Width parameter (controls spread)
        falloff: Controls shape (2.0=Gaussian, >2=steeper, <2=flatter)

    Returns:
        Weight between 0 and 1
    """
    return np.exp(-((distance / scale) ** falloff))


def exponential_weight(distance, scale=100.0, falloff=1.0):
    """
    Exponential decay function.

    Args:
        distance: Distance from center point
        scale: Distance scaling factor (controls width)
        falloff: Controls decay rate

    Returns:
        Weight between 0 and 1
    """
    return np.exp(-falloff * distance / scale)


def logistic_weight(distance, scale=100.0, falloff=1.0):
    """
    Logistic-based weight function.

    Args:
        distance: Distance from center point
        scale: Position of transition midpoint
        falloff: Controls transition steepness

    Returns:
        Weight between 0 and 1
    """
    return 1.0 / (1.0 + np.exp(falloff * (distance - scale)))


# we need a function to cut a pysam alignment to a shorter interval
