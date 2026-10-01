"""In-process port of lamassemble: multiple alignment and consensus of long reads.

lamassemble was written by Martin C. Frith and is described in

    Frith MC, Mitsuhashi S, Katoh K. lamassemble: Multiple Alignment and
    Consensus Sequence of Long Reads. In: Katoh K (ed.), Multiple Sequence
    Alignment. Methods Mol Biol 2231:135-145 (2021).
    https://doi.org/10.1007/978-1-0716-1036-7_9
    https://gitlab.com/mcfrith/lamassemble

This module ports lamassemble 1.7.2. The algorithm is unchanged: LAST
(``lastdb``/``lastal``) aligns all sequence pairs, the alignments are laid
out greedily by score, and MAFFT builds the multiple alignment on that guide
tree, restricted by anchors from the pairwise alignments. The consensus is
called column by column. For identical inputs it returns the same consensus as
the ``lamassemble`` command.

What differs is how the tools are run. lamassemble starts a Python
interpreter and then the ``mafft`` driver, a shell script that launches about
60 helper processes (grep, awk, cat, ...) around a single alignment binary.
For the ~1 kb loci svirlpool assembles, that overhead was ~0.25 s of a
~0.3 s call, while the alignment itself (``disttbfast``) took ~10 ms. Here
the steps run in the calling process, and ``disttbfast`` is called directly
with the arguments the driver would pass for lamassemble's options. That
direct call is checked once per process against the ``mafft`` driver on a
small fixture (:func:`direct_mafft_available`); if the results differ, e.g.
under another MAFFT version, the driver is used instead.

Original copyright and licence of lamassemble:

    Copyright 2019 Martin C. Frith
    SPDX-License-Identifier: MIT

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to permit
    persons to whom the Software is furnished to do so, subject to the
    following conditions:

    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
    OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
    NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
    DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
    OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
    USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from __future__ import annotations

import collections
import dataclasses
import functools
import itertools
import logging
import math
import os
import random
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from operator import itemgetter
from pathlib import Path

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class LamassembleParams:
    """lamassemble's options (same names and defaults as its command line)."""

    #: -s: omit consensus flanks with < this many sequences ("N" or "N%")
    seq_min: str = "1"
    #: -g: use alignment columns with <= this % gaps
    gap_max: float = 50.0
    #: --end: count gaps past the ends of the sequences
    end: bool = False
    #: -p: use pairwise restrictions with error probability <= this
    prob: float = 0.002
    #: -d: max change in alignment diagonal between pairwise alignments
    diagonal_max: int = 1000
    #: --all: use all of each sequence, not just the aligning part
    use_all: bool = False
    #: LAST -u: use ~1 per this many initial matches (None: -uNEAR -W)
    u: int | None = None
    #: LAST -W: minimum positions in length-W windows
    W: int = 19
    #: LAST -m: max initial matches per query position
    m: int = 5
    #: LAST -z: max gap length
    z: int = 30
    #: Run lastal on both strands (lamassemble's behaviour). With False only
    #: same-strand alignments are found, so the input must be strand-oriented.
    both_strands: bool = True


class LamassembleTimeout(Exception):
    """The wall-clock budget of one assembly was used up."""


def read_sequences(path: str | Path) -> list[tuple[str, str]]:
    """(name, sequence) of a FASTA or FASTQ file, parsed as lamassemble does."""
    out = []
    record: list[str] = []

    def flush():
        title = record[0]
        name = title[1:].split()[0]
        seq = "".join(record[1:]) if title[0] == ">" else record[1]
        out.append((name, seq))

    with open(path) as f:
        for line in f:
            if (line[0] == ">" and record and record[0][0] == ">") or (
                line[0] == "@" and len(record) == 4
            ):
                flush()
                record = []
            if not line.isspace():
                record.append(line.rstrip())
    if record:
        flush()
    return out


# --------------------------------------------------------------------------
# Score parameters from a last-train file

BUILTIN_TRAIN_FILES = {
    "promethion-2019": """
# scale of score parameters: 4.5512
# delOpenProb: 0.0369615
# insOpenProb: 0.0340916
# delExtendProb: 0.439744
# insExtendProb: 0.403943
# probability matrix (query letters = columns, reference letters = rows):
#   A              C              G              T
# A 0.278908       0.00109169     0.012899       0.00107855
# C 0.00154506     0.197869       0.00044938     0.00272089
# G 0.0272926      0.000508552    0.177789       0.000859018
# T 0.00126293     0.00328421     0.000545484    0.291896
""",
}


@functools.cache
def parameters_from_last_train(train_file: str):
    builtin = BUILTIN_TRAIN_FILES.get(train_file.lower())
    if builtin:
        lines = builtin.splitlines()
    else:
        with open(train_file) as f:
            lines = f.readlines()
    alphabet_size = 4
    scale = del_open = del_extend = ins_open = ins_extend = -1
    prob_matrix: list = list(range(alphabet_size))
    for line in lines:
        fields = line.split()
        if "scale of score parameters" in line:
            scale = float(fields[5])
        if "delOpenProb" in line:
            del_open = float(fields[2])
        if "insOpenProb" in line:
            ins_open = float(fields[2])
        if "delExtendProb" in line:
            del_extend = float(fields[2])
        if "insExtendProb" in line:
            ins_extend = float(fields[2])
        if "probability matrix" in line:
            prob_matrix = []
        if len(fields) == 2 + alphabet_size and len(prob_matrix) < alphabet_size:
            prob_matrix.append([float(i) for i in fields[2:]])
    gap_probs = del_open, del_extend, ins_open, ins_extend
    if scale < 0 or -1 in gap_probs or len(prob_matrix) != alphabet_size:
        raise RuntimeError("can't read the last-train data")
    prob_matrix = _matrix_with_complement_symmetric_rows(prob_matrix)
    return scale, prob_matrix, gap_probs


def _score_from_prob(scale, prob):
    assert prob > 0
    return scale * math.log(prob)


def _cost_from_prob(scale, prob):
    return -_score_from_prob(scale, prob)


def _matrix_with_complement_symmetric_rows(matrix):
    out = []
    for i, j in zip(matrix, reversed(matrix), strict=True):
        avg = (sum(i) + sum(j)) * 0.5
        out.append([k * avg / sum(i) for k in i])
    return out


def _substitution_prob_matrix(ref_matrix, qry_matrix):
    r = range(4)
    s = [sum(row) for row in ref_matrix]
    m = [[0.0 for j in r] for i in r]
    for i in r:
        for j in r:
            m[i][j] = sum(ref_matrix[k][i] * qry_matrix[k][j] / s[k] for k in r)
    return m


def _score_matrix_from_prob_matrix(scale, prob_matrix):
    r = range(len(prob_matrix))
    row_sums = [sum(prob_matrix[i][j] for j in r) for i in r]
    col_sums = [sum(prob_matrix[i][j] for i in r) for j in r]
    out = [[0.0 for j in r] for i in r]
    for i in r:
        for j in r:
            ratio = prob_matrix[i][j] / (row_sums[i] * col_sums[j])
            out[i][j] = _score_from_prob(scale, ratio)
    return out


def _gap_scores_from_probs(scale, gap_probs):
    del_open, del_extend, ins_open, ins_extend = gap_probs
    gap_open = 1 - (1 - del_open) * (1 - ins_open)
    del_frac = del_open / (1 - del_extend)
    ins_frac = ins_open / (1 - ins_extend)
    gap_extend = (del_frac * del_extend + ins_frac * ins_extend) / (
        del_frac + ins_frac
    )
    gap_close = 1 - gap_extend
    first_gap = gap_open * gap_close
    gap_extend += first_gap
    gap_exist = first_gap / gap_extend
    return _cost_from_prob(scale, gap_exist), _cost_from_prob(scale, gap_extend)


@functools.cache
def alignment_scores(train_file: str):
    """(fwd matrix, rev matrix, gap exist cost, gap extend cost)."""
    scale, prob_matrix, gap_probs = parameters_from_last_train(train_file)
    comp_matrix = [row[::-1] for row in reversed(prob_matrix)]
    fwd = _score_matrix_from_prob_matrix(
        scale, _substitution_prob_matrix(prob_matrix, prob_matrix)
    )
    rev = _score_matrix_from_prob_matrix(
        scale, _substitution_prob_matrix(prob_matrix, comp_matrix)
    )
    gap_exist, gap_extend = _gap_scores_from_probs(scale, gap_probs)
    return fwd, rev, gap_exist, gap_extend


# --------------------------------------------------------------------------
# Score parameters in MAFFT format


def _mafft_num(x):
    return format(x, ".3")


def mafft_gap_options(scores) -> list[str]:
    fwd, _rev, gap_exist, gap_extend = scores
    alph_size = len(fwd)
    # lamassemble: assume a, c, g, t have equal frequency
    mean_match = sum(fwd[i][i] for i in range(alph_size)) / alph_size
    mean_unrelated = sum(map(sum, fwd)) / (alph_size**2)
    mean_diff = mean_match - mean_unrelated
    gap_exist_mafft = gap_exist / mean_diff
    gap_extend_mafft = (gap_extend * 2 + mean_unrelated) / mean_diff
    op = _mafft_num(gap_exist_mafft)
    gop = _mafft_num(-gap_exist_mafft)
    ep = _mafft_num(gap_extend_mafft)
    gep = _mafft_num(-gap_extend_mafft)
    lep = _mafft_num(-mean_unrelated / mean_diff)
    lexp = _mafft_num(-gap_extend / mean_diff)
    return [
        "--op", op, "--gop", gop, "--ep", ep, "--gep", gep, "--gexp", "0",
        "--lop", gop, "--lep", lep, "--lexp", lexp,
    ]  # fmt: skip


def _mafft_matrix_text(fwd, rev) -> str:
    alphabet = "ARNDCQEGHILKMFPSTWYV"
    d = {"A": (0, 0), "C": (0, 1), "G": (0, 2), "T": (0, 3),  # + strand bases
         "E": (1, 0), "Q": (1, 1), "F": (1, 2), "P": (1, 3)}  # - strand bases  # fmt: skip
    lines = [" ".join([" "] + [format(i, ">5") for i in alphabet])]
    for i, x in enumerate(alphabet):
        row = x
        for y in alphabet[: i + 1]:
            score = 0.0
            if x in d and y in d:
                xm, xi = d[x]
                ym, yi = d[y]
                if xm == 0 and ym == 0:
                    score = fwd[xi][yi]
                if xm == 0 and ym == 1:
                    score = rev[xi][yi]
                if xm == 1 and ym == 0:
                    score = rev[3 - xi][3 - yi]
                if xm == 1 and ym == 1:
                    score = fwd[3 - xi][3 - yi]
            row += " " + format(score, "5.3")
        lines.append(row)
    return "\n".join(lines) + "\n"


def mafft_params_text(scores) -> str:
    fwd, rev, _, _ = scores
    alphabet = "ARNDCQEGHILKMFPSTWYV"
    freqs = [(0.0, 0.25)[i in "ACGT"] for i in alphabet]
    return (
        " ".join(["# MAFFT cost:"] + mafft_gap_options(scores))
        + "\n"
        + _mafft_matrix_text(fwd, rev)
        + " ".join(["frequency"] + [str(f) for f in freqs])
        + "\n"
    )


def _last_score_text(matrix, gap_exist, gap_extend) -> str:
    def to_int(x):
        return int(round(x))

    lines = [f"#last -a {to_int(gap_exist)} -b {to_int(gap_extend)}"]
    lines.append("  A C G T")
    for i, row in zip("ACGT", matrix, strict=True):
        lines.append(" ".join([i] + [str(to_int(v)) for v in row]))
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Consensus from aligned sequences


def _alignment_columns_for_consensus(params: LamassembleParams, rows):
    alignment_length = max(map(len, rows))
    seq_num_change = [0] * (alignment_length + 1)
    for i in rows:
        beg = len(i) - len(i.lstrip("-"))
        end = len(i.rstrip("-"))
        seq_num_change[beg] += 1
        seq_num_change[end] -= 1
    if params.seq_min.endswith("%"):
        min_seqs = len(rows) * (float(params.seq_min[:-1]) / 100)
    else:
        min_seqs = int(params.seq_min)
    beg = end = alignment_length
    seq_num = 0
    for i, x in enumerate(seq_num_change):
        seq_num += x
        if seq_num >= min_seqs:
            if beg > i:
                beg = i
            end = i + 1
    columns = zip(*rows, strict=False)
    seq_num = len(rows) if params.end else 0
    for i, (x, y) in enumerate(zip(seq_num_change, columns, strict=False)):
        if not params.end:
            seq_num += x
        if beg <= i < end:
            gap_num = seq_num + y.count("-") - len(y)
            if gap_num * 100 <= seq_num * params.gap_max:
                yield "".join(y)


def _strand_score(score_row, column, bases):
    return sum(score_row[i] * column.count(x) for i, x in enumerate(bases))


def _column_score(prior_scores, score_matrix, column, fwd_base_index):
    rev_base_index = 3 - fwd_base_index
    fwd_score = _strand_score(score_matrix[fwd_base_index], column, "ACGT")
    rev_score = _strand_score(score_matrix[rev_base_index], column, "tgca")
    return prior_scores[fwd_base_index] + fwd_score + rev_score


def _consensus_col(prior_scores, score_matrix, column):
    column = column.replace("U", "T")
    scores = [_column_score(prior_scores, score_matrix, column, i) for i in range(4)]
    j, _m = max(enumerate(scores), key=itemgetter(1))
    return "acgt"[j]


def consensus_sequence(
    params: LamassembleParams, prob_matrix, aligned_rows: list[str]
) -> str:
    """lamassemble's consensus of aligned rows (upper case: + strand; lower: -)."""
    if not aligned_rows:
        raise RuntimeError("can't make a consensus of zero sequences")
    base_probs = [sum(row) for row in prob_matrix]
    prior_scores = [math.log(i) for i in base_probs]
    score_matrix = [[math.log(x / sum(row)) for x in row] for row in prob_matrix]
    # Columns repeat a lot; the consensus base depends on the column only.
    cache: dict[str, str] = {}
    out = []
    for col in _alignment_columns_for_consensus(params, aligned_rows):
        base = cache.get(col)
        if base is None:
            base = _consensus_col(prior_scores, score_matrix, col)
            cache[col] = base
        out.append(base)
    return "".join(out)


# --------------------------------------------------------------------------
# Pairwise alignments and their layout

_X_CRAZY = "ACGT" "RYKMBDHV" "U"
_Y_CRAZY = "PFQE" "YRMKVHDB" "E"
_CRAZY_TABLE = str.maketrans(_X_CRAZY, _Y_CRAZY)
_UNCRAZY_TABLE = str.maketrans("EQFP", "acgt")


def _crazy_reverse_complement(seq: str) -> str:
    # Reverse-strand bases go to MAFFT as the "amino acids" E, Q, F, P.
    return seq[::-1].translate(_CRAZY_TABLE)


def _maf_alignments(lines):
    """lamassemble's alignmentInput: (-score, qry, ref, probCodes) per alignment."""
    score = 0
    seq_records: list = []
    for line in lines:
        if not line or line[0] == "#":
            continue
        fields = line.split()
        if line[0] == "a":
            for i in fields:
                if i.startswith("score="):
                    score = int(i[6:])
                    seq_records = []
        elif line[0] == "s":
            seq_num = int(fields[1])
            seq_beg = int(fields[2])
            seq_end = seq_beg + int(fields[3])
            strand = fields[4]
            seq_len = int(fields[5])
            seq_records.append(
                (seq_num, seq_len, strand, seq_beg, seq_end, fields[6])
            )
        elif line[0] == "p":
            ref, qry = seq_records
            if ref[0] != qry[0]:
                yield -score, qry, ref, fields[1]


def _find_range(old_ranges, new_range):
    qry_beg, qry_end, ref_beg, ref_end = new_range
    for i, (qb, qe, rb, re) in enumerate(old_ranges):
        if qry_beg < qb and qry_end < qe and ref_beg < rb and ref_end < re:
            return i
        if qry_beg <= qb or qry_end <= qe or ref_beg <= rb or ref_end <= re:
            return -1
    return len(old_ranges)


def _is_big_diagonal_change(params, old_ranges, new_range, idx):
    qry_beg, qry_end, ref_beg, ref_end = new_range
    if idx > 0:
        qb, qe, rb, re = old_ranges[idx - 1]
        if abs((re - qe) - (ref_beg - qry_beg)) > params.diagonal_max:
            return True
    if idx < len(old_ranges):
        qb, qe, rb, re = old_ranges[idx]
        if abs((rb - qb) - (ref_end - qry_end)) > params.diagonal_max:
            return True
    return False


def _layout_of_seqs(params, num_seqs, alignments_sorted_by_score):
    data_per_seq = [(i, False) for i in range(num_seqs)]
    alignment_order = []
    kept_alignments = []
    aligned_ranges = collections.defaultdict(list)
    for aln in alignments_sorted_by_score:
        _neg_score, qry, ref, _prob_codes = aln
        if qry[0] > ref[0]:
            qry, ref = ref, qry
        ref_num, ref_len, ref_strand, ref_beg, ref_end, _ = ref
        qry_num, qry_len, qry_strand, qry_beg, qry_end, _ = qry
        ref_group, is_rev_ref = data_per_seq[ref_num]
        qry_group, is_rev_qry = data_per_seq[qry_num]
        is_opposite = is_rev_qry != is_rev_ref
        is_rev_aln = qry_strand != ref_strand
        is_flip = is_rev_aln != is_opposite
        old_ranges = aligned_ranges[qry_num, ref_num]
        if ref_strand == "-":
            ref_beg, ref_end = ref_len - ref_end, ref_len - ref_beg
            qry_beg, qry_end = qry_len - qry_end, qry_len - qry_beg
        new_range = qry_beg, qry_end, ref_beg, ref_end
        rpos = _find_range(old_ranges, new_range)
        if ref_group == qry_group:
            if (
                is_flip
                or rpos < 0
                or _is_big_diagonal_change(params, old_ranges, new_range, rpos)
            ):
                continue
        else:
            min_group = min(ref_group, qry_group)
            max_group = max(ref_group, qry_group)
            alignment_order.append((min_group, max_group))
            for i, (group, is_rev) in enumerate(data_per_seq):
                if group == max_group:
                    data_per_seq[i] = min_group, is_rev != is_flip
        kept_alignments.append(aln)
        old_ranges.insert(rpos, new_range)
    return data_per_seq, alignment_order, kept_alignments


def _update_range(beg_per_seq, end_per_seq, is_rev_per_seq, rec):
    seq_num, seq_len, seq_strand, seq_beg, seq_end, _ = rec
    if (seq_strand == "-") != is_rev_per_seq[seq_num]:
        seq_beg, seq_end = seq_len - seq_end, seq_len - seq_beg
    beg_per_seq[seq_num] = min(beg_per_seq[seq_num], seq_beg)
    end_per_seq[seq_num] = max(end_per_seq[seq_num], seq_end)


def _aligned_range_per_seq(kept_alignments, is_rev_per_seq):
    beg = [2**63 - 1] * len(is_rev_per_seq)  # sys.maxsize in lamassemble
    end = [0] * len(is_rev_per_seq)
    for _neg, qry, ref, _p in kept_alignments:
        _update_range(beg, end, is_rev_per_seq, ref)
        _update_range(beg, end, is_rev_per_seq, qry)
    return beg, end


def _pairwise_anchors(params, min_prob_code, alns, seq_ranks, is_rev, beg_per_seq):
    for _neg, qry, ref, prob_codes in alns:
        if qry[0] > ref[0]:
            qry, ref = ref, qry
        ref_num, ref_len, ref_strand, ref_beg, ref_end, ref_aln = ref
        qry_num, qry_len, qry_strand, qry_beg, qry_end, qry_aln = qry
        ref_rank = seq_ranks[ref_num]
        qry_rank = seq_ranks[qry_num]
        if ref_rank < 1 or qry_rank < 1:
            continue
        if is_rev[ref_num] != (ref_strand == "-"):
            ref_aln = ref_aln[::-1]
            qry_aln = qry_aln[::-1]
            prob_codes = prob_codes[::-1]
            ref_beg = ref_len - ref_end
            qry_beg = qry_len - qry_end
        if not params.use_all:
            ref_beg -= beg_per_seq[ref_num]
            qry_beg -= beg_per_seq[qry_num]
        q_beg = r_beg = size = 0
        for p, q, r in zip(prob_codes, qry_aln, ref_aln, strict=False):
            if q != "-":
                if r != "-":
                    if p >= min_prob_code:
                        if q_beg + size < qry_beg or r_beg + size < ref_beg:
                            if size:
                                yield qry_rank, ref_rank, q_beg, r_beg, size
                            q_beg = qry_beg
                            r_beg = ref_beg
                            size = 0
                        size += 1
                    ref_beg += 1
                qry_beg += 1
            elif r != "-":
                ref_beg += 1
        if size:
            yield qry_rank, ref_rank, q_beg, r_beg, size


# --------------------------------------------------------------------------
# Running the tools


class _Deadline:
    def __init__(self, timeout: float | None):
        self.end = None if timeout is None else time.monotonic() + timeout
        self.timeout = timeout

    def remaining(self) -> float | None:
        if self.end is None:
            return None
        left = self.end - time.monotonic()
        if left <= 0:
            raise LamassembleTimeout(f"timed out after {self.timeout} s")
        return left

    def run(self, cmd, **kw) -> subprocess.CompletedProcess:
        try:
            return subprocess.run(cmd, check=True, timeout=self.remaining(), **kw)
        except subprocess.TimeoutExpired as e:
            raise LamassembleTimeout(f"timed out after {self.timeout} s") from e


def _pairwise_alignments(params, scores, sequences, tmpdir, threads, deadline):
    fwd, rev, gap_exist, gap_extend = scores
    fwd_mat = os.path.join(tmpdir, "fwd.mat")
    rev_mat = os.path.join(tmpdir, "rev.mat")
    with open(fwd_mat, "w") as f:
        f.write(_last_score_text(fwd, gap_exist, gap_extend))
    with open(rev_mat, "w") as f:
        f.write(_last_score_text(rev, gap_exist, gap_extend))
    flipped = [False] * len(sequences)
    if not params.both_strands:
        flipped = orient_by_kmers([seq for _name, seq in sequences])
    seq_file = os.path.join(tmpdir, "x.fa")
    with open(seq_file, "w") as f:
        for i, (_name, seq) in enumerate(sequences):
            f.write(f">{i}\n{_reverse_complement(seq) if flipped[i] else seq}\n")
    db = os.path.join(tmpdir, "db")
    opt_p = f"-P{threads}"
    if params.u:
        cmd = ["lastdb", "-c", opt_p, f"-uRY{params.u}", db, seq_file]
    else:
        cmd = ["lastdb", "-c", opt_p, "-uNEAR", f"-W{params.W}", db, seq_file]
    deadline.run(cmd)
    alignments = []
    strands = ((1, fwd_mat), (0, rev_mat)) if params.both_strands else ((1, fwd_mat),)
    for strand, mat in strands:
        cmd = [
            "lastal", f"-s{strand}", "-j4", "-D1e9", opt_p, f"-m{params.m}",
            f"-z{params.z}g", "-p", mat, db, seq_file,
        ]  # fmt: skip
        out = deadline.run(cmd, stdout=subprocess.PIPE, text=True).stdout
        alignments.extend(_maf_alignments(out.splitlines()))
    if any(flipped):
        alignments = [_unflip(aln, flipped) for aln in alignments]
    alignments.sort()
    return alignments


# --------------------------------------------------------------------------
# One-strand mode (not in lamassemble): orient the reads first, then align
# them on the forward strand only.

_COMPLEMENT = str.maketrans("ACGTUMRWSYKVHDBNacgtumrwsykvhdbn", "TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _kmers(seq: str, k: int) -> set[str]:
    seq = seq.upper()
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def orient_by_kmers(seqs: list[str], k: int = 13) -> list[bool]:
    """Which sequences to reverse-complement so that all share one strand.

    Greedy, longest first: each sequence takes the orientation that shares more
    k-mers with the sequences oriented before it. A sequence that shares none
    either way keeps its orientation (it would not link to the others anyway).
    """
    order = sorted(range(len(seqs)), key=lambda i: -len(seqs[i]))
    flipped = [False] * len(seqs)
    seen: set[str] = set()
    for i in order:
        fwd = _kmers(seqs[i], k)
        rev = _kmers(_reverse_complement(seqs[i]), k)
        if len(rev & seen) > len(fwd & seen):
            flipped[i] = True
            fwd = rev
        seen |= fwd
    return flipped


def _unflip(aln, flipped):
    """Alignment of oriented copies -> the same alignment of the original reads.

    A '+' row of a reverse-complemented copy is the '-' strand of the original
    read, with the same coordinates and aligned text (MAF counts '-' strand
    coordinates on the reverse complement).
    """
    neg_score, qry, ref, prob_codes = aln

    def fix(rec):
        if not flipped[rec[0]]:
            return rec
        num, length, strand, beg, end, text = rec
        return num, length, "-" if strand == "+" else "+", beg, end, text

    return neg_score, fix(qry), fix(ref), prob_codes


def _awk_number(x: float) -> str:
    """Format a number as awk's ``print`` does (OFMT/CONVFMT %.6g)."""
    if x == int(x) and abs(x) < 1e16:
        return str(int(x))
    return format(x, ".6g")


def disttbfast_args(gap_options: list[str], num_seqs: int) -> list[str]:
    """The disttbfast arguments of ``mafft --amino --quiet --aamatrix M
    --treein T --anchors A <gap_options>``, as the MAFFT 7.526 driver builds
    them (single-threaded, FFT-NS-2, no reordering).

    Only --op and --ep reach disttbfast; --gop/--gep/--gexp/--lop/--lep/--lexp
    are parameters of other MAFFT strategies. The driver does one cycle
    instead of two for exactly two sequences.
    """
    opts = dict(zip(gap_options[::2], gap_options[1::2], strict=True))
    gop = opts["--op"]
    aof = _awk_number(0.0 + float(_awk_number(-1.0 * float(opts["--ep"]))))
    return [
        "-q", "0", "-E", "1" if num_seqs == 2 else "2", "-V", "-" + gop, "-s", "0.0", "-W", "6", "-O",
        "-C", "0-0", "-U", "-P", "-b", "-1", "-g", "0", "-f", "-" + gop,
        "-Q", "100.0", "-h", aof, "-F", "-X", "0.1", "-l", "-x", "100",
    ]  # fmt: skip


@functools.cache
def _disttbfast_path() -> str | None:
    candidates = []
    if os.environ.get("MAFFT_BINARIES"):
        candidates.append(Path(os.environ["MAFFT_BINARIES"]) / "disttbfast")
    mafft = shutil.which("mafft")
    if mafft:
        prefix = Path(os.path.realpath(mafft)).parent.parent
        candidates.append(prefix / "libexec" / "mafft" / "disttbfast")
    for c in candidates:
        if c.is_file() and os.access(c, os.X_OK):
            return str(c)
    return None


def _write_mafft_inputs(tmpdir, fasta_text, tree_text, anchors_text, mtx_text):
    paths = {}
    for name, text in (
        ("fa", fasta_text),
        ("tree", tree_text),
        ("pair", anchors_text),
        ("mtx", mtx_text),
    ):
        paths[name] = os.path.join(tmpdir, "mafft." + name)
        with open(paths[name], "w") as f:
            f.write(text)
    return paths


def _run_mafft_driver(paths, gap_options, deadline) -> str:
    cmd = ["mafft", "--amino", *gap_options, "--quiet"]
    cmd += ["--aamatrix", paths["mtx"], "--treein", paths["tree"]]
    cmd += ["--anchors", paths["pair"], paths["fa"]]
    return deadline.run(cmd, stdout=subprocess.PIPE, text=True).stdout


def _run_disttbfast(binary, tmpdir, fasta_text, tree_text, anchors_text,
                    mtx_text, gap_options, deadline) -> str:  # fmt: skip
    """Stage the files as the mafft driver does and run disttbfast on them."""
    work = os.path.join(tmpdir, "mafft_direct")
    os.mkdir(work)

    def nonblank(text):
        return "".join(line + "\n" for line in text.splitlines() if line)

    with open(os.path.join(work, "infile"), "w") as f:
        f.write(fasta_text + "\n")
    with open(os.path.join(work, "_aamtx"), "w") as f:
        f.write(nonblank(mtx_text))
    with open(os.path.join(work, "_guidetree"), "w") as f:
        f.write(nonblank(tree_text))
    with open(os.path.join(work, "_externalanchors"), "w") as f:
        f.write(nonblank(anchors_text))
    with open(os.path.join(work, "infile")) as fin:
        return deadline.run(
            [binary, *disttbfast_args(gap_options, fasta_text.count(">"))],
            stdin=fin,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            cwd=work,
            text=True,
        ).stdout


def _aligned_rows(fasta_text: str) -> list[str]:
    """Rows of a MAFFT FASTA alignment, - strand rows back in lower case acgt."""
    rows: list[str] = []
    is_rev = False
    for line in fasta_text.splitlines():
        if not line:
            continue
        if line[0] == ">":
            is_rev = line.startswith(">_R_")
            rows.append("")
        elif is_rev:
            rows[-1] += line.translate(_UNCRAZY_TABLE).lower()
        else:
            rows[-1] += line
    return rows


def _mafft_alignment(tmpdir, fasta_text, tree_text, anchors_text, scores,
                     deadline, direct: bool) -> list[str]:  # fmt: skip
    gap_options = mafft_gap_options(scores)
    mtx_text = mafft_params_text(scores)
    binary = _disttbfast_path() if direct else None
    if binary is not None:
        out = _run_disttbfast(binary, tmpdir, fasta_text, tree_text,
                              anchors_text, mtx_text, gap_options, deadline)  # fmt: skip
    else:
        paths = _write_mafft_inputs(tmpdir, fasta_text, tree_text, anchors_text, mtx_text)
        out = _run_mafft_driver(paths, gap_options, deadline)
    return _aligned_rows(out)


def multiple_alignment(
    sequences: list[tuple[str, str]],
    train_file: str,
    params: LamassembleParams,
    threads: int = 1,
    timeout: float | None = None,
    tmp_dir: str | Path | None = None,
    direct: bool | None = None,
) -> list[str]:
    """lamassemble's multiple alignment of the largest linked group of sequences.

    Returns the aligned rows (MAFFT order; reverse-strand rows in lower case).
    ``direct`` None: call disttbfast directly if :func:`direct_mafft_available`.
    """
    if not sequences:
        raise RuntimeError("zero sequences")
    ok = set("ACGTUMRWSYKVHDBN-acgtumrwsykvhdbn")
    for _name, seq in sequences:
        diff = set(seq).difference(ok)
        if diff:
            raise RuntimeError(
                "not allowed in nucleotide sequences: " + "".join(sorted(diff))
            )
    sequences = [(name, seq.replace("-", "")) for name, seq in sequences]
    if direct is None:
        direct = direct_mafft_available()
    deadline = _Deadline(timeout)
    scores = alignment_scores(str(train_file))
    with tempfile.TemporaryDirectory(prefix="lamassemble", dir=tmp_dir) as tmpdir:
        alignments = _pairwise_alignments(
            params, scores, sequences, tmpdir, threads, deadline
        )
        deadline.remaining()
        data_per_seq, alignment_order, kept = _layout_of_seqs(
            params, len(sequences), alignments
        )
        if not params.both_strands and len(alignment_order) < len(sequences) - 1:
            # One strand left reads unlinked. In tandem repeats, orienting all
            # reads alike multiplies each seed's matches past LAST's -m limit,
            # which reads on the other strand stay under: redo with both.
            log.debug("lamassemble: one-strand layout incomplete, using both strands")
            params = dataclasses.replace(params, both_strands=True)
            alignments = _pairwise_alignments(
                params, scores, sequences, tmpdir, threads, deadline
            )
            deadline.remaining()
            data_per_seq, alignment_order, kept = _layout_of_seqs(
                params, len(sequences), alignments
            )
        group_per_seq, is_rev_per_seq = zip(*data_per_seq, strict=True)
        seqs_per_group = [0] * len(sequences)
        for g in group_per_seq:
            seqs_per_group[g] += 1
        max_seqs = max(seqs_per_group)
        if max_seqs < len(sequences):
            log.debug(
                f"lamassemble: using {max_seqs} out of {len(sequences)} sequences"
                " (linked by pairwise alignments)"
            )
        kept_group = seqs_per_group.index(max_seqs)
        seq_ranks = []
        rank = 0
        for g in group_per_seq:
            if g == kept_group:
                rank += 1  # the first rank is 1, not 0, for MAFFT
                seq_ranks.append(rank)
            else:
                seq_ranks.append(0)
        beg_per_seq, end_per_seq = _aligned_range_per_seq(kept, is_rev_per_seq)

        fasta = []
        for seq_rank, (name, seq), is_rev, beg, end in zip(
            seq_ranks, sequences, is_rev_per_seq, beg_per_seq, end_per_seq,
            strict=True,
        ):  # fmt: skip
            if seq_rank < 1:
                continue
            seq = seq.upper()
            if is_rev:
                name = "_R_" + name
                seq = _crazy_reverse_complement(seq)
            else:
                seq = seq.replace("U", "T")
            if not params.use_all:
                if beg >= end:
                    beg = end = 0
                seq = seq[beg:end]
            fasta.append(f">{name}\n{seq}\n")
        fasta_text = "".join(fasta)

        tree = []
        for i, j in alignment_order:
            if seq_ranks[i] > 0 and seq_ranks[j] > 0:
                tree.append(f"{seq_ranks[i]} {seq_ranks[j]} 0.1 0.1\n")
        tree_text = "".join(tree)

        min_prob_score = int(math.ceil(-10 * math.log10(params.prob)))
        min_prob_code = chr(33 + min_prob_score)
        anchors = []
        for neg_score, alns in itertools.groupby(kept, itemgetter(0)):
            p = _pairwise_anchors(
                params, min_prob_code, alns, seq_ranks, is_rev_per_seq, beg_per_seq
            )
            for k, _ in itertools.groupby(sorted(p)):
                qry_rank, ref_rank, q_beg, r_beg, size = k
                anchors.append(
                    f"{qry_rank} {ref_rank} {q_beg + 1} {q_beg + size} "
                    f"{r_beg + 1} {r_beg + size} {-neg_score}\n"
                )
        anchors_text = "".join(anchors)
        deadline.remaining()
        return _mafft_alignment(
            tmpdir, fasta_text, tree_text, anchors_text, scores, deadline, direct
        )


def assemble(
    sequences: list[tuple[str, str]],
    train_file: str | Path,
    params: LamassembleParams,
    threads: int = 1,
    timeout: float | None = None,
    tmp_dir: str | Path | None = None,
    direct: bool | None = None,
) -> str:
    """Consensus sequence of ``sequences`` = ``lamassemble`` with ``params``.

    Raises LamassembleTimeout when ``timeout`` (seconds, whole call) runs out,
    subprocess.CalledProcessError when a tool fails, RuntimeError on bad input.
    """
    rows = multiple_alignment(
        sequences, str(train_file), params, threads, timeout, tmp_dir, direct
    )
    _scale, prob_matrix, _gap_probs = parameters_from_last_train(str(train_file))
    return consensus_sequence(params, prob_matrix, rows)


# --------------------------------------------------------------------------
# One-time check of the direct disttbfast call against the mafft driver


def _self_check_sequences() -> list[tuple[str, str]]:
    """A few noisy copies (both strands) of a random sequence with an insertion."""
    rng = random.Random(7)
    base = "".join(rng.choice("ACGT") for _ in range(400))
    ins = "".join(rng.choice("ACGT") for _ in range(60))
    comp = str.maketrans("ACGT", "TGCA")
    seqs = []
    for i in range(6):
        s = base[:200] + (ins if i % 2 else "") + base[200:]
        out = []
        for c in s:
            r = rng.random()
            if r < 0.03:
                continue
            if r < 0.06:
                out.append(rng.choice("ACGT"))
            if r < 0.08:
                out.append(rng.choice("ACGT"))
            out.append(c)
        s = "".join(out)[rng.randrange(20) :]
        if i % 3 == 2:
            s = s[::-1].translate(comp)
        seqs.append((f"s{i}", s))
    return seqs


@functools.cache
def direct_mafft_available() -> bool:
    """True if disttbfast, called directly, reproduces the mafft driver here."""
    if _disttbfast_path() is None:
        log.warning("lamassemble: disttbfast not found; using the mafft driver")
        return False
    seqs = _self_check_sequences()
    params = LamassembleParams(seq_min="2", gap_max=67, m=50)
    try:
        direct = multiple_alignment(seqs, "promethion-2019", params, direct=True)
        driver = multiple_alignment(seqs, "promethion-2019", params, direct=False)
    except (OSError, subprocess.CalledProcessError, RuntimeError) as e:
        log.warning(f"lamassemble: direct MAFFT self-check failed ({e}); using the mafft driver")
        return False
    if direct != driver:
        log.warning(
            "lamassemble: disttbfast called directly differs from the mafft "
            "driver (another MAFFT version?); using the mafft driver"
        )
        return False
    return True
