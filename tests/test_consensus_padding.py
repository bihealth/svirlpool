"""The padded consensus continues its padding reads across both junctions.

create_padding_for_consensus joins the flank of one read to each end of the
consensus. The flank must begin exactly where that read's cut read stops
aligning to the consensus; otherwise the padded consensus duplicates or skips
sequence at the core boundary, and its alignment to the reference shows an
insertion (a tandem copy of the adjacent reference, or foreign sequence) or a
deletion right at the core start or end -- a false SV call.

The synthetic reads here are error-free pieces of a random reference, so a
correct padded consensus of a reference-like locus is an exact substring of the
reference.
"""

import random

import pysam
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from svirlpool.localassembly import consensus_class
from svirlpool.localassembly.consensus import (
    PaddingFlank,
    aligned_interval_on_read,
    choose_padding_flank,
    create_padding_for_consensus,
    read_anchoring,
    trim_consensus_inplace,
)

REF_LEN = 80_000
CORE = (40_000, 42_000)  # the cut extent of the reads on the reference


def random_seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


REF = random_seq(REF_LEN, 1)


def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())


class Reads:
    """Full reads and their cut reads, as process_consensus_container has them."""

    def __init__(self) -> None:
        self.read_records: dict[str, SeqRecord] = {}
        self.cutreads: dict[str, SeqRecord] = {}
        self.intervals: list[tuple[int, int, str, bool]] = []

    def add(
        self,
        name: str,
        seq: str,
        cut: tuple[int, int],
        forward: bool = True,
        interval: tuple[int, int] = (0, 0),
    ) -> None:
        """`seq` in the read's original orientation; `cut` on it; `forward`:
        the read's strand relative to the consensus (as the cut reads align)."""
        self.read_records[name] = SeqRecord(
            Seq(seq), id=name, name=name, description=""
        )
        start, end = cut
        self.cutreads[name] = SeqRecord(
            Seq(seq[start:end]),
            id=name,
            name=name,
            description=(
                f"crID=0,start={start},end={end},ref_start_chr=chr1,ref_start=0,"
                f"ref_end_chr=chr1,ref_end=0"
            ),
        )
        self.intervals.append((*interval, name, forward))

    def add_ref_read(
        self, name: str, ref_start: int, ref_end: int, forward: bool = True
    ) -> None:
        """A read of the reference [ref_start, ref_end), cut to CORE."""
        piece = REF[ref_start:ref_end]
        cut_ref = (max(CORE[0], ref_start), min(CORE[1], ref_end))
        if forward:
            seq = piece
            cut = (cut_ref[0] - ref_start, cut_ref[1] - ref_start)
        else:
            seq = revcomp(piece)
            cut = (ref_end - cut_ref[1], ref_end - cut_ref[0])
        self.add(name, seq, cut, forward)


def make_consensus(seq: str, reads: Reads) -> consensus_class.Consensus:
    return consensus_class.Consensus(
        ID="0.0",
        crIDs=[0],
        original_regions=[("chr1", *CORE)],
        consensus_sequence=seq,
        intervals_cutread_alignments=[
            (s, e if e > s else len(seq), rn, fw) for s, e, rn, fw in reads.intervals
        ],
        clustering_meta_data={},
    )


def pad(cons: consensus_class.Consensus, reads: Reads, anchoring=None):
    return create_padding_for_consensus(
        consensus_object=cons,
        cutreads=reads.cutreads,
        read_records=reads.read_records,
        max_padding_size=100_000,
        anchoring=anchoring,
    )


def assert_padded_is_reference(padding, cons) -> tuple[int, int]:
    """The padded consensus is one piece of the reference; returns where the
    core lies on the reference."""
    padded = padding.sequence.upper()
    pos = REF.find(padded)
    assert pos >= 0, "padded consensus is not contiguous reference sequence"
    assert REF.find(padded, pos + 1) < 0
    start, end = padding.consensus_interval_on_sequence_with_padding
    assert padded[start:end] == cons.consensus_sequence.upper()
    assert padding.padding_size_left == start
    assert padding.padding_size_right == len(padded) - end
    return pos + start, pos + end


def spanning_reads(reads: Reads, n: int = 6) -> None:
    """n reads crossing CORE with different flank lengths, alternating strand."""
    for i in range(n):
        left = 3_000 + 1_500 * i
        right = 12_000 - 1_200 * i
        reads.add_ref_read(
            f"span{i}", CORE[0] - left, CORE[1] + right, forward=i % 2 == 0
        )


# --- whole-function behaviour ------------------------------------------------------


def test_reference_locus_forward_and_reverse_reads():
    reads = Reads()
    spanning_reads(reads)
    cons = make_consensus(REF[CORE[0] : CORE[1]], reads)
    padding = pad(cons, reads)
    core = assert_padded_is_reference(padding, cons)
    assert core == CORE
    # the longest flank on each side: span5 (left, 10.5 kb), span0 (right, 12 kb)
    assert padding.readname_left == "span5"
    assert padding.readname_right == "span0"
    assert padding.padding_size_left == 3_000 + 1_500 * 5
    assert padding.padding_size_right == 12_000


@pytest.mark.parametrize("forward", [True, False])
def test_only_reverse_or_only_forward_padding_reads(forward):
    reads = Reads()
    reads.add_ref_read("long", CORE[0] - 20_000, CORE[1] + 20_000, forward=forward)
    reads.add_ref_read("short", CORE[0] - 2_000, CORE[1] + 2_000, forward=not forward)
    cons = make_consensus(REF[CORE[0] : CORE[1]], reads)
    padding = pad(cons, reads)
    assert assert_padded_is_reference(padding, cons) == CORE
    assert padding.readname_left == padding.readname_right == "long"


def test_consensus_end_reached_by_unanchored_reads_is_trimmed():
    """25130.0-like: two reads are cut 1.7 kb further (their reference
    alignment ends inside the cut; the cut continues into the clipped rest);
    the consensus follows them, the reads with the longest flank stop at
    CORE[1]. Joining such a read's flank to the consensus end (the old code)
    duplicates 1.7 kb."""
    reads = Reads()
    spanning_reads(reads)
    for i in range(2):
        name = f"long_cut{i}"
        seq = REF[CORE[0] - 1_000 : CORE[1] + 3_000]
        reads.add(name, seq, (1_000, 1_000 + (CORE[1] - CORE[0]) + 1_700))
    cons = make_consensus(REF[CORE[0] : CORE[1] + 1_700], reads)
    anchoring = dict.fromkeys(reads.read_records, (True, True))
    anchoring["long_cut0"] = anchoring["long_cut1"] = (True, False)
    padding = pad(cons, reads, anchoring=anchoring)
    core = assert_padded_is_reference(padding, cons)
    assert core == CORE
    assert cons.clustering_meta_data["padding_trim_right"] == 1_700
    assert cons.clustering_meta_data["padding_trim_left"] == 0
    # the cut-read intervals now refer to the trimmed core
    assert all(
        0 <= s < e <= len(cons.consensus_sequence)
        for s, e, _, _ in cons.intervals_cutread_alignments
    )


def test_foreign_consensus_tail_from_chimeric_reads_is_trimmed():
    """6779.0-like: two chimeric reads continue with foreign sequence (another
    chromosome) at CORE[1]; their clipped part is cut with them and the
    consensus ends in 464 bp of it. They have the longest flanks."""
    foreign = random_seq(20_000, 7)
    reads = Reads()
    spanning_reads(reads)
    for i in range(2):
        seq = REF[CORE[0] - 5_000 : CORE[1]] + foreign[: 15_000 + 1_000 * i]
        reads.add(f"chimeric{i}", seq, (5_000, 5_000 + CORE[1] - CORE[0] + 464))
    cons = make_consensus(REF[CORE[0] : CORE[1]] + foreign[:464], reads)
    anchoring = dict.fromkeys(reads.read_records, (True, True))
    anchoring["chimeric0"] = anchoring["chimeric1"] = (True, False)
    padding = pad(cons, reads, anchoring=anchoring)
    assert assert_padded_is_reference(padding, cons) == CORE
    assert cons.clustering_meta_data["padding_trim_right"] == 464
    assert padding.readname_right == "span0"


def test_clipped_sequence_at_the_consensus_start_is_trimmed():
    """22967.1-like: three reads end (physically) in 141 bp of clipped
    sequence before CORE[0]; the consensus starts with it."""
    junk = random_seq(141, 11)
    reads = Reads()
    spanning_reads(reads)
    for i in range(3):
        seq = junk + REF[CORE[0] : CORE[1] + 4_000 + 1_000 * i]
        cut = (0, 141 + CORE[1] - CORE[0])
        if i == 1:  # sequenced on the other strand
            seq = revcomp(seq)
            cut = (len(seq) - cut[1], len(seq) - cut[0])
        reads.add(f"clipped{i}", seq, cut, forward=i != 1)
    cons = make_consensus(junk + REF[CORE[0] : CORE[1]], reads)
    padding = pad(cons, reads)
    assert assert_padded_is_reference(padding, cons) == CORE
    assert cons.clustering_meta_data["padding_trim_left"] == 141


def test_foldback_read_does_not_pad():
    """5042.1-like: a foldback read (the locus, then its reverse complement)
    has the longest flanks but aligns to the consensus on both strands."""
    reads = Reads()
    spanning_reads(reads)
    fwd = REF[CORE[0] - 30_000 : CORE[1] + 300]
    seq = (
        fwd
        + revcomp(REF[CORE[0] - 300 : CORE[1] + 300])
        + revcomp(REF[CORE[0] - 30_000 : CORE[0] - 300])
    )
    # its cut spans both copies of the locus: read [30000, 32000) is the
    # locus, read [32600, 34600) its reverse complement
    reads.add("foldback", seq, (30_000, len(fwd) + 300 + (CORE[1] - CORE[0])))
    cons = make_consensus(REF[CORE[0] : CORE[1]], reads)
    padding = pad(cons, reads)
    assert assert_padded_is_reference(padding, cons) == CORE
    assert "foldback" not in (padding.readname_left, padding.readname_right)


def test_reads_ending_inside_the_core():
    """Reads that start or end inside the core have no flank on that side;
    the longest reads of all are such reads."""
    reads = Reads()
    spanning_reads(reads, n=3)
    reads.add_ref_read("ends_inside", CORE[0] - 30_000, CORE[0] + 900, forward=True)
    reads.add_ref_read("starts_inside", CORE[1] - 700, CORE[1] + 30_000, forward=False)
    cons = make_consensus(REF[CORE[0] : CORE[1]], reads)
    padding = pad(cons, reads)
    assert assert_padded_is_reference(padding, cons) == CORE
    assert padding.readname_left == "ends_inside"
    assert padding.readname_right == "starts_inside"
    assert padding.padding_size_left == 30_000
    assert padding.padding_size_right == 30_000


def test_soft_clipped_unanchored_flank_is_avoided():
    """A read whose reference alignment ends at the cut (soft clip: the rest
    of the read is foreign) has the longest flank, but an anchored read pads."""
    foreign = random_seq(40_000, 13)
    reads = Reads()
    spanning_reads(reads)
    seq = REF[CORE[0] - 4_000 : CORE[1]] + foreign
    reads.add("clipped", seq, (4_000, 4_000 + CORE[1] - CORE[0]))
    cons = make_consensus(REF[CORE[0] : CORE[1]], reads)
    anchoring = dict.fromkeys(reads.read_records, (True, True))
    anchoring["clipped"] = (True, False)
    padding = pad(cons, reads, anchoring=anchoring)
    assert assert_padded_is_reference(padding, cons) == CORE
    assert padding.readname_right != "clipped"


def test_consensus_shorter_than_the_cut_reads():
    """lamassemble may end the consensus inside the reads' cut: the flank
    then starts with the read's unaligned cut bases, without loss."""
    reads = Reads()
    spanning_reads(reads)
    cons = make_consensus(REF[CORE[0] + 20 : CORE[1] - 30], reads)
    padding = pad(cons, reads)
    assert assert_padded_is_reference(padding, cons) == (CORE[0] + 20, CORE[1] - 30)
    assert cons.clustering_meta_data == {}


def test_consensus_with_a_deletion_and_reads_of_both_haplotypes():
    """The consensus of the haplotype with a 300 bp deletion, padded by reads
    of either haplotype: only the deletion separates it from the reference."""
    m = CORE[0] + 900
    alt = REF[CORE[0] : m] + REF[m + 300 : CORE[1]]
    reads = Reads()
    # haplotype with the deletion (shorter flanks)
    for i in range(3):
        left, right = 3_000 + 500 * i, 4_000 + 500 * i
        seq = REF[CORE[0] - left : m] + REF[m + 300 : CORE[1] + right]
        forward = i % 2 == 0
        cut = (left, left + len(alt))
        if not forward:
            seq = revcomp(seq)
            cut = (len(seq) - cut[1], len(seq) - cut[0])
        reads.add(f"del{i}", seq, cut, forward=forward)
    # reference haplotype reads with the longest flanks (aligned to this consensus
    # with the deletion in their alignment)
    reads.add_ref_read("ref_hap_left", CORE[0] - 25_000, CORE[1] + 2_000, forward=False)
    reads.add_ref_read("ref_hap_right", CORE[0] - 2_000, CORE[1] + 25_000, forward=True)
    cons = make_consensus(alt, reads)
    padding = pad(cons, reads)
    padded = padding.sequence.upper()
    start, end = padding.consensus_interval_on_sequence_with_padding
    assert padded[start:end] == alt
    assert padding.readname_left == "ref_hap_left"
    assert padding.readname_right == "ref_hap_right"
    assert padded == REF[CORE[0] - 25_000 : m] + REF[m + 300 : CORE[1] + 25_000]


def test_no_flanking_read_leaves_the_side_unpadded():
    reads = Reads()
    reads.add("exact", REF[CORE[0] : CORE[1]], (0, CORE[1] - CORE[0]))
    cons = make_consensus(REF[CORE[0] : CORE[1]], reads)
    padding = pad(cons, reads)
    assert padding.sequence == REF[CORE[0] : CORE[1]]
    assert padding.consensus_interval_on_sequence_with_padding == (0, CORE[1] - CORE[0])
    assert padding.readname_left == padding.readname_right == ""


def test_padding_is_capped_next_to_the_core():
    reads = Reads()
    spanning_reads(reads)
    cons = make_consensus(REF[CORE[0] : CORE[1]], reads)
    padding = create_padding_for_consensus(
        consensus_object=cons,
        cutreads=reads.cutreads,
        read_records=reads.read_records,
        max_padding_size=5_000,
    )
    core = assert_padded_is_reference(padding, cons)
    assert core == CORE
    assert padding.padding_size_left == padding.padding_size_right == 5_000


# --- the pieces -----------------------------------------------------------------------


def flank(name, junction, length, anchored=True):
    return PaddingFlank(name, junction, 0, length, True, anchored)


def test_choose_padding_flank_takes_the_outermost_anchored_junction():
    flanks = [flank(f"r{i}", 900, 5_000 + i) for i in range(5)]
    tails = [
        flank("tail_a", 1_364, 50_000, anchored=False),
        flank("tail_b", 1_364, 60_000, anchored=False),
    ]
    assert choose_padding_flank(flanks + tails, "right").readname == "r4"
    # one anchored read reaching farther decides the junction
    far = flank("far", 1_364, 2_000)
    assert choose_padding_flank([*flanks, far], "right").readname == "far"
    flanks = [flank(f"l{i}", 141, 5_000 + i) for i in range(5)]
    head = flank("head", 0, 9_000, anchored=False)
    assert choose_padding_flank([*flanks, head], "left").readname == "l4"
    head = flank("head", 0, 1_200)
    assert choose_padding_flank([*flanks, head], "left").readname == "head"
    # within the tolerance of the outermost junction the longest flank pads
    flanks = [flank("end", 2_070, 3_000), flank("near_end", 2_045, 36_000)]
    assert choose_padding_flank(flanks, "right").readname == "near_end"


def test_choose_padding_flank_without_anchored_flanks_takes_the_median():
    flanks = [flank(f"r{i}", 900, 5_000 + i, anchored=False) for i in range(3)]
    flanks += [
        flank("tail_a", 1_364, 50_000, anchored=False),
        flank("tail_b", 1_364, 60_000, anchored=False),
    ]
    assert choose_padding_flank(flanks, "right").readname == "r2"
    # on a tie the junction farther inside the consensus
    flanks = [
        flank("a", 900, 5_000, anchored=False),
        flank("b", 1_400, 9_000, anchored=False),
    ]
    assert choose_padding_flank(flanks, "right").readname == "a"
    assert choose_padding_flank(flanks, "left").readname == "b"


def test_choose_padding_flank_prefers_anchored_and_long_flanks():
    flanks = [
        flank("unanchored", 900, 90_000, anchored=False),
        flank("anchored", 905, 3_000),
    ]
    assert choose_padding_flank(flanks, "right").readname == "anchored"
    # flanks shorter than PADDING_MIN_FLANK only count when no longer one exists
    flanks = [
        flank("short", 2_000, 300),
        flank("long1", 900, 4_000),
        flank("long2", 905, 3_000),
    ]
    assert choose_padding_flank(flanks, "right").readname == "long1"
    assert (
        choose_padding_flank([flank("short", 2_000, 300)], "right").readname == "short"
    )
    assert choose_padding_flank([], "right") is None
    assert choose_padding_flank([flank("empty", 10, 0)], "left") is None
    with pytest.raises(ValueError):
        choose_padding_flank(flanks, "up")


def test_trim_consensus_inplace():
    cons = consensus_class.Consensus(
        ID="1.0",
        crIDs=[1],
        original_regions=[("chr1", 0, 10)],
        consensus_sequence="AACCGGTTAA",
        intervals_cutread_alignments=[
            (0, 10, "a", True),
            (0, 2, "b", False),
            (5, 9, "c", True),
        ],
    )
    trim_consensus_inplace(cons, 2, 8)
    assert cons.consensus_sequence == "CCGGTT"
    assert cons.intervals_cutread_alignments == [(0, 6, "a", True), (3, 6, "c", True)]
    with pytest.raises(ValueError):
        trim_consensus_inplace(cons, 4, 4)


def make_aln(name, cigar, reverse, query_len):
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 100_000}]})
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.reference_id = 0
    a.reference_start = 1_000
    a.cigarstring = cigar
    a.is_reverse = reverse
    a.query_sequence = "A" * query_len
    return a


def test_aligned_interval_on_read():
    assert aligned_interval_on_read(make_aln("r", "100S800M100S", False, 1_000)) == (
        100,
        900,
    )
    assert aligned_interval_on_read(make_aln("r", "100S850M50S", True, 1_000)) == (
        50,
        900,
    )


def test_read_anchoring():
    # read of 10 kb; alignment covers read [0, 6000) forward; cut [2000, 5000)
    alns = {0: [make_aln("a", "6000M4000S", False, 10_000)]}
    reads = Reads()
    reads.add("a", "A" * 10_000, (2_000, 5_000))
    assert read_anchoring(alns, reads.cutreads, reads.read_records) == {
        "a": (True, True)
    }
    # cut ends where the alignment ends: the rest of the read is a clip
    reads = Reads()
    reads.add("a", "A" * 10_000, (2_000, 6_000))
    assert read_anchoring(alns, reads.cutreads, reads.read_records) == {
        "a": (True, False)
    }
    # reverse alignment covering read [4000, 10000) in the original orientation
    alns = {0: [make_aln("a", "6000M4000S", True, 10_000)]}
    reads = Reads()
    reads.add("a", "A" * 10_000, (4_000, 8_000))
    assert read_anchoring(alns, reads.cutreads, reads.read_records) == {
        "a": (False, True)
    }
