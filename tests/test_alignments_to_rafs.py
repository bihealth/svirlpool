import json
import tempfile
from gzip import open as gzopen
from pathlib import Path

import cattrs
from pysam import AlignedSegment, AlignmentFile

from svirlpool.localassembly.consensus import \
    parse_ReadAlignmentSignals_from_alignment
from svirlpool.signalprocessing.alignments_to_rafs import get_start_end
from svirlpool.util import util
from svirlpool.util.datatypes import Alignment, ReadAlignmentSignals, SVsignal

DATA_DIR = Path(__file__).parent / "data" / "signalprocessing"


# This is used to generate the test data. It's not part of the test run, it has to be run manually!
def generate_test_data() -> None:
    # --- INV for signal parsing tests --- #
    path_alignments = DATA_DIR / "alignments_to_rafs.dummy_inversion.json.gz"
    reference: list[str] = util.generate_sequence(size=6_000, seed=1)
    reference[1996:2004] = tuple("ACGTACGT")
    reference[3996:4004] = tuple("TGCATGCA")
    read = (
        reference[:2000]
        + util.reverse_complement(reference[2000:4000])
        + reference[4000:]
    )
    tmp_alignment = tempfile.NamedTemporaryFile(delete=True, suffix=".bam")
    alignments: list[Alignment] = [
        Alignment.from_pysam(aln)
        for aln in util.create_alignments_to_reference(
            reads=[read], reference=reference, alignments=Path(tmp_alignment.name)
        )
    ]
    with gzopen(path_alignments, "wt") as f:
        json.dump(
            {"alignments": [aln.unstructure() for aln in alignments]},
            f,
            indent=4,
        )


# This is used in the tests
def load_alignments(path: Path) -> list[Alignment]:
    with gzopen(path, "rt") as f:
        data = json.load(f)
    alignments = cattrs.structure(data["alignments"], list[Alignment])
    return alignments


def test_get_start_end() -> None:
    # load alignments
    alignments: list[Alignment] = load_alignments(
        DATA_DIR / "alignments_to_rafs.dummy_inversion.json.gz"
    )
    ref_start, ref_end, read_start, read_end = get_start_end(alignments[0].to_pysam())
    assert (0, 2000, 0, 2000) == (ref_start, ref_end, read_start, read_end)
    ref_start, ref_end, read_start, read_end = get_start_end(alignments[1].to_pysam())
    assert (2000, 4000, 2000, 4000) == (ref_start, ref_end, read_start, read_end)
    ref_start, ref_end, read_start, read_end = get_start_end(alignments[2].to_pysam())
    assert (4000, 6000, 4000, 6000) == (ref_start, ref_end, read_start, read_end)


# def test_parse_SVsignals_from_alignment() -> None:
#     # load alignments
#     alignments:list[Alignment] = [a.to_pysam() for a in load_alignments(
#         DATA_DIR / "alignments_to_rafs.dummy_inversion.json.gz"
#     )]
#     results = [parse_SVsignals_from_alignment(
#         aln, *get_start_end(aln), 10, 1000) for aln in alignments]

#     expected = [
#         [SVsignal(ref_start=2000, ref_end=2000, read_start=2000, read_end=2000, size=4000, sv_type=4)],
#         [SVsignal(ref_start=2000, ref_end=2000, read_start=4000, read_end=4000, size=2000, sv_type=3), SVsignal(ref_start=4000, ref_end=4000, read_start=2000, read_end=2000, size=2000, sv_type=4)],
#         [SVsignal(ref_start=4000, ref_end=4000, read_start=4000, read_end=4000, size=4000, sv_type=3)]]
#     for res, exp in zip(results, expected):
#         assert res == exp


# this tests parse_SVsignals_from_alignment and get_start_end
def test_consensus_parse_ReadAlignmentSignals_from_alignment() -> None:
    # dummy data
    alignments: list[AlignedSegment] = [
        a.to_pysam()
        for a in load_alignments(
            DATA_DIR / "alignments_to_rafs.dummy_inversion.json.gz"
        )
    ]
    results = [
        parse_ReadAlignmentSignals_from_alignment(
            samplename="test", alignment=aln, min_signal_size=10, min_bnd_size=1000
        )
        for aln in alignments
    ]
    expected = [
        ReadAlignmentSignals(
            samplename="test",
            read_name="read-0",
            reference_name="ref",
            alignment_forward=True,
            SV_signals=[
                SVsignal(
                    ref_start=2000,
                    ref_end=2000,
                    read_start=2000,
                    read_end=2000,
                    size=4000,
                    sv_type=4,
                )
            ],
        ),
        ReadAlignmentSignals(
            samplename="test",
            read_name="read-0",
            reference_name="ref",
            alignment_forward=False,
            SV_signals=[
                SVsignal(
                    ref_start=2000,
                    ref_end=2000,
                    read_start=4000,
                    read_end=4000,
                    size=2000,
                    sv_type=3,
                ),
                SVsignal(
                    ref_start=4000,
                    ref_end=4000,
                    read_start=2000,
                    read_end=2000,
                    size=2000,
                    sv_type=4,
                ),
            ],
        ),
        ReadAlignmentSignals(
            samplename="test",
            read_name="read-0",
            reference_name="ref",
            alignment_forward=True,
            SV_signals=[
                SVsignal(
                    ref_start=4000,
                    ref_end=4000,
                    read_start=4000,
                    read_end=4000,
                    size=4000,
                    sv_type=3,
                )
            ],
        ),
    ]
    for res, exp in zip(results, expected, strict=True):
        assert res == exp
