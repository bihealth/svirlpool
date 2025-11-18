# %%

import tempfile
from pathlib import Path

import pysam

from svirlpool.localassembly import consensus
from svirlpool.util import util

DATA_DIR = Path(__file__).resolve().parent / "data" / "read_cutting"

# %%

# =============================================================================
#  test read cutting
# =============================================================================


def test_trim_reads_multiple():
    filepath = DATA_DIR / "platinum.bam"
    try:
        alns = list(pysam.AlignmentFile(filepath, "rb"))
    except Exception as e:
        print(f"Error opening BAM: {e}")
        print(f"pysam version: {pysam.__version__}")
        raise

    region_start = 21112695
    regions_end = 21122675
    dict_all_intervals = consensus.get_read_alignment_intervals_in_region(
        alignments=alns,
        buffer_clipped_length=0,
        region_start=region_start,
        regions_end=regions_end,
    )
    max_intervals = consensus.get_max_extents_of_read_alignments_on_cr(
        dict_all_intervals=dict_all_intervals
    )
    read_records = consensus.get_full_read_sequences_of_alignments(
        dict_alignments={0: alns}, path_alignments=DATA_DIR / "platinum.bam"
    )
    cut_reads = consensus.trim_reads(
        dict_alignments={0: alns},
        intervals=max_intervals,
        read_records=read_records,
    )
    result = {readname: len(rc) for readname, rc in cut_reads.items()}

    expected = {
        "06f1d59d-2559-49a0-ac27-3c7bc8dfa723": 9936,
        "408fe495-893f-4cfd-bf2f-15953e3bf41f": 9982,
        "e32c0759-77d1-4ae9-8a17-f8ad228d2f53": 9893,
        "4bcc7456-8df6-4472-9c64-cc82ca642b38": 9900,
        "e98b61e1-f273-49a1-b74e-f8f965a890fd": 9885,
        "4e977c76-e639-4887-8b1c-71ce98432603": 9925,
        "1f07a9a1-7bf0-4120-b47d-9d1e54b47412": 9904,
        "6c40fb61-ad17-4e63-ac9d-026a1c82dab1": 9891,
        "d69878f8-4ee1-41f0-b472-a952b7b8dcd6": 9900,
        "277483bd-edbc-4ca4-bf40-f9e5be0cb63b": 9893,
        "a513e5f1-1011-4a3d-a3b4-b9c184e9f083": 9869,
        "226ef612-42f1-4b00-85fa-5800824eb44e": 9898,
        "6b6fcdda-997a-4406-9832-3782abe54504": 9946,
        "4809cafa-6202-4981-a810-a752c96b183b": 9900,
        "796ad9b0-8325-4774-96c8-8dce86526b25": 9860,
        "21985f85-24eb-4461-bab7-98f2f70516d8": 9862,
        "44824ca7-4ef5-42e7-9756-3275157850ac": 9929,
        "d9d47dfd-f913-475b-b3f0-24a5e9be045c": 9905,
        "886c55ec-a679-40af-9d8f-58e43f4c8a7a": 5747,
        "7bfd65fa-58b5-4ade-8565-fd3bb08fe4db": 2099,
        "143d83e8-b375-449d-8564-7c05fea7b07f": 196,
    }
    for key in expected:
        assert result[key] == expected[key]


def test_get_interval_on_read_in_region__multiple():
    aln = None
    for a in pysam.AlignmentFile(DATA_DIR / "platinum.bam"):
        if a.query_name == "886c55ec-a679-40af-9d8f-58e43f4c8a7a":
            aln = a
            break
    region_start = 21112695
    regions_end = 21122675
    result = util.get_interval_on_read_in_region(
        a=aln, start=region_start, end=regions_end
    )
    expected = (93, 5840)
    assert result == expected


# test cutting reads that have an open end
def test_cut_reads_open_end_fwd():
    # create ref and 4 reads that have a large insertion at the end
    with tempfile.TemporaryDirectory() as tdir:
        ref = util.generate_sequence(10_000, seed=31)
        reads = [
            ref[:6000] + util.generate_sequence(4000, seed=137),
            ref[:6000] + util.generate_sequence(2000, seed=138),
        ]
        tmp_alignments = tempfile.NamedTemporaryFile(
            delete=True, suffix=".bam", dir=tdir
        )
        alns = util.create_alignments_to_reference(
            reference=ref,
            reads=reads,
            alignments=Path(tmp_alignments.name),
            aln_args=" -Y --secondary=no --sam-hit-only",
        )
        # util.display_ascii_alignments(alns)
        region_start = 1000
        regions_end = 8000
        dict_all_intervals = consensus.get_read_alignment_intervals_in_region(
            alignments=alns,
            buffer_clipped_length=4000,
            region_start=region_start,
            regions_end=regions_end,
        )
        max_intervals = consensus.get_max_extents_of_read_alignments_on_cr(
            dict_all_intervals=dict_all_intervals
        )
        read_records = consensus.get_full_read_sequences_of_alignments(
            dict_alignments={0: alns}, path_alignments=tmp_alignments.name
        )
        cut_reads = consensus.trim_reads(
            dict_alignments={0: alns},
            intervals=max_intervals,
            read_records=read_records,
        )
        # check lengths of cut reads
        assert len(cut_reads) == 2
        assert len(cut_reads["read-0"]) == 9_000
        assert len(cut_reads["read-1"]) == 7_000


# now do the same with reversed reads
def test_cut_reads_open_end_rvs():
    with tempfile.TemporaryDirectory() as tdir:
        ref = util.generate_sequence(10_000, seed=31)
        reads = [
            util.reverse_complement(ref[6000:]) + util.generate_sequence(6000, 137),
            util.reverse_complement(ref[6000:]) + util.generate_sequence(3000, 138),
        ]
        tmp_alignments = tempfile.NamedTemporaryFile(
            delete=True, suffix=".bam", dir=tdir
        )
        alns = util.create_alignments_to_reference(
            reference=ref,
            reads=reads,
            alignments=Path(tmp_alignments.name),
            aln_args=" -Y --secondary=no --sam-hit-only",
        )
        region_start = 1000
        regions_end = 8000
        buffer_clipped_length = 4000
        dict_all_intervals = consensus.get_read_alignment_intervals_in_region(
            alignments=alns,
            buffer_clipped_length=buffer_clipped_length,
            region_start=region_start,
            regions_end=regions_end,
        )
        max_intervals = consensus.get_max_extents_of_read_alignments_on_cr(
            dict_all_intervals=dict_all_intervals
        )
        read_records = consensus.get_full_read_sequences_of_alignments(
            dict_alignments={0: alns}, path_alignments=tmp_alignments.name
        )
        cut_reads = consensus.trim_reads(
            dict_alignments={0: alns},
            intervals=max_intervals,
            read_records=read_records,
        )
        # check lengths of cut reads
        assert len(cut_reads) == 2
        assert len(cut_reads["read-0"]) == 6000
        assert len(cut_reads["read-1"]) == 5000


def test_get_ref_pitx_on_read__real_long_rev():
    alns = DATA_DIR / "long_rev.bam"
    aln = next(pysam.AlignmentFile(alns, "rb"))
    positions = [819_386, 932_891]
    expecteds = [113338, 68]
    for position, expected in zip(positions, expecteds, strict=True):
        result = util.get_ref_pitx_on_read(
            alignment=aln, position=position, direction=util.Direction.NONE
        )[0]
        assert result == expected


def test_get_interval_on_read_in_region__real_long_rev():
    alns = DATA_DIR / "long_rev.bam"
    aln = next(pysam.AlignmentFile(alns, "rb"))
    regions = [(920481, 927769), (865112, 876018)]
    expectets = [(5322, 12566), (58324, 68267)]
    for region, expected in zip(regions, expectets, strict=True):
        interval = util.get_interval_on_read_in_region(
            a=aln, start=region[0], end=region[1]
        )
        assert interval == expected


def test_get_interval_on_read_in_region__real_long_rev2():
    alns = DATA_DIR / "long_rev.2.bam"
    aln = next(pysam.AlignmentFile(alns, "rb"))
    regions = [(920481, 927769), (865112, 876018)]
    expectets = [(31, 6300), (51939, 61888)]
    for region, expected in zip(regions, expectets, strict=True):
        interval = util.get_interval_on_read_in_region(
            a=aln, start=region[0], end=region[1]
        )
        assert interval == expected


# =============================================================================
#  forward read tests
# =============================================================================


# %%


def test_util_get_interval_on_read_in_region__simple():
    # tests a read as long as the reference on the region (2000,6000).
    # ==|====|==== #
    # --|----|---- #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = (0, reflength)
    expected_read_region = 2000, 6000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__longer_clipped_right():
    # tests a read that is longer than the region (2000,6000) and is clipped on the right side.
    # ==|====|==== #
    #  -|---#|###  # 4000M 4000S
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 1000, 5000
    reference: list = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    read.extend(
        util.generate_sequence(4000, seed=137)
    )  # append more random DNA to the right
    expected_read_region = (
        1000,
        4000,
    )  # left is the position on the read, where the region starts, right is the min(region, read).
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__longer_clipped_left():
    # tests a read that is longer than the region (2000,6000) and is clipped on the left side.
    # ==|====|==== #
    #  #|#---|---  #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 3000, 9000
    reference: list = util.generate_sequence(reflength, seed=0)
    read = util.generate_sequence(2000, seed=137)
    read.extend(reference.copy()[read_subsequence[0] : read_subsequence[1]])
    expected_read_region = (
        2000,
        5000,
    )  # left is the start of the read, right is the min(read,region).
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__shorter():
    # tests a read that is shorter than the ref, but covers the region (2000,6000).
    # ==|====|==== #
    #  -|----|---  #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 1000, 9000
    expected_read_region = 1000, 5000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__half_overlap():
    # tests a read that is shorter than the ref and covers the region only partially.
    # ==|====|==== #
    #     ---|--   #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 3000, 8000
    expected_read_region = 0, 3000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__full_underlap():
    # tests a read that covers only the middle of the region.
    # ==|====|==== #
    #     --       #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 3000, 5000
    expected_read_region = 0, 2000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__clipped_start():
    # tests a read that is shorter than the ref, but covers the region (2000,6000).
    # ==|====|==== #
    #  x|x---|---  #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 1000, 9000
    expected_read_region = 2000, 5000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    # replace first 2kb with rnd dna
    read[:2000] = util.generate_sequence(2000)
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


# =============================================================================
#  reverse read tests
# =============================================================================


def test_util_get_interval_on_read_in_region__reverse_simple():
    # tests a read as long as the reference on the region (2000,6000).
    # ==|====|==== #
    # --|----|---- #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = (0, reflength)
    expected_read_region = 4000, 8000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    read = util.reverse_complement(read)
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__reverse_shorter():
    # tests a read that is shorter than the ref, but covers the region (2000,6000).
    # ==|====|==== #
    #  -|----|---  #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 1000, 9000
    expected_read_region = 3000, 7000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    read = util.reverse_complement(read)
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__reverse_half_overlap():
    # tests a read that is shorter than the ref and covers the region only partially.
    # ==|====|==== #
    #     ---|--   #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 3000, 8000
    expected_read_region = 2000, 5000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    read = util.reverse_complement(read)
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__reverse_full_underlap():
    # tests a read that covers only the middle of the region.
    # ==|====|==== #
    #     --       #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 3000, 5000
    expected_read_region = 0, 2000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    read = util.reverse_complement(read)
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__reverse_clipped_start():
    # tests a read that is shorter than the ref, but covers the region (2000,6000).
    # contains left max(region, read).
    # ==|====|==== #
    #  x|x---|---  #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 1000, 9000
    expected_read_region = 3000, 6000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    # replace first 2kb with rnd dna
    read[:2000] = util.generate_sequence(2000, seed=1)
    read = util.reverse_complement(read)
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__reverse_clipped_half_overlap():
    # tests a read that is shorter than the ref and covers the region only partially.
    # ==|====|==== #
    #     x--|--   #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 3000, 8000
    expected_read_region = 2000, 4000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    read[:1000] = util.generate_sequence(1000, seed=1)
    read = util.reverse_complement(read)
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


def test_util_get_interval_on_read_in_region__reverse_clipped_half_overlap_end():
    # tests a read that is shorter than the ref and covers the region only partially.
    # ==|====|==== #
    #  -|--x       #
    reflength = 10000
    region = 2000, 6000
    read_subsequence = 1000, 5000
    expected_read_region = 1000, 3000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[read_subsequence[0] : read_subsequence[1]]
    read[3000:] = util.generate_sequence(1000, seed=1)
    read = util.reverse_complement(read)
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        interval = util.get_interval_on_read_in_region(
            a=alignment, start=region[0], end=region[1]
        )
        # assert if interval is within 10 bp of expected_read_region
        assert (
            expected_read_region[0] - 10 <= interval[0] <= expected_read_region[0] + 10
        )
        assert (
            expected_read_region[1] - 10 <= interval[1] <= expected_read_region[1] + 10
        )


# %% =============================================================================
# test get_ref_pitx_on_read
def test_get_ref_pitx_on_read_simple():
    # ==|====|==== #
    #   |----|     #
    reflength = 10000
    region = 2000, 6000
    reference = util.generate_sequence(reflength, seed=0)
    read = reference.copy()[region[0] : region[1]]
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get interval on read
        positions_to_test = [
            1998,
            1999,
            2000,
            2001,
            5999,
            6000,
            6001,
            6002,
        ]  # positions on reference
        expected = [0, 0, 0, 1, 3999, 4000, 4000, 4000]  # positions on read
        for i, pos in enumerate(positions_to_test):
            ref_pitx = util.get_ref_pitx_on_read(
                alignment=alignment, position=pos, direction=util.Direction.NONE
            )
            assert ref_pitx[0] == expected[i]


def test_get_read_pitx_on_ref_simple():
    # ==|====|==== #
    #  #|----|#    #
    reflength = 10000
    region = 2000, 6000
    reference = util.generate_sequence(reflength, seed=4)
    read = reference.copy()[region[0] : region[1]]
    # add 1kb of random sequence to the start and another 1kb of random sequence to the end of the read
    read = (
        util.generate_sequence(1000, seed=1)
        + read
        + util.generate_sequence(1000, seed=2)
    )
    # make temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        path_reference = Path(tmpdirname) / "reference.fasta"
        util.write_sequences_to_fasta([reference], path_reference, True)
        util.index_reference(path_reference)
        # write read to fasta
        path_read = Path(tmpdirname) / "read.fasta"
        util.write_sequences_to_fasta([read], path_read)
        # align read to reference
        path_bam = Path(tmpdirname) / "alignment.bam"
        util.align_reads_with_minimap(path_reference, path_read, path_bam, "map-ont")
        # load alignment
        alignment = list(pysam.AlignmentFile(path_bam, "rb"))[0]
        # get positions on read
        positions_to_test = [
            1,
            999,
            1000,
            1001,
            4999,
            5000,
            5001,
            5999,
        ]  # positions on read
        expected = [
            2000,
            2000,
            2000,
            2001,
            5999,
            6000,
            6000,
            6000,
        ]  # positions on reference
        for i, pos in enumerate(positions_to_test):
            read_pitx = util.get_read_pitx_on_ref(
                alignment=alignment, position=pos, direction=util.Direction.NONE
            )
            assert read_pitx[0] == expected[i]
