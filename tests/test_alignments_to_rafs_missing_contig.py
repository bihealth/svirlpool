"""Signal extraction (alignments_to_rafs) with reference contigs that the
alignments do not have.

The workflow writes one region per contig of the reference .fai, so a reference
with more contigs than the alignment file's header (decoys, alts, a reference
the reads were not aligned to in full) gives regions on unknown contigs.
pysam's fetch raises a ValueError for these."""

import csv
import gzip
import json
import logging
import random

import pysam
import pytest

from svirlpool.signalprocessing import alignments_to_rafs
from svirlpool.util.datatypes import ReadAlignmentFragment

READ_START, READ_END = 1_000, 3_500
DEL_START, DEL_SIZE = 2_000, 100
N_READS = 5


def _write_reference(tmp_path, ref: dict[str, str]):
    fa = tmp_path / "ref.fa"
    with open(fa, "w") as f:
        for name, seq in ref.items():
            f.write(f">{name}\n{seq}\n")
    pysam.faidx(str(fa))
    return fa


def _regions_from_fai(fai, bed):
    """as create_regions_from_fai in workflows/main.smk"""
    with open(fai) as f, open(bed, "w") as out:
        for line in f:
            chrom, length = line.split("\t")[:2]
            out.write(f"{chrom}\t0\t{length}\n")
    return bed


def _write_bam(tmp_path, ref: dict[str, str], header_contigs: list[str]):
    """N_READS reads on chr1, each with the same DEL_SIZE bp deletion."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": len(ref[name])} for name in header_contigs],
    }
    chr1 = ref["chr1"]
    bam = tmp_path / "alignments.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as f:
        for i in range(N_READS):
            start = READ_START + 10 * i
            seq = chr1[start:DEL_START] + chr1[DEL_START + DEL_SIZE : READ_END]
            a = pysam.AlignedSegment(f.header)
            a.query_name = f"read-{i}"
            a.flag = 0
            a.reference_id = f.get_tid("chr1")
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = [
                (0, DEL_START - start),
                (2, DEL_SIZE),
                (0, READ_END - DEL_START - DEL_SIZE),
            ]
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            f.write(a)
    pysam.index(str(bam))
    return bam


def _read_rafs(path) -> list[ReadAlignmentFragment]:
    with gzip.open(path, "rt") as f:
        return [
            ReadAlignmentFragment.from_unstructured(json.loads(row[3]))
            for row in csv.reader(f, delimiter="\t")
        ]


# the non-separated read filter (cache size > 0) fetches before process_region does
@pytest.mark.parametrize("cache_size", [1000, 0])
def test_contigs_missing_from_the_alignments_are_skipped(tmp_path, caplog, cache_size):
    rng = random.Random(1)
    ref = {
        name: "".join(rng.choice("ACGT") for _ in range(size))
        for name, size in [
            ("chr1", 5_000),
            ("chrNoReads", 3_000),  # in the header, but no reads
            ("chrExtra1", 2_000),  # not in the header
            ("chrExtra2", 2_000),  # not in the header
        ]
    }
    fa = _write_reference(tmp_path, ref)
    bed = _regions_from_fai(f"{fa}.fai", tmp_path / "regions.bed")
    bam = _write_bam(tmp_path, ref, header_contigs=["chr1", "chrNoReads"])
    out = tmp_path / "rafs.tsv.gz"
    args = alignments_to_rafs.get_parser().parse_args([
        "-a", str(bam),
        "-s", "sample",
        "-r", str(bed),
        "-o", str(out),
        "-t", "1",
        "--filter-nonseparated__cache-size-alignments-filtering", str(cache_size),
    ])  # fmt: skip

    with caplog.at_level(logging.INFO, logger=alignments_to_rafs.__name__):
        alignments_to_rafs.run(args)

    rafs = _read_rafs(out)
    assert sorted(raf.read_name for raf in rafs) == [
        f"read-{i}" for i in range(N_READS)
    ]
    assert {raf.reference_name for raf in rafs} == {"chr1"}
    for raf in rafs:
        assert [
            (s.sv_type, s.ref_start, s.ref_end, s.size) for s in raf.SV_signals
        ] == [(1, DEL_START, DEL_START + DEL_SIZE, DEL_SIZE)]

    # one message for all skipped contigs, none for the contig without reads
    skipped = [r.getMessage() for r in caplog.records if "chrExtra1" in r.getMessage()]
    assert len(skipped) == 1
    assert "chrExtra2" in skipped[0]
    assert "chrNoReads" not in skipped[0]
