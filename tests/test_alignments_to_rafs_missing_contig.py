"""Signal extraction (alignments_to_rafs) when the contigs of the reference and
of the alignments differ.

The workflow writes one region per contig of the reference .fai.
- Reference contigs that the alignment file does not have (decoys, alts, a
  reference the reads were not aligned to in full) are skipped silently.
  pysam's fetch raises a ValueError for them.
- Reads on contigs that the reference does not have would be lost without
  notice, so they are an error. Contigs of the alignment header without reads
  may be missing from the reference; they are skipped silently too."""

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


def _sequences(sizes: dict[str, int]) -> dict[str, str]:
    rng = random.Random(1)
    return {
        name: "".join(rng.choice("ACGT") for _ in range(size))
        for name, size in sizes.items()
    }


def _write_reference(path, seqs: dict[str, str], contigs: list[str]):
    with open(path, "w") as f:
        for name in contigs:
            f.write(f">{name}\n{seqs[name]}\n")
    pysam.faidx(str(path))
    return path


def _regions_from_fai(fai, bed):
    """as create_regions_from_fai in workflows/main.smk"""
    with open(fai) as f, open(bed, "w") as out:
        for line in f:
            chrom, length = line.split("\t")[:2]
            out.write(f"{chrom}\t0\t{length}\n")
    return bed


def _write_alignments(
    path,
    seqs: dict[str, str],
    header_contigs: list[str],
    read_contigs: list[str],
    cram_reference=None,
):
    """N_READS reads on each of read_contigs, each with the same DEL_SIZE bp
    deletion. A CRAM if cram_reference is given, else a BAM."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": len(seqs[name])} for name in header_contigs],
    }
    mode, kwargs = "wb", {}
    if cram_reference:
        mode, kwargs = "wc", {"reference_filename": str(cram_reference)}
    with pysam.AlignmentFile(str(path), mode, header=header, **kwargs) as f:
        for chrom in header_contigs:  # in coordinate order
            if chrom not in read_contigs:
                continue
            for i in range(N_READS):
                start = READ_START + 10 * i
                seq = (
                    seqs[chrom][start:DEL_START]
                    + seqs[chrom][DEL_START + DEL_SIZE : READ_END]
                )
                a = pysam.AlignedSegment(f.header)
                a.query_name = f"{chrom}-read-{i}"
                a.flag = 0
                a.reference_id = f.get_tid(chrom)
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
    pysam.index(str(path))
    return path


def _run(tmp_path, reference, alignments, cache_size: int = 1000):
    """alignments_to_rafs as the workflow runs it, on one region per contig of
    the reference. Returns the written ReadAlignmentFragments."""
    bed = _regions_from_fai(f"{reference}.fai", tmp_path / "regions.bed")
    out = tmp_path / "rafs.tsv.gz"
    args = alignments_to_rafs.get_parser().parse_args([
        "-a", str(alignments),
        "-s", "sample",
        "-r", str(bed),
        "--reference", str(reference),
        "-o", str(out),
        "-t", "1",
        "--filter-nonseparated__cache-size-alignments-filtering", str(cache_size),
    ])  # fmt: skip
    alignments_to_rafs.run(args)
    with gzip.open(out, "rt") as f:
        return [
            ReadAlignmentFragment.from_unstructured(json.loads(row[3]))
            for row in csv.reader(f, delimiter="\t")
        ]


def _assert_chr1_signals(rafs: list[ReadAlignmentFragment]) -> None:
    assert sorted(raf.read_name for raf in rafs) == [
        f"chr1-read-{i}" for i in range(N_READS)
    ]
    assert {raf.reference_name for raf in rafs} == {"chr1"}
    for raf in rafs:
        assert [
            (s.sv_type, s.ref_start, s.ref_end, s.size) for s in raf.SV_signals
        ] == [(1, DEL_START, DEL_START + DEL_SIZE, DEL_SIZE)]


# --- reference contigs that the alignments do not have --- #


# the non-separated read filter (cache size > 0) fetches before process_region does
@pytest.mark.parametrize("cache_size", [1000, 0])
def test_contigs_missing_from_the_alignments_are_skipped(tmp_path, caplog, cache_size):
    seqs = _sequences({
        "chr1": 5_000,
        "chrNoReads": 3_000,  # in the header, but no reads
        "chrExtra1": 2_000,  # not in the header
        "chrExtra2": 2_000,  # not in the header
    })
    reference = _write_reference(tmp_path / "ref.fa", seqs, list(seqs))
    bam = _write_alignments(
        tmp_path / "alignments.bam",
        seqs,
        header_contigs=["chr1", "chrNoReads"],
        read_contigs=["chr1"],
    )

    with caplog.at_level(logging.DEBUG, logger=alignments_to_rafs.__name__):
        rafs = _run(tmp_path, reference, bam, cache_size=cache_size)

    _assert_chr1_signals(rafs)
    # the skipped contigs are not reported at the default log levels
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert not [
        r
        for r in caplog.records
        if r.levelno >= logging.INFO and "chrExtra" in r.getMessage()
    ]


# --- alignment contigs that the reference does not have --- #


def test_reads_on_a_contig_missing_from_the_reference_are_an_error(tmp_path):
    seqs = _sequences({"chr1": 5_000, "chrUnknown": 5_000})
    reference = _write_reference(tmp_path / "ref.fa", seqs, ["chr1"])
    bam = _write_alignments(
        tmp_path / "alignments.bam",
        seqs,
        header_contigs=["chr1", "chrUnknown"],
        read_contigs=["chr1", "chrUnknown"],
    )

    with pytest.raises(ValueError, match="chrUnknown") as error:
        _run(tmp_path, reference, bam)

    assert str(bam) in str(error.value)
    assert str(reference) in str(error.value)
    # checked before any region is processed
    assert not (tmp_path / "rafs.tsv.gz").exists()


def test_the_error_lists_the_first_ten_missing_contigs(tmp_path):
    unknown = [f"chrUn{i:02d}" for i in range(12)]
    seqs = _sequences({"chr1": 5_000} | dict.fromkeys(unknown, 4_000))
    reference = _write_reference(tmp_path / "ref.fa", seqs, ["chr1"])
    bam = _write_alignments(
        tmp_path / "alignments.bam",
        seqs,
        header_contigs=list(seqs),
        read_contigs=list(seqs),
    )

    with pytest.raises(ValueError) as error:
        _run(tmp_path, reference, bam)

    message = str(error.value)
    assert all(name in message for name in unknown[:10])
    assert not any(name in message for name in unknown[10:])
    assert "2 more" in message


def test_header_contigs_without_reads_may_be_missing_from_the_reference(
    tmp_path, caplog
):
    seqs = _sequences({"chr1": 5_000, "chrDecoy": 3_000})
    reference = _write_reference(tmp_path / "ref.fa", seqs, ["chr1"])
    bam = _write_alignments(
        tmp_path / "alignments.bam",
        seqs,
        header_contigs=["chr1", "chrDecoy"],
        read_contigs=["chr1"],
    )

    with caplog.at_level(logging.DEBUG, logger=alignments_to_rafs.__name__):
        rafs = _run(tmp_path, reference, bam)  # no error

    _assert_chr1_signals(rafs)
    # and not reported at the default log levels
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert not [
        r
        for r in caplog.records
        if r.levelno >= logging.INFO and "chrDecoy" in r.getMessage()
    ]


def test_a_cram_index_has_no_read_counts_so_all_header_contigs_must_be_in_the_reference(
    tmp_path,
):
    # pysam reports 0 mapped reads for every contig of a CRAM index
    seqs = _sequences({"chr1": 5_000, "chrDecoy": 3_000})
    cram_reference = _write_reference(tmp_path / "cram_ref.fa", seqs, list(seqs))
    reference = _write_reference(tmp_path / "ref.fa", seqs, ["chr1"])
    cram = _write_alignments(
        tmp_path / "alignments.cram",
        seqs,
        header_contigs=["chr1", "chrDecoy"],
        read_contigs=["chr1"],
        cram_reference=cram_reference,
    )

    with pytest.raises(ValueError, match="chrDecoy"):
        _run(tmp_path, reference, cram)
