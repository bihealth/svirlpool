# One of the modules of svirlpool to be used in the command line interface

import argparse
import gzip
import logging
from pathlib import Path

import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval

from ..localassembly.consensus import (
    get_full_read_sequences_of_alignments,
    get_max_extents_of_read_alignments_on_cr,
    get_read_alignment_intervals_in_region,
    trim_reads,
)
from ..util import util

log = logging.getLogger(__name__)


def cut_reads(
    dict_alignments: dict[str, list[pysam.AlignedSegment]],
    intervals: dict[str, tuple[int, int, int, int]],
    read_records: dict[str, SeqRecord],
) -> dict[str, SeqRecord]:
    """Returns a dict of the form {readname:SeqRecord} with all reads cut to the provided intervals"""
    # dict_alignments is of the form: dict_alignments[cr.crID][sampleID]=[pysam.AlignedSegment]
    # returns a list of read sequence records
    # check if dict_alignments is of the correct form
    for _readname, alignments in dict_alignments.items():
        for aln in alignments:
            assert isinstance(aln, pysam.AlignedSegment), (
                f"aln {aln} is not a pysam.AlignedSegment"
            )
    # check if intervals is of the correct form
    for readname in intervals.keys():
        assert isinstance(readname, str), f"readname {readname} is not a string"
        assert isinstance(intervals[readname], tuple), (
            f"intervals[readname] {intervals[readname]} is not a tuple"
        )
        assert len(intervals[readname]) == 4, (
            f"intervals[readname] {intervals[readname]} is not of length 4"
        )
        for i in range(4):
            assert isinstance(intervals[readname][i], int), (
                f"intervals[readname][{i}] {intervals[readname][i]} is not an integer"
            )

    # log.info(f"cutting reads from alignments")
    cut_reads: dict[str, SeqRecord] = {}
    for alignments in dict_alignments.values():
        for aln in alignments:
            if aln.query_name in intervals:
                if aln.query_name not in read_records:
                    log.warning(
                        f"read {aln.query_name} not found in read records. Skipping."
                    )
                    continue
                if aln.query_name in cut_reads:
                    continue
                if aln.infer_read_length() != len(read_records[aln.query_name]):
                    raise ValueError(
                        f"read length of {aln.query_name} does not match the length of the read record. read length from alignment: {aln.infer_read_length()}, read length from read record: {len(read_records[aln.query_name])}"
                    )
                start, end, ref_start, ref_end = intervals[aln.query_name]
                record = read_records[aln.query_name][start:end]
                record.name = aln.query_name
                record.id = aln.query_name
                # log.info(f"{aln.query_name} cut from {start} to {end} (length={end-start}) on ref: {ref_start} to {ref_end}")
                record.description = f"start={start}, end={end}, ref_start={min(ref_start, ref_end)}, ref_end={max(ref_start, ref_end)}"
                cut_reads[aln.query_name] = record
    return cut_reads


def get_read_alignments_from_region(
    region: tuple[str, int, int], alignments: Path, no_secondary: bool = True
) -> dict[int, list[pysam.AlignedSegment]]:
    """Returns a dict of the form readname:[pysam.alignedSegment]"""
    result = {0: []}
    region_interval = Interval(region[1], region[2])
    with pysam.AlignmentFile(str(alignments), "rb") as f:
        for aln in f.fetch(*region):
            if no_secondary and aln.is_secondary:
                log.info(f"skipping secondary alignment {aln.query_name}")
                continue
            aln_interval = Interval(aln.reference_start, aln.reference_end)
            if aln_interval.overlaps(region_interval):
                result[0].append(aln)
    return result


def get_all_supplementary_positions_with_full_read_sequence(
    aln: pysam.AlignedSegment,
) -> list[tuple[str, int]]:
    """
    get the next supplementary alignment with a full read sequence. Only apply this to alignments with hard clippings!
    Returns None if no such alignment is found.
    """
    if not aln.has_tag("SA"):
        raise ValueError(
            f"Alignment {aln.query_name} in {aln.reference_name}:{aln.reference_start}-{aln.reference_end} \
has no supplementary alignments. This should not happen. For each aligned read should exist at least one alignment \
that brings the whole read sequence with it - identified by absence of 'H' or hard clippings in the cigar string. \
This is only tested with minimap2 alignments so far. Consider reporting to the author if you have used a different \
read aligner so we can fix it."
        )
    result = []
    lines = aln.get_tag("SA").split(";")
    for line in lines:
        fields = line.split(",")
        if len(fields) < 6:
            continue
        # check if H not in cigar string, if so, return chrom,start
        result.append((fields[0], int(fields[1])))
    return result


def write_cut_reads_to_output_file(
    dict_cut_reads: dict[str, SeqRecord], output: Path
) -> None:
    extension = Path(output).suffix
    is_fastq = False
    if any(ext in extension for ext in ["fastq", "fq"]):
        if all(r.letter_annotations for r in dict_cut_reads.values()):
            is_fastq = True
        else:
            raise ValueError(
                f"output file {output} has a fastq extension but no quality scores are present. Please provide a fasta extension for these alignments."
            )
    # check if extension is .gz and if so, open the file with gzip
    if Path(output).suffix == ".gz":
        with gzip.open(output, "wt") as f:
            if is_fastq:
                SeqIO.write(dict_cut_reads.values(), f, "fastq")
            else:
                SeqIO.write(dict_cut_reads.values(), f, "fasta")
    else:
        with open(output, "w") as f:
            if is_fastq:
                SeqIO.write(dict_cut_reads.values(), f, "fastq")
            else:
                SeqIO.write(dict_cut_reads.values(), f, "fasta")


def cut_reads_from_alignments(
    alignments: Path, region: tuple[str, int, int], buffer_clipped_length: int
) -> dict[str, SeqRecord]:
    # fetch all read alignments in region
    dict_alignments = get_read_alignments_from_region(
        region=region, alignments=alignments, no_secondary=True
    )
    # log.info(f"found {len(set(dict_alignments.keys()))} reads with {len(dict_alignments)} read alignments in region {region}")
    # print(sorted(dict_alignments.keys()))
    # filter for non-separated reads
    # get all read alignment intervals in region
    dict_intervals = get_read_alignment_intervals_in_region(
        region_start=region[1],
        regions_end=region[2],
        alignments=[aln for alns in dict_alignments.values() for aln in alns],
        buffer_clipped_length=buffer_clipped_length,
    )
    log.debug(f"found intervals: {dict_intervals}")
    dict_max_intervals = get_max_extents_of_read_alignments_on_cr(dict_intervals)
    log.debug(f"max intervals: {dict_max_intervals}")
    # log.info(f"found intervals:\n{'\n'.join([f'{k}:{v}' for k,v in dict_max_intervals.items()])}")
    # get all full read sequences of alignments
    # log.info(f"fetching full read sequences of alignments...")
    read_records = get_full_read_sequences_of_alignments(
        dict_alignments=dict_alignments, path_alignments=alignments
    )
    # then cut the reads
    # log.info(f"cutting the read sequences to their according sizes...")
    return trim_reads(
        dict_alignments=dict_alignments,
        intervals=dict_max_intervals,
        read_records=read_records,
    )


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, start_end = region.split(":")
    start, end = start_end.split("-")
    start = start.replace(",", "")
    end = end.replace(",", "")
    return (chrom, int(start), int(end))


def run(args, **kwargs):
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        force=True,
    )
    # parsed region can include , separating thousands
    parsed_region = parse_region(args.region)
    dict_cut_reads = cut_reads_from_alignments(
        alignments=args.input,
        region=parsed_region,
        buffer_clipped_length=args.buffer_clipped_length,
    )
    write_cut_reads_to_output_file(dict_cut_reads=dict_cut_reads, output=args.output)
    if args.reference is not None:
        bam_output = args.output.parent / (args.output.name.split(".")[0] + ".bam")
        util.align_reads_with_minimap(
            reference=args.reference,
            reads=args.output,
            bamout=bam_output,
        )


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the alignments file (bam/sam).",
    )
    parser.add_argument(
        "-r",
        "--region",
        type=str,
        required=True,
        help="Region to cut reads from. Format: chr:start-end",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output fasta/fastq (.gz) file. If the filename extension is .gz, the output file will be written as a gzip compressed file.",
    )
    parser.add_argument(
        "--buffer-clipped-length",
        type=int,
        default=1000,
        help="The maximum number of bases that are included in the cut sequences, if they have been hard or soft clipped within the region bounds. Defaults to 1000.",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=None,
        help="Optional path to a reference genome. If provided, the output reads will be aligned to the reference using minimap2 and written as a BAM file with the same base name as the output file.",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="WARNING",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level. Defaults to WARNING.",
    )


def get_parser():
    parser = argparse.ArgumentParser(description="")
    add_arguments(parser)
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
