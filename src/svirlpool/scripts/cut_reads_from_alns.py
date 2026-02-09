# One of the modules of svirlpool to be used in the command line interface

import argparse
import gzip
import logging
from pathlib import Path

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval, IntervalTree

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
) -> dict[str, list[pysam.AlignedSegment]]:
    """Returns a dict of the form readname:[pysam.alignedSegment]"""
    result = {}
    region_interval = Interval(region[1], region[2])
    with pysam.AlignmentFile(str(alignments), "rb") as f:
        for aln in f.fetch(*region):
            if no_secondary and aln.is_secondary:
                log.info(f"skipping secondary alignment {aln.query_name}")
                continue
            aln_interval = Interval(aln.reference_start, aln.reference_end)
            if aln_interval.overlaps(region_interval):
                if aln.query_name not in result:
                    result[aln.query_name] = []
                result[aln.query_name].append(aln)
    return result


def get_read_alignment_intervals_in_region(
    region: tuple[str, int, int],
    dict_alignments: dict[str, list[pysam.AlignedSegment]],
    buffer_clipped_length: int,
) -> dict[str, tuple[int, int, int, int]]:
    # check if all elements in alignments are of type pysam.AlignedSegment
    for alignments in dict_alignments.values():
        for aln in alignments:
            assert isinstance(aln, pysam.AlignedSegment), (
                f"aln {aln} is not a pysam.AlignedSegment"
            )
    result = {}
    for alignments in dict_alignments.values():
        for aln in alignments:
            if aln.query_name not in result:
                result[aln.query_name] = []
            read_start, read_end = util.get_interval_on_read_in_region(
                a=aln,
                start=region[1],
                end=region[2],
                buffer_clipped_length=buffer_clipped_length,
            )
            ref_start, ref_end = util.get_interval_on_ref_in_region(
                a=aln, start=read_start, end=read_end
            )
            # check read_start and read_end if they are numerical and end-start is greater than 0
            if (
                read_start is not None
                and read_end is not None
                and read_end > read_start
            ):
                result[aln.query_name].append((
                    read_start,
                    read_end,
                    ref_start,
                    ref_end,
                ))
    return result


def get_max_extents_of_read_alignments_in_region(
    dict_all_intervals: dict[str, list[tuple[int, int, int, int]]],
) -> dict[str, tuple[int, int, int, int]]:
    """Returns a dict of the form {readname:(read_start,read_end,ref_start,ref_end)}"""
    dict_max_extents = {}
    for readname in dict_all_intervals.keys():
        intervals = dict_all_intervals[readname]
        min_start = min(
            (start, ref_start) for (start, end, ref_start, ref_end) in intervals
        )
        max_end = max((end, ref_end) for (start, end, ref_start, ref_end) in intervals)
        dict_max_extents[readname] = (
            min_start[0],
            max_end[0],
            min_start[1],
            max_end[1],
        )
    return dict_max_extents


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


def get_full_read_sequences_of_alignments(
    dict_alignments: dict[str, list[pysam.AlignedSegment]], alignments: Path
) -> dict[str, SeqRecord]:
    """Retrieves all DNA sequences of all given read alignments. The read DNA is in original orientation, as it was given in the original fasta file."""
    # iterate all alignments of a sampleID across all crIDs
    dict_supplementary_positions: dict[str, list[tuple[str, int]]] = {}
    for alns in dict_alignments.values():
        for aln in alns:
            if not aln.query_sequence:
                continue
            if (
                aln.infer_read_length() != len(aln.query_sequence)
                and aln.query_name not in dict_supplementary_positions
            ):
                dict_supplementary_positions[aln.query_name] = (
                    get_all_supplementary_positions_with_full_read_sequence(aln=aln)
                )

    n_reads_with_full_sequences = len(set(dict_supplementary_positions.keys()))
    n_reads_without_full_sequences = (
        len(set(dict_alignments.keys())) - n_reads_with_full_sequences
    )
    log.info(
        f"{n_reads_with_full_sequences} reads have full read sequences, {n_reads_without_full_sequences} reads have no full read sequences"
    )

    dict_positions = {
        chrom: []
        for chrom in {
            chr for _l in dict_supplementary_positions.values() for chr, pos in _l
        }
    }
    for _readname, alns in dict_supplementary_positions.items():
        for chr, pos in alns:
            dict_positions[chr].append(pos)
    # then sort the list of start positions in each chromosome
    for chrom in dict_positions.keys():
        dict_positions[chrom] = sorted(dict_positions[chrom])
    # then merge any that are within 50kb of the next one
    # to do so, build an intervaltree, where 50_000 is subtracted from start and added to end
    for chrom in dict_positions.keys():
        it = IntervalTree()
        for start in dict_positions[chrom]:
            it[max(0, start - 50_000) : start + 50_000] = 1
        # then merge the intervaltree
        it.merge_overlaps(data_reducer=lambda a, b: a + b)
        # then extract the merged intervals and fill dict_positions with (start,end) tuples
        dict_positions[chrom] = [
            (interval.begin + 50_000 - 1, interval.end - 50_000 + 1) for interval in it
        ]
    # the regions are then used to fetch the alignments from the bam file
    regions = [
        (chrom, min(start, end), max(start, end))
        for chrom in dict_positions.keys()
        for start, end in dict_positions[chrom]
    ]
    log.info(
        "To find the full read sequences, the following regions are fetched from the alignments file:"
    )
    dict_read_sequences = {}
    with pysam.AlignmentFile(alignments, "rb") as f:
        for region in regions:
            log.info(f"fetching from {region}")
            for aln in f.fetch(region[0], region[1], region[2]):
                if aln.query_name in dict_supplementary_positions.keys():
                    if not aln.query_sequence:
                        continue
                    if aln.query_name not in dict_read_sequences.keys():
                        dict_read_sequences[aln.query_name] = None
                    if not aln.infer_read_length() == len(aln.query_sequence):
                        # log.info(f"read seq for {aln.query_name} is not full. continue.")
                        continue
                    else:
                        seq = (
                            Seq(aln.query_sequence).reverse_complement()
                            if aln.is_reverse
                            else Seq(aln.query_sequence)
                        )
                        qualities = None
                        if aln.query_qualities:
                            qualities = {
                                "phred_quality": aln.query_qualities[::-1]
                                if aln.is_reverse
                                else aln.query_qualities
                            }
                        else:
                            log.warning(
                                f"read {aln.query_name} has no quality scores. This is unexpected for a full read sequence. Please check your alignments file."
                            )
                        dict_read_sequences[aln.query_name] = SeqRecord(
                            seq=seq,
                            letter_annotations=qualities,
                            id=aln.query_name,
                            name=aln.query_name,
                        )
                        # log.info(f"found a full read sequence for {aln.query_name}")
            # check if all keys of dict_read_sequences have a SeqRecord
            # if so, break this loop
            if all(v is not None for v in dict_read_sequences.values()):
                break
    # for all other alignments, just add the sequence and reverse flag
    for alignments in dict_alignments.values():
        for aln in alignments:
            if not aln.query_sequence:
                continue
            if (
                aln.query_name not in dict_read_sequences.keys()
                and aln.infer_read_length() == len(aln.query_sequence)
            ):
                seq = (
                    Seq(aln.query_sequence).reverse_complement()
                    if aln.is_reverse
                    else Seq(aln.query_sequence)
                )
                qualities = (
                    {"phred_quality": aln.query_qualities}
                    if aln.query_qualities
                    else None
                )
                dict_read_sequences[aln.query_name] = SeqRecord(
                    seq=seq,
                    letter_annotations=qualities,
                    id=aln.query_name,
                    name=aln.query_name,
                )
    # check if all elements in dict_read_sequences are of type SeqRecord
    for readname, record in dict_read_sequences.items():
        assert isinstance(record, SeqRecord), (
            f"record {record} is not a SeqRecord. dict_read_sequences without record are: {' '.join({readname for readname, record in dict_read_sequences.items() if record is None})}"
        )
    # filter dict_read_sequences for reads that have a full sequence
    dict_read_sequences = {
        readname: record
        for readname, record in dict_read_sequences.items()
        if record is not None
    }
    # report all reads for which no full reads could be found
    reads_without_full_sequences = set(dict_alignments.keys()) - set(
        dict_read_sequences.keys()
    )
    log.warning(
        f"the following reads have no full read sequences: {reads_without_full_sequences}"
    )
    # # check if for each readname in dict_alignments if there is a SeqRecord in dict_read_sequences
    # for readname in dict_alignments.keys():
    #     if readname not in dict_read_sequences.keys():
    #         raise ValueError(f"readname {readname} in dict_alignments has no SeqRecord in dict_read_sequences")
    return dict_read_sequences


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
        region=region,
        dict_alignments=dict_alignments,
        buffer_clipped_length=buffer_clipped_length,
    )
    dict_max_intervals = get_max_extents_of_read_alignments_in_region(dict_intervals)
    # log.info(f"found intervals:\n{'\n'.join([f'{k}:{v}' for k,v in dict_max_intervals.items()])}")
    # get all full read sequences of alignments
    # log.info(f"fetching full read sequences of alignments...")
    read_records = get_full_read_sequences_of_alignments(
        dict_alignments=dict_alignments, alignments=alignments
    )
    # then cut the reads
    # log.info(f"cutting the read sequences to their according sizes...")
    return cut_reads(
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
    # parsed region can include , separating thousands
    parsed_region = parse_region(args.region)
    dict_cut_reads = cut_reads_from_alignments(
        alignments=args.input,
        region=parsed_region,
        buffer_clipped_length=args.buffer_clipped_length,
    )
    write_cut_reads_to_output_file(dict_cut_reads=dict_cut_reads, output=args.output)


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
