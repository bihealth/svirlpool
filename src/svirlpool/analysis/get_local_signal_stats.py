# this script plots a summed indel histogram, the sorted signal types and sizes
# and for each line the sorted readnames for a given region and a given bam file
# to the console using ascii characters

import argparse
from pathlib import Path, PosixPath

import pysam
from intervaltree import Interval

from ..scripts import alignments_to_rafs, datatypes


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, start_end = region.split(":")
    start, end = start_end.split("-")
    start = start.replace(",", "")
    end = end.replace(",", "")
    return (chrom, int(start), int(end))


def parse_alignments(
    path_alignments: PosixPath, region: tuple[str, int, int]
) -> list[pysam.AlignedSegment]:
    bamfile = pysam.AlignmentFile(path_alignments, "rb")
    return list(bamfile.fetch(*region))


# parse breakends to tuples start, end, size (ref coordinates)
def parse_signals(
    alignments: list[pysam.AlignedSegment],
    min_insertion_size: int,
    min_deletion_size: int,
    min_bnd_size: int,
    region: tuple[str, int, int],
) -> dict[str, dict[str, int]]:
    ri = Interval(region[1], region[2])
    rafs: list[datatypes.ReadAlignmentFragment] = [
        alignments_to_rafs.parse_ReadAlignmentFragment_from_alignment(
            alignment=aln,
            samplename="sample",
            min_bnd_size=min_bnd_size,
            min_signal_size=min(min_insertion_size, min_deletion_size),
            filter_density_min_bp=min(min_insertion_size, min_deletion_size),
            filter_density_radius=100,
        )
        for aln in alignments
    ]
    results = {
        raf.read_name: {"insertions": 0, "deletions": 0, "breakends": 0} for raf in rafs
    }
    for raf in rafs:
        for signal in raf.SV_signals:
            if signal.sv_type == 0 and ri.overlaps(
                Interval(signal.ref_start, signal.ref_end)
            ):
                results[raf.read_name]["insertions"] += signal.size
            elif signal.sv_type == 1 and ri.overlaps(
                Interval(signal.ref_start, signal.ref_end)
            ):
                results[raf.read_name]["deletions"] += signal.size
            elif signal.sv_type >= 4 and ri.overlaps(
                Interval(signal.ref_start, signal.ref_end)
            ):
                results[raf.read_name]["breakends"] += signal.size
    return results


def print_results(
    counts: dict[str, dict[str, int]], figpath: Path, region: tuple[str, int, int]
):
    indel_sums = {
        readname: counts[readname]["insertions"] - counts[readname]["deletions"]
        for readname in counts
    }
    sorted_readnames = sorted(indel_sums.keys(), key=lambda x: indel_sums[x])
    # create a plot. iterate each readname and print the indel_sum of it
    from matplotlib import pyplot as plt

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)
    ax.plot([indel_sums[readname] for readname in sorted_readnames])
    ax.set_title(f"Indel sums for region {region[0]}:{region[1]}-{region[2]}")
    # label x axis "reads"
    ax.set_xlabel("reads")
    ax.set_ylabel("sum( insertions - deletions )")
    plt.savefig(figpath)


def run(
    alignments: Path,
    region_str: str,
    min_insertion_size: int,
    min_deletion_size: int,
    min_breakend_size: int,
    output: Path,
):
    region = parse_region(region=region_str)
    alignments = parse_alignments(path_alignments=alignments, region=region)
    signals: dict[str, dict[str, int]] = parse_signals(
        alignments=alignments,
        min_insertion_size=min_insertion_size,
        min_deletion_size=min_deletion_size,
        min_bnd_size=min_breakend_size,
        region=region,
    )
    print_results(counts=signals, figpath=output, region=region)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Get local signal stats and print them to a figure"
    )
    parser.add_argument(
        "-a",
        "--alignments",
        type=Path,
        required=True,
        help="Path to the alignments (bam) file",
    )
    parser.add_argument(
        "-r",
        "--region",
        type=str,
        required=True,
        help="Region to get the signals from, format: chr:start-end",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output figure in .png format",
    )
    parser.add_argument(
        "-mi",
        "--min-insertion-size",
        type=int,
        default=6,
        help="Minimum size of an insertion to be considered",
    )
    parser.add_argument(
        "-md",
        "--min-deletion-size",
        type=int,
        default=6,
        help="Minimum size of a deletion to be considered",
    )
    parser.add_argument(
        "-mb",
        "--min-breakend-size",
        type=int,
        default=150,
        help="Minimum size of a breakend to be considered",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(
        alignments=args.alignments,
        region_str=args.region,
        min_insertion_size=args.min_insertion_size,
        min_deletion_size=args.min_deletion_size,
        min_breakend_size=args.min_breakend_size,
        output=args.output,
    )
    return


if __name__ == "__main__":
    main()
