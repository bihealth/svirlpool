# filter non-separated reads
# N number of maximum aligned segments in cache
# P path to bam
# R region chr,start,end to fetch from bam file and iterate over
# S set of read names that are filtered out
# %%
import argparse
import multiprocessing as mp
import subprocess
import tempfile
import typing
from math import ceil
from pathlib import Path
from shlex import split

import pysam
from logzero import logger as log

# %%


def find_candidates(
    region: typing.Tuple[str, int, int],
    bam_path: Path,
    cache_size: int,
    min_overlap: float,
    min_fragment_size: int,
    max_inversion_coverage: float,
) -> set:
    region_size = region[2] - region[1]
    if region_size < 1:
        return set()
    # create cache that contains elements: [hash(read_name), readname, ref_start, ref_end, reverse]
    # store all alignments that are candidates for filtering
    # after building the cache, identify regions of several overlaps.
    cache = []
    # alignment sizes are stored to compute the average alignment coverage
    # if a significant portion of the coverage is by filter candidates, then don't filter them,
    # because they might be real inversions
    alignment_sizes = []
    candidates = []
    # log.info("finding candidate non-separated reads..")
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(*region):
            alignment_sizes.append(read.reference_end - read.reference_start)
            item = [
                hash(read.query_name),
                read.query_name,
                read.reference_start,
                read.reference_end,
                read.is_reverse,
            ]
            # check if there is another alignment of this read in the cache by comparing the hashes.
            for other in cache:
                # test: hash, readname, direction, then overlap
                if item[0] == other[0] and item[1] == other[1] and item[4] != other[4]:
                    # compute the overlap between the two alignments
                    overlap = min(item[3], other[3]) - max(item[2], other[2])
                    size_item = item[3] - item[2]
                    size_other = other[3] - other[2]
                    overlapping_ratio = overlap / min(size_item, size_other)
                    # it is a candidate if the overlap is at least 50% of the smaller alignment
                    # and the smaller alignment is bigger than 500 bp
                    if (
                        overlapping_ratio >= min_overlap
                        and min(size_item, size_other) >= min_fragment_size
                    ):
                        candidates.append(item)
            # add item to cache
            cache.append(item)
            # pop items until cache is back to cache_size
            while len(cache) > cache_size:
                cache.pop(0)
    if len(candidates) == 0:
        return set()
    # empty cache
    cache = []
    # calc mean depth of coverage
    # find regions that contain potential true inversions
    # write intervals with coverage by candidates
    # sort candidates by start position
    candidates = sorted(candidates, key=lambda x: x[2])
    # create a list of starting and ending positions of all candidates. each position is a tuple (position, +1/-1)
    positions = []
    # log.info("finding regions of coverage by candidates..")
    for candidate in candidates:
        positions.append((candidate[2], 1))
        positions.append((candidate[3], -1))
    # sort positions by position
    positions = sorted(positions, key=lambda x: x[0])
    # then iterate this list. keep track of current coverage. If a candidate starts, add to coverage,
    # if a candidate ends, subtract from coverage.
    # each time, the current position changes, create a new interval with the current coverage
    intervals = []
    # coverage of candidates
    current_candidates_coverage = 0
    last_pos = 0
    for pos, cov_change in positions:
        # add coverage counts as long as current_pos doesn't change
        if last_pos != pos and current_candidates_coverage > 0:
            # if pos changes, then create a new interval with the current coverage
            # if current coverage is 0, then don't create an interval
            intervals.append([last_pos, pos, current_candidates_coverage])
        current_candidates_coverage += cov_change
    # iterate candidates again, and check if any interval they overlap with
    # has a ratio current_candidates_coverage / mean_coverage of less than 10%. If so,
    # then add to filtered_reads
    # log.info("reducing filter to candidates with low coverage..")
    filtered_reads = set()
    mean_coverage = sum(alignment_sizes) / region_size
    max_coverage = ceil(mean_coverage * max_inversion_coverage)
    for _, o_name, o_start, o_end, _ in candidates:
        for start, end, cov in intervals:
            left = max(o_start, start)
            right = min(o_end, end)
            overlap = right - left
            if overlap > 0 and cov < max_coverage:
                filtered_reads.add(o_name)
    bam.close()
    return filtered_reads


def filter_reads(
    region: typing.Tuple[str, int, int], readnames: set, bam_path: Path, out_path: Path
):
    log.info(
        f"filtering reads from {region[0]}:{region[1]}-{region[2]} in file {bam_path} and writing to file {out_path}."
    )
    bam = pysam.AlignmentFile(bam_path, "rb")
    with pysam.AlignmentFile(out_path, "wb", header=bam.header) as out:
        for alignment in bam.fetch(*region):
            if alignment.query_name not in readnames:
                out.write(alignment)
    bam.close()
    return


# to parallelize, split the regons provided in a bed file into N regions
# to determine if the regions are to be split, first compute the number of total bases in the
# regions and then divide the total number by the number of provided threads to get the chunk size.
# if a chromosome is much longer than chunksize, then split it into regions of chunksize.
# otherwise, use the provided regions.
def split_regions(
    bedfile: Path, threads: int
) -> typing.List[typing.Tuple[str, int, int]]:
    # first, compute the total number of bases in the regions
    total_bases = 0
    with open(bedfile, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            total_bases += int(line[2]) - int(line[1])
    # then, compute the chunk size
    chunk_size = int(round(total_bases / threads))
    # then, iterate over the regions and split if necessary
    # split to sizes of chunk_size. Split only, if the region is 30% longer than chunk_size.
    regions = []
    with open(bedfile, "r") as f:
        for l in f:
            line: list = l.strip().split("\t")
            # split if region is longer than chunk_size * 1.3
            if int(line[2]) - int(line[1]) < chunk_size * 1.3:
                regions.append((line[0], int(line[1]), int(line[2])))
            else:
                # split region
                start: int = int(line[1])
                while int(line[2]) > start:
                    step: int = min(chunk_size, int(line[2]) - start)
                    regions.append((line[0], start, start + step))
                    start: int = start + step
    return regions


# %%


def process_func(args):
    (
        region,
        path_bam,
        tmp_path,
        cache_size,
        min_overlap,
        min_fragment_size,
        max_inversion_coverage,
    ) = args
    # filter reads
    filtered_reads = find_candidates(
        region,
        path_bam,
        cache_size,
        min_overlap,
        min_fragment_size,
        max_inversion_coverage,
    )
    # write to tmp file
    filter_reads(region, filtered_reads, path_bam, tmp_path)
    return tmp_path


# %%


def filter_nonseparated(
    path_bam: Path,
    out_path: Path,
    path_bed: Path,
    cache_size: int,
    threads: int,
    min_overlap: float,
    min_fragment_size: int,
    max_inversion_coverage: float,
):
    # create regions
    regions = split_regions(path_bed, threads * 2)
    # sort regions by size, so that the biggest regions are first
    regions = sorted(regions, key=lambda x: x[2] - x[1], reverse=True)
    # create temp dir and tmp file paths for each region. The idea is that each thread
    # writes its output to another tmp file.
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_bam_paths = [Path(tmp_dir.name) / f"tmp_{i}.bam" for i in range(len(regions))]
    jobs_args = [
        (
            region,
            path_bam,
            tmp_bam_paths[i],
            cache_size,
            min_overlap,
            min_fragment_size,
            max_inversion_coverage,
        )
        for i, region in enumerate(regions)
    ]
    pool = mp.Pool(threads)
    pool.map(process_func, jobs_args)
    # cleanup tmp bam files. Keep only the ones that are not empty
    for i in range(-1, -len(tmp_bam_paths), -1):
        if tmp_bam_paths[i].stat().st_size == 0:
            tmp_bam_paths.pop(i)
    # merge tmp files and sort with samtools sort
    tmp_merged = Path(tmp_dir.name) / "tmp_merged.bam"
    cmd_merge = split(
        f"samtools merge --threads {threads} -f -o {str(tmp_merged)} {' '.join(list(map(str,tmp_bam_paths)))}"
    )
    cmd_sort = split(
        f"samtools sort --threads {threads} -o {str(out_path)} {str(tmp_merged)}"
    )
    cmd_index = split(f"samtools index -@ {threads} {str(out_path)}")
    log.info(f"merging bam files {str(tmp_bam_paths)}..")
    subprocess.run(cmd_merge)
    log.info(f"sorting bam file {str(out_path)}..")
    subprocess.run(cmd_sort)
    log.info(f"indexing bam file {str(out_path)}..")
    subprocess.run(cmd_index)


def test_argument_validity(args):
    # check params from args. min-overlap, max-inversion-coverage must be within the interval [0,1]
    if args.min_overlap < 0 or args.min_overlap > 1:
        raise ValueError("min-overlap must be within the interval [0,1].")
    if args.max_inversion_coverage < 0 or args.max_inversion_coverage > 1:
        raise ValueError("max-inversion-coverage must be within the interval [0,1].")
    # cache size must be at least 100
    if args.cache_size < 100:
        raise ValueError("cache size must be at least 100.")
    # min fragment size must be at least 100
    if args.min_fragment_size < 100:
        raise ValueError("min fragment size must be at least 100.")


def run(args, **kwargs):
    test_argument_validity(args)
    filter_nonseparated(
        path_bam=args.input,
        out_path=args.output,
        path_bed=args.bed,
        cache_size=args.cache_size,
        threads=args.threads,
        min_overlap=args.min_overlap,
        min_fragment_size=args.min_fragment_size,
        max_inversion_coverage=args.max_inversion_coverage,
    )


def get_parser():
    parser = argparse.ArgumentParser(description=".")
    parser.add_argument("-i", "--input", type=Path, help="path to input bam file")
    parser.add_argument("-o", "--output", type=Path, help="path to output bam file")
    parser.add_argument(
        "-b", "--bed", type=Path, help="path to bed file with regions to filter"
    )
    parser.add_argument("-t", "--threads", type=int, help="number of threads")
    parser.add_argument(
        "-c",
        "--cache-size",
        type=int,
        default=300,
        help="cache size. Shoul dbe around 10 times the depth of coverage.",
    )
    parser.add_argument(
        "-m",
        "--min-overlap",
        type=float,
        default=0.5,
        help="minimum overlap between two alignments to be considered as candidates.",
    )
    parser.add_argument(
        "-s",
        "--min-fragment-size",
        type=int,
        default=500,
        help="minimum size of fragment to be considered as candidate.",
    )
    parser.add_argument(
        "-v",
        "--max-inversion-coverage",
        type=float,
        default=0.1,
        help="maximum coverage by candidates to be considered as candidate. Increase if a greater proportion of the depth of coverage is attributed to candidates.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
