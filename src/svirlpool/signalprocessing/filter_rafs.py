# %%
# this script takes as input the bgzipped rafs file (columns: chr, start, end, structured datatypes.ReadAlignmentFragment)
# it reduces all RAFs to intervals that are flanked by unique regions
# the procedure is:
#  1) filter all rafs that don't overlap any unique region and write filtered rafs to a new (tmp) file that is bgzipped and indexed
#  2) for all rafs, create a BED file of start (chr,start,start+1) positions. (=RAF_starts)
#  3) for all rafs, create a BED file of end (chr,end-1,end) positions. (=RAF_ends)
#  4) use bedtools closest to find for each RAF start the closest unique downstream region.
#  5) use bedtools closest to find for each RAF end the closest unique upstream region.
#  6) reduce each RAF effective region to start+distance_to_closest_downstream_unique_region and end-distance_to_closest_upstream_unique_region
#  7) filter all SV signals that are out of the reduced RAF region
#  8) write the reduced RAFs to a new file

# %%

import argparse
import csv
import gzip
import json
import multiprocessing as mp
import shlex
import subprocess
import tempfile
from pathlib import Path

from logzero import logger as log
from tqdm import tqdm

from ..signalprocessing import alignments_to_rafs
from ..util import datatypes, util

# %%

#  1) filter all rafs that don't overlap any unique region and write filtered rafs to a new (tmp) file that is bgzipped and indexed


def sv_signals_densities(svSignals: list[datatypes.SVsignal], radius: int) -> list[int]:
    """Returns"""
    if len(svSignals) == 0:
        return []
    # check if list is sorted
    if not all(
        svSignals[i].ref_start <= svSignals[i + 1].ref_start
        for i in range(len(svSignals) - 1)
    ):
        raise ValueError("List of SVsignals is not sorted by ref_start.")

    # build a mask. The mask keeps the density of each signal. The mask has the length of svSignals
    mask = [0] * len(svSignals)
    left: int = 0
    right: int = 0
    for mid, midSvSignal in enumerate(svSignals):
        # move right until right is out of radius (svSignals[right].ref_start - midSvSignal.ref_end > radius)
        while (
            right < len(svSignals)
            and svSignals[right].ref_start - midSvSignal.ref_end <= radius
        ):
            right += 1
        # right is now the the first signal that is out of radius
        # move left until left is in radius (midSvSignal.ref_start - svSignals[left].ref_end <= radius)
        while left < mid and midSvSignal.ref_start - svSignals[left].ref_end > radius:
            left += 1
        # left is now the first signal that is in the radius
        # sum all signals sizes of the same type as midSvSignal that are in the radius
        mask[mid] = sum(
            s.size for s in svSignals[left:right] if s.sv_type == midSvSignal.sv_type
        )
    # filter all signals that have a mask value below min_signal_bp
    return mask


def filter_signals_for_minimum_density(
    svSignals: list[datatypes.SVsignal], min_signal_bp: int, radius: int
) -> list[datatypes.SVsignal]:
    """filters any signal that fails to accumulate another min_signal_bp within radius. The returned list is sorted by ref_start."""
    densities = sv_signals_densities(
        svSignals=sorted(svSignals, key=lambda x: x.ref_start), radius=radius
    )
    return [
        svSignals[i] for i in range(len(svSignals)) if densities[i] >= min_signal_bp
    ]


# TODO: test rafs_filtered
def filter_rafs_for_fully_non_unique_overlaps(
    rafs_path: Path,
    genome: Path,
    unique_regions: Path,
    output: Path,
    threads: int,
    rafs_filtered: Path | None = None,
) -> None:
    cmd_filter = f"bedtools intersect -u -a {rafs_path} -b {unique_regions}"
    cmd_sort = f"bedtools sort -g {str(genome)} -i -"
    cmd_bgzip = f"bgzip -@ {threads} -f -c"
    cmd_index = f"tabix -p bed {output}"
    with open(output, "wb") as f:
        p0 = subprocess.Popen(shlex.split(cmd_filter), stdout=subprocess.PIPE)
        p1 = subprocess.Popen(
            shlex.split(cmd_sort), stdin=p0.stdout, stdout=subprocess.PIPE
        )
        p2 = subprocess.Popen(shlex.split(cmd_bgzip), stdin=p1.stdout, stdout=f)
        p2.communicate()
    subprocess.check_call(shlex.split(cmd_index))
    # check if output exists. Else raise error
    if not Path(output).exists():
        raise FileNotFoundError(
            f"Output file {output} not found. Something went wrong."
        )
    if rafs_filtered:
        # subtract output from rafs_path to get the filtered rafs
        cmd_subtract = f"bedtools subtract -a {rafs_path} -b {output}"
        with open(rafs_filtered, "wb") as f:
            p0 = subprocess.Popen(shlex.split(cmd_subtract), stdout=f)
            p0.communicate()


#  2) for all rafs, create a BED file of start (chr,start,start+1) positions. (=RAF_starts)
#  3) for all rafs, create a BED file of end (chr,end-1,end) positions. (=RAF_ends)


def create_raf_starts_ends_bed(
    rafs_path: Path, starts_out: Path, ends_out: Path, genome: Path | None
) -> None:
    """returns the end positions in sorted order!"""
    tmp_ends = tempfile.NamedTemporaryFile(delete=True, suffix=".bed")
    with gzip.open(rafs_path, "rt") as f:
        with open(starts_out, "w") as out_starts:
            with open(tmp_ends.name, "w") as out_ends:
                for line in f:
                    chr, start, end, _ = line.strip().split("\t")
                    out_starts.write(f"{chr}\t{start}\t{int(start) + 1}\n")
                    out_ends.write(f"{chr}\t{int(end) - 1}\t{end}\n")
    # the ends file is not necessarily sorted. Sort it.
    if genome:
        cmd_sort = f"bedtools sort -g {str(genome)} -i {tmp_ends.name}"
    else:
        cmd_sort = f"sort -k1,1 -k2,2n {tmp_ends.name}"
    with open(ends_out, "w") as f:
        subprocess.check_call(shlex.split(cmd_sort), stdout=f)
    return


#  4) use bedtools closest to find for each RAF start the closest unique downstream region.
#  5) use bedtools closest to find for each RAF end the closest unique upstream region.


def mp_process_compress_and_index_bedlike(kwargs) -> None:
    alignments_to_rafs.compress_and_index_bedlike(**kwargs)


def split_bed_file(
    input: Path, dirpath: Path, n_splits: int, threads: int, genome: Path | None = None
) -> list[Path]:
    """Split a BED file into n_splits files."""
    if not 0 <= threads <= n_splits:
        threads = n_splits
        log.warning(
            f"Number of threads must be between 0 and n_splits. Setting threads to n_splits ({n_splits})."
        )
    if n_splits == 1:
        return [input]
    # log.info(f"Splitting input ({str(input)}) file into {n_splits} files..")
    output_files_uncompressed = [
        Path(dirpath) / f"split.{i}.bed" for i in range(n_splits)
    ]
    with gzip.open(input, "rt") as f:
        output_handles = [open(f, "w") for f in output_files_uncompressed]
        for i, line in tqdm(enumerate(f)):
            output_handles[i % n_splits].write(line)
        for handle in output_handles:
            handle.close()
    # compress and gzip each output file
    jobs_args = [
        {
            "genome": Path(genome) if genome else None,
            "sort_numerically": True if not genome else False,
            "input": f,
            "output": str(f) + ".gz",
            "threads": threads,
        }
        for f in output_files_uncompressed
    ]
    with mp.Pool(threads) as pool:
        pool.map(mp_process_compress_and_index_bedlike, jobs_args, chunksize=1)
    return [Path(str(outf) + ".gz") for outf in output_files_uncompressed]


# def adjust_effective_interval_starts(
#         rafs_path:Path,
#         output_rafs:Path,
#         genome:Path,
#         unique_regions:Path,
#         tmp_dir_path:Path) -> None:
#     """This function updates the effective interval starts of RAFs to the closest downstream unique region.

#     Args:
#         rafs_path (Path): path to the sorted rafs file (columns: chr, ref_start, ref_end, RAF)
#         output_rafs (Path): path to sorted output RAFs file
#         genome (Path): bedtools genome file
#         unique_regions (Path): unique regions bed file (sorted)
#         tmp_dir_path (Path): path to a temporary directory
#     """
#     # workflow:
#     # get start positions of all RAFs and calc distance to closest unique region downstream ("-iu -D ref")
#     # update effective_interval start to new start position (start+distance)
#     # save updated rafs to output_rafs
#     pass


def find_closest_unique_to_starts_and_ends(
    starts_in: Path,
    ends_in: Path,
    unique_regions: Path,
    tmp_dir_path: Path,
    genome: Path,
) -> list[int | None, int | None]:
    """Find for each RAF start the closest unique downstream region\
and for each RAF end the closest unique upstream region and report their distances."""
    cmd_closest_starts = f"bedtools closest -g {str(genome)} -iu -D ref -a {starts_in} -b {Path(unique_regions)}"
    cmd_closest_ends = f"bedtools closest -g {str(genome)} -id -D ref -a {ends_in} -b {Path(unique_regions)}"
    # the reulting table has columns chr, start, end, chr_target, start_target, end_target, distance
    tmp_closest_starts = tempfile.NamedTemporaryFile(
        dir=Path(tmp_dir_path), delete=False, suffix=".closest_starts.bed"
    )
    tmp_closest_ends = tempfile.NamedTemporaryFile(
        dir=Path(tmp_dir_path), delete=False, suffix=".closest_ends.bed"
    )
    # log.info(f"Finding closest unique regions to RAF starts..")
    with open(tmp_closest_starts.name, "w") as f:
        subprocess.check_call(shlex.split(cmd_closest_starts), stdout=f)
    # log.info(f"Finding closest unique regions to RAF ends..")
    with open(tmp_closest_ends.name, "w") as f:
        subprocess.check_call(shlex.split(cmd_closest_ends), stdout=f)
    # check if output exists. Else raise error
    if not Path(tmp_closest_starts.name).exists():
        raise FileNotFoundError(
            f"Output file {tmp_closest_starts.name} not found. Something went wrong."
        )
    if not Path(tmp_closest_ends.name).exists():
        raise FileNotFoundError(
            f"Output file {tmp_closest_ends.name} not found. Something went wrong."
        )
    # generate a table with the lett and right distances (colum 7) to the closest unique region
    # log.info(f"parsing distance from RAF start end RAF end to list..")
    starts_distances = [
        (
            abs(int(line.strip().split("\t")[-1]))
            if line.strip().split("\t")[-2] != "-1"
            else None
        )
        for line in open(tmp_closest_starts.name, "r").readlines()
    ]
    ends_distances = [
        (
            abs(int(line.strip().split("\t")[-1]))
            if line.strip().split("\t")[-2] != "-1"
            else None
        )
        for line in open(tmp_closest_ends.name, "r").readlines()
    ]
    return list(zip(starts_distances, ends_distances, strict=True))


#  6) reduce each RAF effective region to start+distance_to_closest_downstream_unique_region and end-distance_to_closest_upstream_unique_region


def get_effective_start(
    reference_start: int, reference_end: int, start_distance: int | None
) -> int:
    if start_distance is None:
        new_start = int(reference_end)
    else:
        new_start = min((int(reference_start) + start_distance, reference_end))
    return new_start


def get_effective_ends(
    reference_start: int, reference_end: int, end_distance: int | None
) -> int:
    if end_distance is None:
        new_end = int(reference_start)
    else:
        new_end = max((int(reference_end) - end_distance, reference_start))
    return new_end


# def get_effective_start_end(
#         reference_start:int,
#         reference_end:int,
#         start_distance:int|None,
#         end_distance:int|None) -> tuple[int,int]:
#     """Reduce the RAF effective region to start+distance_to_closest_downstream_unique_region \
# and end-distance_to_closest_upstream_unique_region."""
#     if start_distance is None:
#         new_start = int(reference_end)
#     else:
#         new_start = min((int(reference_start) + start_distance, reference_end))
#     if end_distance is None:
#         new_end = int(reference_start)
#     else:
#         new_end = max((int(reference_end) - end_distance, reference_start))
#     return (new_start,new_end)


def rafs_set_effective_sizes(
    input_rafs: Path,
    output_rafs: Path,
    distances: list[tuple[int | None, int | None]],
    min_raf_size: int,
    genome: Path,
    threads: int,
    tmp_dir_path: Path,
    rafs_filtered: Path | None = None,
) -> None:
    """Reduce each RAF effective region to start+distance_to_closest_downstream_unique_region \
and end-distance_to_closest_upstream_unique_region. Filters any RAF that is reduced to min_raf_size or below. \
This operation only adjusts the effective_interval parameter.
output_rafs is sorted, bgzip compressed, and tabix indexed"""
    if min_raf_size < 0:
        raise ValueError("min_raf_size must be 0 or greater")
    # iterate over all rafs and distances and reduce the effective region
    # write directly to the output file
    counter_filtered_rafs = 0
    counter_passed_rafs = 0
    # start coords are in the same order as rafs.
    # end coords are in the order of raf end coords. Sort rafs again by end coords and bedtools with a trick:
    # set the columns in the file to chr,end,end+1,raf
    # then sort with bedtools and update effective intervals. write the updated rafs again to file (chr,start,end,raf)
    tmp_output_starts = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=False, suffix=".tmp_output_starts.tsv"
    )
    # log.debug(f"Writing effective intervals to {tmp_output_starts.name}..")
    with open(tmp_output_starts.name, "w") as out:
        writer = csv.writer(out, delimiter="\t")
        for raf, (start_distance, _) in zip(
            util.yield_from_raf(Path(input_rafs)), distances, strict=True
        ):
            chr: str = str(raf.reference_name)
            start: int = int(raf.reference_alignment_start)
            end: int = int(raf.reference_alignment_end)
            effective_start = get_effective_start(
                start_distance=start_distance, reference_end=end, reference_start=start
            )
            raf.effective_interval = (chr, effective_start, end)
            # write chr,end,end+1,raf to allow to sort with bedtools by end coordinates
            writer.writerow([chr, end, end + 1, json.dumps(raf.unstructure())])

    tmp_output_starts_sorted = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=False, suffix=".tmp_output_starts_sorted.tsv.gz"
    )
    alignments_to_rafs.compress_and_index_bedlike(
        genome=genome,
        input=Path(tmp_output_starts.name),
        output=Path(tmp_output_starts_sorted.name),
        threads=threads,
    )

    # cmd_sort = f"bedtools sort -g {str(genome)} -i {tmp_output_starts.name}"
    # log.debug(f"Sorting RAFs {tmp_output_starts_sorted.name} by end coordinates..")
    # with open(tmp_output_starts_sorted.name,'w') as f:
    #     subprocess.check_call(shlex.split(cmd_sort), stdout=f)

    rafs_filtered_handle = None
    if rafs_filtered:
        rafs_filtered_handle = open(rafs_filtered, "w")

    # log.debug(f"Reading sorted effective intervals from {tmp_output_starts_sorted.name}..")
    tmp_output = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=False, suffix=".tmp_output.tsv"
    )
    with open(tmp_output.name, "w") as out:
        writer = csv.writer(out, delimiter="\t")
        writer_rafs_filtered = None
        if rafs_filtered_handle:
            writer_rafs_filtered = csv.writer(rafs_filtered_handle, delimiter="\t")
        for raf, (_, end_distance) in zip(
            util.yield_from_raf(Path(tmp_output_starts_sorted.name)),
            distances,
            strict=True,
        ):
            chr: str = str(raf.reference_name)
            start: int = int(raf.reference_alignment_start)
            end: int = int(raf.reference_alignment_end)
            effective_start: int = raf.effective_interval[1]
            effective_end: int = get_effective_ends(
                end_distance=end_distance, reference_end=end, reference_start=start
            )
            raf.effective_interval = (chr, effective_start, effective_end)
            if effective_end - effective_start >= min_raf_size:
                writer.writerow([chr, start, end, json.dumps(raf.unstructure())])
                counter_passed_rafs += 1
            else:
                if writer_rafs_filtered:
                    writer_rafs_filtered.writerow([
                        chr,
                        start,
                        end,
                        json.dumps(raf.unstructure()),
                    ])
                counter_filtered_rafs += 1
                continue
    # file will get sorted by alignments_to_rafs.compress_and_index_bedlike
    # check if output exists. Else raise error
    if not Path(tmp_output.name).exists():
        # compare length of distances and line count of input_rafs
        line_count = sum(1 for _ in gzip.open(input_rafs, "rt"))
        if len(distances) != line_count:
            raise ValueError(
                f"Length of distances ({len(distances)}) does not match line count of input rafs ({line_count})."
            )
        raise FileNotFoundError(
            f"Not able to generate output file {tmp_output.name}. Please check the input file {input_rafs}."
        )
    # compress and index with bgzip and tabix
    alignments_to_rafs.compress_and_index_bedlike(
        genome=genome,
        input=Path(tmp_output.name),
        output=Path(output_rafs),
        threads=threads,
    )
    tmp_output.close()
    # log.info(f"Done. {counter_passed_rafs} passed RAFs. {counter_filtered_rafs} filtered RAFs.")


#  7) filter all SV signals that are out of the reduced RAF region (?)
#  8) write the reduced RAFs to a new file


def process_single_rafs_file(
    input: Path,
    genome: Path,
    unique_regions: Path,
    output: Path,
    min_raf_size: int,
    tmp_dir_path: Path,
    rafs_filtered_none_unique: Path | None = None,
    rafs_filtered_too_small: Path | None = None,
) -> None:
    """_summary_

    Args:
        input (Path): path to the sorted, compressed, and indexed RAFs file (TSV: chr, start, end, RAF)
        unique_regions (Path): path to the sorted unique regions BED file
        output (Path): path to the sorted output RAFs file
        min_raf_size (int): minimum span of a RAF's effective interval on the reference
        tmp_dir_path (Path): path to a temporary directory
    """
    # log.info(f"Filtering RAFs for fully non-unique overlaps..")
    tmp_rafs_with_overlap_unique = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=True, suffix=".rafs.with_overlap_unique.bed.gz"
    )
    # filter_rafs_for_fully_non_unique_overlaps
    filter_rafs_for_fully_non_unique_overlaps(
        rafs_path=Path(input),
        genome=genome,
        threads=1,
        unique_regions=unique_regions,
        output=tmp_rafs_with_overlap_unique.name,
        rafs_filtered=rafs_filtered_none_unique,
    )
    # create_raf_starts_ends_bed
    # log.info(f"Creating RAF starts and ends BED files..")
    tmp_starts = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=True, suffix=".starts.bed"
    )
    tmp_ends = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=True, suffix=".ends.bed"
    )
    create_raf_starts_ends_bed(
        rafs_path=tmp_rafs_with_overlap_unique.name,
        starts_out=tmp_starts.name,
        ends_out=tmp_ends.name,
        genome=genome,
    )
    # find_closest_unique_to_starts_and_ends
    # log.info(f"Finding closest unique regions to RAF starts and ends..")
    distances = find_closest_unique_to_starts_and_ends(
        starts_in=tmp_starts.name,
        ends_in=tmp_ends.name,
        unique_regions=unique_regions,
        tmp_dir_path=tmp_dir_path,
        genome=genome,
    )
    # rafs_set_effective_sizes
    # log.info(f"Setting effective sizes of RAFs..")
    rafs_set_effective_sizes(
        input_rafs=tmp_rafs_with_overlap_unique.name,
        output_rafs=output,
        distances=distances,
        min_raf_size=min_raf_size,
        genome=genome,
        tmp_dir_path=tmp_dir_path,
        threads=1,
        rafs_filtered=rafs_filtered_too_small,
    )


def mp_process_single_bed(kwargs):
    process_single_rafs_file(**kwargs)


def adjust_rafs_effective_intervals(
    input_rafs: Path,
    reference: Path,
    output_rafs: Path,
    unique_regions: Path,
    threads: int,
    tmp_dir_path: Path,
    min_raf_size: int,
    rafs_dropped: Path | None = None,
) -> None:
    custom_tmp_dir: bool = True
    if tmp_dir_path is None:
        custom_tmp_dir = False
        tmp_dir = tempfile.TemporaryDirectory()
        tmp_dir_path = Path(tmp_dir.name)
        log.info(f"Created temporary directory {str(tmp_dir_path)}")
    # -1) make bedtools genome file
    log.info(f"creating genome file for bedtools at {str(tmp_dir_path)}..")
    tmp_genome = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=False, suffix=".genome"
    )
    util.genome_file_for_bedtools(reference=reference, output=tmp_genome.name)
    # 0) sort unique regions
    tmp_unique_regions = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=False, suffix=".unique_regions.sorted.bed.gz"
    )
    tmp_unique_regions_path = Path(tmp_unique_regions.name)
    log.info(
        f"Sorting, compressing, and indexing unique regions at {str(tmp_unique_regions_path)}.."
    )
    alignments_to_rafs.compress_and_index_bedlike(
        genome=tmp_genome.name,
        input=unique_regions,
        output=tmp_unique_regions_path,
        threads=threads,
    )
    # cmd_sort = f"sort -k1,1 -k2,2n {unique_regions}"
    # with open(tmp_unique_regions_path,'w') as f:
    #     subprocess.check_call(shlex.split(cmd_sort), stdout=f)
    # 1) split bed files
    # split_paths are sorted, bgzip compressed, and tabix indexed
    log.info(f"Splitting RAFs into {threads} files..")
    split_paths = split_bed_file(
        input=input_rafs,
        genome=tmp_genome.name,
        dirpath=tmp_dir_path,
        n_splits=threads,
        threads=threads,
    )
    log.info(f"split RAFs into files at: {' '.join([str(p) for p in split_paths])}..")
    # 2) create jobs to be executed in parallel
    input_output_files = list(
        zip(
            split_paths,
            [tmp_dir_path / f"split.{i}.output.bed" for i in range(threads)],
            strict=True,
        )
    )
    jobs_args = [
        {
            "input": split_in,
            "unique_regions": tmp_unique_regions_path,
            "output": split_out,
            "min_raf_size": min_raf_size,
            "tmp_dir_path": Path(tmp_dir_path) / f"split.{i}",
            "genome": Path(tmp_genome.name),
            "rafs_filtered_none_unique": (
                Path(tmp_dir_path) / f"split.{i}.rafs.none_unique.bed"
            ),
            "rafs_filtered_too_small": (
                Path(tmp_dir_path) / f"split.{i}.rafs.too_small.bed"
            ),
        }
        for i, (split_in, split_out) in enumerate(input_output_files)
    ]
    # create temporary directory for each job if it does not exist
    for job_args in jobs_args:
        if not Path(job_args["tmp_dir_path"]).exists():
            Path(job_args["tmp_dir_path"]).mkdir()
    with mp.Pool(threads) as pool:
        pool.map(mp_process_single_bed, jobs_args, chunksize=1)
    # zcat all output files into one
    tmp_output = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=True, suffix=".tsv"
    )
    cmd_zcat = f"zcat {' '.join([str(ff[1]) for ff in input_output_files])}"
    log.debug(
        f"Zcatting all output files into one: {tmp_output.name}. cmd is: {cmd_zcat}"
    )
    with open(tmp_output.name, "w") as f:
        subprocess.check_call(shlex.split(cmd_zcat), stdout=f)
    # compress and index with bgzip and tabix to output
    log.info("Compressing and indexing final output file..")
    alignments_to_rafs.compress_and_index_bedlike(
        genome=tmp_genome.name,
        input=Path(tmp_output.name),
        output=output_rafs,
        threads=threads,
    )
    if rafs_dropped:
        # cat all rafs_filtered_none_unique and rafs_filtered_too_small files into one
        # sort,compress,index to rafs_filtered
        log.info("Concatenating all filtered RAFs into one file..")
        tmp_filtered = tempfile.NamedTemporaryFile(
            dir=tmp_dir_path, delete=False, suffix=".filtered.tsv"
        )
        cmd_cat = f"cat {' '.join([str(fp) for fp in [job_args['rafs_filtered_none_unique'] for job_args in jobs_args] + [job_args['rafs_filtered_too_small'] for job_args in jobs_args]])}"
        with open(tmp_filtered.name, "w") as f:
            subprocess.check_call(shlex.split(cmd_cat), stdout=f)
        log.info("Compressing and indexing filtered RAFs file..")
        alignments_to_rafs.compress_and_index_bedlike(
            genome=tmp_genome.name,
            input=Path(tmp_filtered.name),
            output=rafs_dropped,
            threads=threads,
        )
    if not custom_tmp_dir:
        log.info(f"Cleaning up temporary directory {tmp_dir_path}..")
        tmp_dir.cleanup()
    log.info("Done.")


def run(args, **kwargs):
    adjust_rafs_effective_intervals(
        input_rafs=args.input,
        reference=args.reference,
        output_rafs=args.output,
        unique_regions=args.unique_regions,
        threads=args.threads,
        min_raf_size=args.min_raf_size,
        tmp_dir_path=args.tmp_dir,
        rafs_dropped=args.rafs_dropped,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Filter RAFs for fully non-unique overlaps and adjust effective intervals to unique regions."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the input RAFs bgzipped file, e.g. rafs.tsv.gz",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to the reference genome FASTA file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output RAFs bgzipped file, e.g. rafs.filtered.gz",
    )
    parser.add_argument(
        "-u",
        "--unique_regions",
        type=Path,
        required=True,
        help="Path to the unique regions BED file",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=8,
        help="Number of threads to use.",
    )
    parser.add_argument(
        "--min-raf-size",
        type=int,
        required=False,
        default=300,
        help="Minimum size of a RAF's effective interval. RAFs with less will get filtered.",
    )
    parser.add_argument(
        "--tmp-dir",
        required=False,
        default=None,
        help="Path to the temporary directory.",
    )
    parser.add_argument(
        "--rafs-dropped",
        type=Path,
        required=False,
        default=None,
        help="Path to the filtered RAFs bgzipped file, e.g. rafs.dropped.gz",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
