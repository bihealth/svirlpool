# extract signals from alignments directly
# %%
import argparse
import csv
import json
import multiprocessing as mp
import subprocess
import tempfile
import typing
from pathlib import Path
from shlex import split

import numpy as np
import pysam
from logzero import logger as log
from tqdm import tqdm

from . import datatypes, filter_nonseparated, filter_rafs, util

# %%

# path_regions = Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/data/d03.F.regions.bed")
# path_alignments = Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/alignments/minimap2.ont-r10.32x.bam")
# path_output = Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d01/HG002.signal")

# %%


def compress_and_index_bedlike(
    input: Path,
    output: Path,
    threads: int,
    genome: Path | None = None,
    sort_numerically: bool = False,
) -> None:
    """compress and index a bedlike file with bgzip and tabix"""
    if genome and sort_numerically:
        raise ValueError("Cannot use both genome and sort_numerically.")
    # if the first line of input has more than 11 columns, raise an error!
    # if len(open(input).readline().strip().split("\t")) > 11:
    #     raise ValueError(f"Input file {input} has more than 11 columns. This is not allowed.")

    cmd_sort = f"cat {Path(input)}"
    if sort_numerically:
        cmd_sort = f"sort -k1,1V -k2,2n -k3,3n {Path(input)}"
    if genome:
        cmd_sort = f"bedtools sort -g {str(genome)} -i {Path(input)}"
    cmd_bgzip = f"bgzip -@ {max(1,threads)} -f -c"
    cmd_index = f"tabix -0 -s 1 -b 2 -e 3 {Path(output)}"
    with open(Path(output), "wb") as file_out:
        p0 = subprocess.Popen(split(cmd_sort), stdout=subprocess.PIPE)
        p1 = subprocess.Popen(split(cmd_bgzip), stdin=p0.stdout, stdout=file_out)
        p1.communicate()
    subprocess.check_call(split(cmd_index))
    # check if both output and index exist, else raise error
    if not Path(output).exists():
        raise FileNotFoundError(
            f"Output file {output} not found. Called cmds: 1) {cmd_bgzip}."
        )
    if not Path(str(Path(output)) + "tbi").exists:
        raise FileNotFoundError(
            f"Index file {output}.tbi not found. Called cmds: 1) {cmd_index}."
        )


# parse regions
def parse_and_split_regions(
    path_regions: Path, desired_region_size: int = 10_000_000
) -> list[tuple[str, int, int]]:
    regions = []
    with open(path_regions, "r") as file:
        for line in file:
            chrom, start, end = line.strip().split("\t")[:3]
            start, end = int(start), int(end)
            regions.append((chrom, start, end))
    # then, evaluate if the regions need to be split further. Desired is a region size of 10_000_000
    # 1) find the regions that are already smaller than desired_region_size
    smaller_regions = [
        region for region in regions if region[2] - region[1] <= desired_region_size
    ]
    larger_regions = [
        region for region in regions if region[2] - region[1] > desired_region_size
    ]
    split_regions = []
    # split the larger regions into smaller regions of size desired_region_size
    for region in larger_regions:
        chrom, start, end = region
        N_regions = (end - start) // desired_region_size
        for i in range(N_regions):
            split_regions.append(
                (
                    chrom,
                    start + i * desired_region_size,
                    start + (i + 1) * desired_region_size,
                )
            )
        # add the last region
        split_regions.append((chrom, start + N_regions * desired_region_size, end))
    # concatenate the smaller regions with the split regions
    regions = smaller_regions + split_regions
    return regions


# %%
# function to extract signals
def get_indels(
    alignment: pysam.AlignedSegment,
) -> tuple[
    list[tuple[int, int]],
    list[tuple[int, int]],
    list[tuple[int, int]],
    list[tuple[int, int]],
]:
    """computes the indels from a pysam AlignedSegment object. That is:
    - deletions in reference
    - deletions in read
    - insertions in reference
    - insertions in read"""
    is_reverse = alignment.is_reverse
    ref_start = alignment.reference_start
    t, x = zip(*alignment.cigartuples)
    x = np.array(x)
    t = np.array(t)
    x_read_starts, x_read_ends, x_ref_starts, x_ref_ends = util.get_starts_ends(
        t=t, x=x, is_reverse=is_reverse, reference_start=ref_start
    )
    # get deletion ref start, ref end tuples
    deletions_ref = list(zip(x_ref_starts[(t == 2)], x_ref_ends[(t == 2)]))
    # get deletion read start, read end tuples
    deletions_read = list(zip(x_read_starts[(t == 2)], x_read_ends[(t == 2)]))
    # get insertion ref start, ref end tuples
    insertions_ref = list(zip(x_ref_starts[(t == 1)], x_ref_ends[(t == 1)]))
    # get insertion read start, read end tuples
    insertions_read = list(zip(x_read_starts[(t == 1)], x_read_ends[(t == 1)]))
    return deletions_ref, deletions_read, insertions_ref, insertions_read


def get_start_end(alignment: pysam.AlignedSegment):
    """computes the start and end positions of the alignment on the read and on the reference.
    - ref_start: start position on the reference
    - ref_end: end position on the reference
    - read_start: start position on the read
    - read_end: end position on the read"""
    if alignment.is_unmapped:
        raise ValueError(
            "Alignment is unmapped. Make sure unly mapped alignments are passed."
        )
    is_reverse = alignment.is_reverse
    ref_start = alignment.reference_start
    ref_end = alignment.reference_end
    t, x = zip(*alignment.cigartuples)
    x = np.array(x)
    t = np.array(t)
    m = sum(x[(t == 0) | (t == 1) | (t == 7) | (t == 8)])
    if is_reverse:
        read_start = (
            alignment.cigartuples[-1][1]
            if alignment.cigartuples[-1][0] in (4, 5)
            else 0
        )
        read_end = read_start + m
    else:
        read_start = (
            alignment.cigartuples[0][1] if alignment.cigartuples[0][0] in (4, 5) else 0
        )
        read_end = read_start + m
    return ref_start, ref_end, read_start, read_end


# %%


def parse_SVsignals_from_alignment(
    alignment: pysam.AlignedSegment,
    ref_start: int,
    ref_end: int,
    read_start: int,
    read_end: int,
    min_signal_size: int,
    min_bnd_size: int,
) -> list[datatypes.SVsignal]:
    """parse the SV signals from a pysam AlignedSegment object"""
    if alignment.is_reverse:
        read_start, read_end = read_end, read_start

    bndl = None
    if (
        alignment.cigartuples[0][0] in (4, 5)
        and alignment.cigartuples[0][1] >= min_bnd_size
    ):
        clipped_left_length = alignment.cigartuples[0][1]
        bndl = datatypes.SVsignal(
            ref_start=int(ref_start),
            ref_end=int(ref_start),
            read_start=int(read_start),
            read_end=int(read_start),
            size=int(clipped_left_length),
            sv_type=int(3),
        )

    bndr = None
    if (
        alignment.cigartuples[-1][0] in (4, 5)
        and alignment.cigartuples[-1][1] >= min_bnd_size
    ):
        clipped_right_length = alignment.cigartuples[-1][1]
        bndr = datatypes.SVsignal(
            ref_start=int(ref_end),
            ref_end=int(ref_end),
            read_start=int(read_end),
            read_end=int(read_end),
            size=int(clipped_right_length),
            sv_type=int(4),
        )

    sv_signals = []
    if bndl:
        sv_signals.append(bndl)
    if bndr:
        sv_signals.append(bndr)

    deletions_ref, deletions_read, insertions_ref, insertions_read = get_indels(
        alignment
    )
    deletions = [
        datatypes.SVsignal(
            ref_start=int(dell),
            ref_end=int(delr),
            read_start=int(rdell),
            read_end=int(rdelr),
            size=int(abs(delr - dell)),
            sv_type=int(1),
        )
        for ((dell, delr), (rdell, rdelr)) in zip(deletions_ref, deletions_read)
        if abs(delr - dell) >= min_signal_size
    ]
    insertions = [
        datatypes.SVsignal(
            ref_start=int(insl),
            ref_end=int(insr),
            read_start=int(rinsl),
            read_end=int(rinsr),
            size=int(abs(rinsr - rinsl)),
            sv_type=int(0),
        )
        for ((insl, insr), (rinsl, rinsr)) in zip(insertions_ref, insertions_read)
        if abs(rinsr - rinsl) >= min_signal_size
    ]
    if len(deletions):
        sv_signals.extend(deletions)
    if len(insertions):
        sv_signals.extend(insertions)
    return sv_signals


def parse_ReadAlignmentFragment_from_alignment(
    alignment: pysam.AlignedSegment,
    samplename: str,
    min_signal_size: int,
    min_bnd_size: int,
    filter_density_radius: int = 0,
    filter_density_min_bp: int = 0,
) -> datatypes.ReadAlignmentFragment:

    ref_start, ref_end, read_start, read_end = get_start_end(alignment)

    sv_signals = sorted(
        parse_SVsignals_from_alignment(
            alignment=alignment,
            ref_start=ref_start,
            ref_end=ref_end,
            read_start=read_start,
            read_end=read_end,
            min_signal_size=min_signal_size,
            min_bnd_size=min_bnd_size,
        )
    )

    sv_signals_density_filtered = sv_signals
    if filter_density_radius > 0 and filter_density_min_bp > 0:
        sv_signals_density_filtered = sorted(
            filter_rafs.filter_signals_for_minimum_density(
                svSignals=sv_signals,
                min_signal_bp=filter_density_min_bp,
                radius=filter_density_radius,
            )
        )

    raf = datatypes.ReadAlignmentFragment(
        samplename=samplename,
        read_name=str(alignment.query_name),
        reference_name=str(alignment.reference_name),
        referenceID=int(alignment.reference_id),
        reference_alignment_start=int(ref_start),
        reference_alignment_end=int(ref_end),
        read_alignment_start=int(read_start),
        read_alignment_end=int(read_end),
        mapping_quality=int(alignment.mapping_quality),
        alignment_forward=bool(not alignment.is_reverse),
        alignment_secondary=bool(alignment.is_secondary),
        alignment_supplementary=bool(alignment.is_supplementary),
        effective_interval=(
            str(alignment.reference_name),
            int(ref_start),
            int(ref_end),
        ),
        SV_signals=sv_signals_density_filtered,
    )
    return raf


def process_region(
    path_alignments: Path,
    region: typing.Tuple[str, int, int],
    samplename: str,
    min_signal_size: int,
    min_bnd_size: int,
    min_mapq: int,
    min_segment_size: int,
    filter_density_radius: int,
    filter_density_min_bp: int,
    filter_nonseparated__cache_size_alignments_filtering: int,
    filter_nonseparated__min_overlap: float,
    filter_nonseparated__min_fragment_size: int,
    filter_nonseparated__max_inversion_coverage: float,
) -> typing.List[datatypes.ReadAlignmentFragment]:
    """parse all alignments in a region and create from each alignment a ReadAlignmentFragment"""
    chrom, start, end = region
    N_alignments = 0
    # log.debug(f"process region {chrom}:{start}-{end}")
    non_separated_reads = set()
    if filter_nonseparated__cache_size_alignments_filtering > 0:
        # at first, find all non-separated ONT reads that should be filtered
        non_separated_reads: set = filter_nonseparated.find_candidates(
            bam_path=path_alignments,
            region=region,
            cache_size=filter_nonseparated__cache_size_alignments_filtering,
            min_overlap=filter_nonseparated__min_overlap,
            min_fragment_size=filter_nonseparated__min_fragment_size,
            max_inversion_coverage=filter_nonseparated__max_inversion_coverage,
        )
    lines = []
    with pysam.AlignmentFile(path_alignments, "rb") as file:
        for alignment in file.fetch(chrom, start, end):
            N_alignments += 1
            if alignment.query_name in non_separated_reads:
                continue
            if (
                alignment.mapping_quality < min_mapq
                or alignment.is_unmapped
                or alignment.is_secondary
                or alignment.reference_end - alignment.reference_start
                < min_segment_size
                or alignment.reference_end < start
                or alignment.reference_start > end
            ):
                continue
            lines.append(
                parse_ReadAlignmentFragment_from_alignment(
                    alignment=alignment,
                    samplename=samplename,
                    min_bnd_size=min_bnd_size,
                    min_signal_size=min_signal_size,
                    filter_density_radius=filter_density_radius,
                    filter_density_min_bp=filter_density_min_bp,
                )
            )
            # filter any sv signal that is not within region
            lines[-1].SV_signals = [
                sv
                for sv in lines[-1].SV_signals
                if sv.ref_start < end and start < sv.ref_end
            ]
            # write line to tmp file
    return lines


def mp_process_region(kwargs) -> list:
    return process_region(**kwargs)


# %%


# not used right now
def filter_signals_for_too_many_aligned_segments(lines: list, max_bnds: int) -> list:
    counter = {}
    c = 0
    for line in lines:
        if line[7] not in counter:
            counter[line[7]] = 1
            c += 1
        else:
            counter[line[7]] += 1
    # count the number of values that are > max_bnds
    large_values = sum([1 for value in counter.values() if value > max_bnds])
    log.info(
        f"among {c} reads there are {large_values} reads with more than {max_bnds} BNDs which are filtered out."
    )
    return [line for line in lines if counter[line[7]] <= max_bnds]


def process_bam(
    path_alignments: Path | str,
    path_regions: Path | str,
    path_output: Path | str,
    samplename: int,
    min_signal_size: int,
    min_bnd_size: int,
    min_segment_size: int,
    min_mapq: int,
    filter_density_radius: int,
    filter_density_min_bp: int,
    filter_nonseparated__cache_size_alignments_filtering: int,
    filter_nonseparated__min_overlap: float,
    filter_nonseparated__min_fragment_size: int,
    filter_nonseparated__max_inversion_coverage: float,
    threads: int = 1,
    tmp_dir_path: Path | None = None,
):

    if threads < 1:
        threads = mp.cpu_count()
    if threads > mp.cpu_count():
        log.warning(
            f"threads {threads} > cpu_count {mp.cpu_count}. Setting threads to {mp.cpu_count()}"
        )
        threads = mp.cpu_count()

    log.info(f"parse regions from {path_regions}")
    regions = parse_and_split_regions(path_regions=path_regions)

    jobs = [
        {
            "path_alignments": path_alignments,
            "region": region,
            "samplename": samplename,
            "min_signal_size": min_signal_size,
            "min_bnd_size": min_bnd_size,
            "min_mapq": min_mapq,
            "min_segment_size": min_segment_size,
            "filter_density_radius": filter_density_radius,
            "filter_density_min_bp": filter_density_min_bp,
            "filter_nonseparated__cache_size_alignments_filtering": (
                filter_nonseparated__cache_size_alignments_filtering
            ),
            "filter_nonseparated__min_overlap": filter_nonseparated__min_overlap,
            "filter_nonseparated__min_fragment_size": (
                filter_nonseparated__min_fragment_size
            ),
            "filter_nonseparated__max_inversion_coverage": (
                filter_nonseparated__max_inversion_coverage
            ),
        }
        for i, region in enumerate(regions)
    ]

    log.info(
        f"execute {len(jobs)} jobs on {len(regions)} regions in parallel on {threads} threads.."
    )
    if threads > 1:
        with mp.Pool(threads) as pool:
            # adjust to saving the list outputs of all jobs
            job_results: typing.List[datatypes.ReadAlignmentFragment] = [
                raf
                for chunk in tqdm(
                    pool.imap(mp_process_region, jobs, chunksize=1), total=len(jobs)
                )
                for raf in chunk
            ]
    else:
        job_results = [
            raf for chunk in [process_region(**job) for job in jobs] for raf in chunk
        ]

    log.info(f"write job results to {path_output}")
    tmp_out = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=(False if tmp_dir_path else True), suffix=".tsv"
    )
    with open(tmp_out.name, "w") as file_out:
        writer = csv.writer(file_out, delimiter="\t")
        for raf in sorted(
            job_results,
            key=lambda raf: (
                raf.referenceID,
                raf.reference_alignment_start,
                raf.reference_alignment_end,
            ),
        ):
            writer.writerow(
                [
                    raf.reference_name,
                    raf.reference_alignment_start,
                    raf.reference_alignment_end,
                    json.dumps(raf.unstructure()),
                ]
            )
    # compress and index the output file
    compress_and_index_bedlike(input=tmp_out.name, output=path_output, threads=threads)
    log.info(f"written tmp output of sample {samplename} to {path_output}")


# %%

# def process_all_samples(
#             path_sampledicts:Path,
#             path_regions:Path,
#             output_prefix:Path,
#             min_signal_size:int,
#             min_bnd_size:int,
#             min_segment_size:int,
#             min_mapq:int,
#             threads:int,
#             filter_density_radius:int,
#             filter_density_min_bp:int,
#             tmp_dir_path:Path|None=None):
#     log.info(f"load sample dicts from {path_sampledicts}")
#     sampleID_bam_dict = util.load_sampledicts(path_sampledicts)[2]

#     # # iterate over all samples and store their outputs in a temporary directory
#     # if not tmp_dir_path:
#     #     tmp_dir = tempfile.TemporaryDirectory()
#     #     tmp_dir_path = Path(tmp_dir.name)

#     log.info(f"process all samples")
#     for sampleID, path_alignments in sampleID_bam_dict.items():
#         path_output = str(output_prefix)+f".{sampleID}.rafs.tsv.gz"
#         process_bam(
#             path_alignments=path_alignments,
#             path_regions=path_regions,
#             path_output=path_output,
#             sampleID=sampleID,
#             min_signal_size=min_signal_size,
#             min_bnd_size=min_bnd_size,
#             min_segment_size=min_segment_size,
#             min_mapq=min_mapq,
#             threads=threads,
#             filter_density_radius=filter_density_radius,
#             filter_density_min_bp=filter_density_min_bp,
#             tmp_dir_path=tmp_dir_path)

# log.info("concatenate and sort all sample outputs")
# # TODO: each sample output is already sorted. So, we can merge them in a sorted manner to avoid the quadratic complexity of sorting.
# cmd_cat = "cat " + ' '.join([str(tmp_dir_path / f"{sampleID}.tsv") for sampleID in sampleID_bam_dict.keys()])
# cmd_sort = "sort -k1,1 -k2,2n"
# tmp_text_out = tmp_dir_path / "tmp_text_out.tsv"
# with open(tmp_text_out, "w") as file_out:
#     p0 = subprocess.Popen(split(cmd_cat),stdout=subprocess.PIPE)
#     p1 = subprocess.Popen(split(cmd_sort),stdin=p0.stdout,stdout=file_out)
#     p1.communicate()
# # remove tmp files ([str(tmp_dir_path / f"{sampleID}.tsv") for sampleID in sampleID_bam_dict.keys()])
# for sampleID in sampleID_bam_dict.keys():
#     (tmp_dir_path / f"{sampleID}.tsv").unlink()

# compress_and_index_bedlike(input=tmp_text_out,output=output,threads=threads)
# # remove tmp file
# tmp_text_out.unlink()
# log.info(f"done. Output written to {output}")


def run(args, **kwargs):
    process_bam(
        path_alignments=args.alignments,
        path_regions=args.regions,
        path_output=args.output,
        samplename=args.samplename,
        min_signal_size=args.min_signal_size,
        min_bnd_size=args.min_bnd_size,
        min_segment_size=args.min_segment_size,
        min_mapq=args.min_mapq,
        filter_density_radius=args.filter_density_radius,
        filter_density_min_bp=args.filter_density_min_bp,
        threads=args.threads,
        filter_nonseparated__cache_size_alignments_filtering=args.filter_nonseparated__cache_size_alignments_filtering,
        filter_nonseparated__min_overlap=args.filter_nonseparated__min_overlap,
        filter_nonseparated__min_fragment_size=args.filter_nonseparated__min_fragment_size,
        filter_nonseparated__max_inversion_coverage=args.filter_nonseparated__max_inversion_coverage,
        tmp_dir_path=args.tmp_dir_path,
    )

    # path_alignments:Path,
    # path_regions:Path,
    # path_output:Path,
    # sampleID:int,
    # min_signal_size:int,
    # min_bnd_size:int,
    # min_segment_size:int,
    # min_mapq:int,
    # filter_density_radius:int,
    # filter_density_min_bp:int,
    # threads:int=1,
    # tmp_dir_path:Path|None=None):

    # process_all_samples(
    #     path_sampledicts=args.sampledicts,
    #     path_regions=args.regions,
    #     output_prefix=args.output_prefix,
    #     min_signal_size=args.min_signal_size,
    #     min_bnd_size=args.min_bnd_size,
    #     min_segment_size=args.min_segment_size,
    #     min_mapq=args.min_mapq,
    #     threads=args.threads,
    #     tmp_dir_path=args.tmp_dir_path,
    #     filter_density_radius=args.filter_density_radius,
    #     filter_density_min_bp=args.filter_density_min_bp)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Extract read alignment fragment objects from a provided bam file."
    )
    parser.add_argument(
        "-a", "--alignments", type=Path, required=True, help="Path to the bam file."
    )
    parser.add_argument(
        "-s",
        "--samplename",
        type=str,
        required=True,
        help="Name of the sample. Should be unique in your whole study.",
    )
    parser.add_argument(
        "-r",
        "--regions",
        type=Path,
        required=True,
        help="Path to the regions bed file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output sorted and compressed tsv file.",
    )
    parser.add_argument(
        "--min-signal-size",
        type=int,
        default=7,
        help="Minimum size of an indel to be considered in any case. Default is 7.",
    )
    parser.add_argument(
        "--min-bnd-size",
        type=int,
        default=150,
        help="Minimum size of a clipped end to be considered. Default is 150.",
    )
    parser.add_argument(
        "--min-segment-size",
        type=int,
        default=300,
        help="Minimum size of a segment to be considered. Default is 300.",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=1,
        help="Minimum mapping quality to consider an alignment.",
    )
    parser.add_argument(
        "--filter-density-radius",
        type=int,
        default=500,
        help="The radius in which the density of similar sv signals is calcualted. Default is 500.",
    )
    parser.add_argument(
        "--filter-density-min-bp",
        type=int,
        default=30,
        help="The minimum number of summed bp of the same SV type in radius. e.g. an insertion of 25bp and no other insertion within a 500bp radius will get filtered.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=8,
        help="Number of threads to use. Default is 8.",
    )
    parser.add_argument(
        "--filter-nonseparated__cache-size-alignments-filtering",
        type=int,
        default=1000,
        help="Cache size for the non-separated read filtering. Default is 1000.",
    )
    parser.add_argument(
        "--filter-nonseparated__min-overlap",
        type=float,
        default=0.3,
        help="minimum overlap between two alignments to be considered as candidates. Default is 0.3",
    )
    parser.add_argument(
        "--filter-nonseparated__min-fragment-size",
        type=int,
        default=500,
        help="minimum fragment size of a read to be considered. Default is 500.",
    )
    parser.add_argument(
        "--filter-nonseparated__max-inversion-coverage",
        type=float,
        default=0.15,
        help="maximum coverage by candidates to be considered as candidate. Increase if a greater proportion of the depth of coverage is attributed to candidates. Default is 0.15",
    )
    parser.add_argument(
        "--tmp-dir-path",
        default=None,
        help="Path to a temporary directory to store intermediate results. Default is ''.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
