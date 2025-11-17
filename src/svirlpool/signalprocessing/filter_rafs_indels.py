# %%

import argparse
import csv
import json
import logging
import multiprocessing as mp
import shlex
import subprocess
import tempfile
from collections.abc import Generator
from pathlib import Path

import cattrs
import numpy as np

from ..signalprocessing import alignments_to_rafs
from ..util import datatypes, util

log = logging.getLogger(__name__)


# %%
def get_n_indels_signals(
    raf: datatypes.ReadAlignmentFragment, min_size: int, max_size: int
) -> int:
    return sum(
        1
        for sv in raf.SV_signals
        if (
            min_size <= sv.size <= max_size
            and sv.sv_type < 2
            and (
                raf.effective_interval[1] <= sv.ref_start <= raf.effective_interval[2]
                or raf.effective_interval[1] <= sv.ref_end <= raf.effective_interval[2]
            )
        )
    )


def get_indels_bp_per_1kb(
    raf: datatypes.ReadAlignmentFragment, min_size: int, max_size: int
) -> float:
    raf_effective_size = raf.effective_interval[2] - raf.effective_interval[1]
    n_indels = get_n_indels_signals(raf, min_size, max_size)
    if raf_effective_size == 0:
        return n_indels
    return n_indels / raf_effective_size * 1000


def read_rafs_per_chr(
    filename: str, chr: str
) -> Generator[datatypes.ReadAlignmentFragment, None, None]:
    cmd_tabix = f"tabix -0 -f -p bed {filename} {chr}"
    process = subprocess.Popen(shlex.split(cmd_tabix), stdout=subprocess.PIPE)
    for line in process.stdout:
        raw = (
            line.decode("utf-8")
            .strip()
            .split("\t")[-1]
            .strip()
            .replace('""', '"')[1:-1]
        )
        yield cattrs.structure(json.loads(raw), datatypes.ReadAlignmentFragment)


# %%
def filter_rafs_by_excessive_indel_counts_per_chr(
    input: Path,
    output: Path,
    chr: str,
    multiplier: float,
    dropped: Path | None,
    window_size: int = 30,
    percentile: int = 80,
) -> None:
    # Step 1: Read all RAFs and create list with (chr, start, end, read_name, n_deletions, raf_object, original_index)
    log.info(f"Reading all rafs from {chr} into memory.")
    rafs_data = []

    for idx, raf in enumerate(read_rafs_per_chr(input, chr)):
        # Count deletions (sv_type == 1)
        n_deletions = sum(1 for sv in raf.SV_signals if sv.sv_type == 1)
        rafs_data.append({
            "chr": raf.reference_name,
            "start": raf.effective_interval[1],
            "end": raf.effective_interval[2],
            "read_name": raf.read_name,
            "n_deletions": n_deletions,
            "raf": raf,
            "original_index": idx,
        })

    # Step 2: Sort by chr, start, end
    log.info(f"Sorting {len(rafs_data)} rafs for {chr}.")
    rafs_data.sort(key=lambda x: (x["chr"], x["start"], x["end"]))

    # Step 3-5: Sliding window to identify rafs exceeding threshold
    log.info(
        f"Filtering rafs in {chr} with window_size={window_size}, percentile={percentile}, multiplier={multiplier}."
    )
    filtered_indices = set()

    for i in range(len(rafs_data)):
        current_raf = rafs_data[i]
        current_start = current_raf["start"]
        current_end = current_raf["end"]

        # Define window boundaries
        window_start = max(0, i - window_size // 2)
        window_end = min(len(rafs_data), i + window_size // 2 + 1)

        # Get n_deletions values only for overlapping rafs in window
        window_deletions = []
        for j in range(window_start, window_end):
            if j == i:
                continue  # Skip the current raf itself
            other_raf = rafs_data[j]
            # Check if rafs overlap
            if current_start < other_raf["end"] and current_end > other_raf["start"]:
                window_deletions.append(other_raf["n_deletions"])

        # Calculate percentile (skip if not enough overlapping rafs)
        if len(window_deletions) < 3:
            continue

        percentile_value = np.percentile(window_deletions, percentile)

        # Check if current raf exceeds threshold
        if current_raf["n_deletions"] > multiplier * percentile_value:
            filtered_indices.add(i)

    # Step 6: Write results
    log.info(
        f"Writing filtered rafs for {chr}. Dropping {len(filtered_indices)} out of {len(rafs_data)} rafs."
    )
    with open(output, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        if dropped is not None:
            dropped_file_handle = open(dropped, "w")
            writer_dropped = csv.writer(dropped_file_handle, delimiter="\t")

        for i, raf_data in enumerate(rafs_data):
            raf = raf_data["raf"]
            chr = raf.reference_name
            start = raf.reference_alignment_start
            end = raf.reference_alignment_end

            if i in filtered_indices:
                if dropped:
                    writer_dropped.writerow([
                        chr,
                        start,
                        end,
                        json.dumps(raf.unstructure()),
                    ])
            else:
                writer.writerow([chr, start, end, json.dumps(raf.unstructure())])

        if dropped is not None:
            dropped_file_handle.close()


# def filter_rafs_by_excessive_indel_counts_per_chr(
#             input:Path,
#             output:Path,
#             chr:str,
#             rafs_cache_size:int,
#             multiplier:float,
#             dropped:Path|None) -> None:
#     rafs_cache_size_radius=int(rafs_cache_size/2)
#     log.info(f"filtering rafs in {chr} with rafs_cache_size_radius {rafs_cache_size_radius}.")
#     cache:list[datatypes.ReadAlignmentFragment] = []
#     cache_counts:list[float] = []
#     cache_abs_counts:list[int] = []
#     with open(output,'w') as f:
#         writer = csv.writer(f,delimiter='\t')
#         if dropped != None:
#             dropped_file_handle = open(dropped,'w')
#             writer_dropped = csv.writer(dropped_file_handle,delimiter='\t')
#         # with open(dropped,'w') as g:
#         #     writer_dropped = csv.writer(g,delimiter='\t')
#             for i_last,new_raf in enumerate(tqdm(read_rafs_per_chr(input,chr))):
#                 # fill cache until it has radius elements
#                 cache.append(new_raf)
#                 cache_counts.append(get_indels_bp_per_1kb(new_raf,8,10**6))
#                 cache_abs_counts = [get_n_indels_signals(raf,8,10**6) for raf in cache]
#                 if i_last >= rafs_cache_size_radius:
#                     selected_raf = cache[-rafs_cache_size_radius]
#                     selected_raf_effective_size = selected_raf.effective_interval[2] - selected_raf.effective_interval[1]
#                     # find all rafs in the cache that overlap with selected_raf
#                     indeces_overlapping_rafs = [i for i,raf in enumerate(cache)
#                         if selected_raf.referenceID == raf.referenceID
#                             and util.overlap((selected_raf.effective_interval[1:]),(raf.effective_interval[1:])) > selected_raf_effective_size*0.5]
#                     n_indels_in_overlapping_rafs = [cache_counts[i] for i in indeces_overlapping_rafs]
#                     # check if selected_raf has more indels than 2* the 90 percentile of the overlapping rafs
#                     if len(indeces_overlapping_rafs) > 5 \
#                             and cache_abs_counts[-rafs_cache_size_radius] > 5 \
#                             and cache_counts[-rafs_cache_size_radius] > 0.4 \
#                             and cache_counts[-rafs_cache_size_radius] > multiplier*np.percentile(n_indels_in_overlapping_rafs,90):
#                         if dropped:
#                             # write to dropped
#                             chr = selected_raf.reference_name
#                             start = selected_raf.reference_alignment_start
#                             end = selected_raf.reference_alignment_end
#                             writer_dropped.writerow([chr,start,end,json.dumps(selected_raf.unstructure())])
#                     else:
#                         # write to output
#                         chr = selected_raf.reference_name
#                         start = selected_raf.reference_alignment_start
#                         end = selected_raf.reference_alignment_end
#                         writer.writerow([chr,start,end,json.dumps(selected_raf.unstructure())])
#                     # pop first element from cache and cache_counts
#                     if len(cache) > 2*rafs_cache_size_radius:
#                         cache.pop(0)
#                         cache_counts.pop(0)
#                         cache_abs_counts.pop(0)
#             # iterate the last radius elements in cache. Don't change the cache, just advance selected_raf until it reaches the end
#             rafs_cache_size_radius = min(rafs_cache_size_radius,len(cache))
#             for i_mid, _ in enumerate(cache[-rafs_cache_size_radius:]):
#                 # # debug start
#                 # debug_index = -rafs_cache_size_radius+i_mid
#                 # if len(cache) < abs(debug_index):
#                 #     raise ValueError(f"cache has size {len(cache)}. debug_index = {debug_index}. i_mid = {i_mid}")
#                 # # debug end

#                 selected_raf = cache[-rafs_cache_size_radius+i_mid]
#                 selected_raf_effective_size = selected_raf.effective_interval[2] - selected_raf.effective_interval[1]
#                 indeces_overlapping_rafs = [i for i,raf in enumerate(cache)
#                     if selected_raf != raf
#                     and selected_raf.referenceID == raf.referenceID
#                         and util.overlap((selected_raf.effective_interval[1:]),(raf.effective_interval[1:])) > selected_raf_effective_size*0.5]
#                 n_indels_in_overlapping_rafs = [cache_counts[i] for i in indeces_overlapping_rafs]


#                 if len(indeces_overlapping_rafs) > 5 \
#                         and cache_abs_counts[-rafs_cache_size_radius+i_mid] > 5 \
#                         and cache_counts[-rafs_cache_size_radius+i_mid] > 0.4 \
#                         and cache_counts[-rafs_cache_size_radius+i_mid] > multiplier*np.percentile(n_indels_in_overlapping_rafs,90):
#                     if dropped:
#                         # write to dropped
#                         chr = selected_raf.reference_name
#                         start = selected_raf.reference_alignment_start
#                         end = selected_raf.reference_alignment_end
#                         writer_dropped.writerow([chr,start,end,json.dumps(selected_raf.unstructure())])
#                 else:
#                     # write to output
#                     chr = selected_raf.reference_name
#                     start = selected_raf.reference_alignment_start
#                     end = selected_raf.reference_alignment_end
#                     writer.writerow([chr,start,end,json.dumps(selected_raf.unstructure())])

#         if dropped != None:
#             dropped_file_handle.close()
#     log.info(f"filtered rafs in {chr}.")


def get_chr_names(reference: Path) -> list[str]:
    fasta_index = str(reference) + ".fai"
    with open(fasta_index) as f:
        return [line.split("\t")[0] for line in f]


def mp_process_filter_rafs_by_excessive_indel_counts_per_chr(args: dict) -> None:
    filter_rafs_by_excessive_indel_counts_per_chr(**args)


def filter_rafs_by_excessive_indel_counts(
    input: Path,
    output: Path,
    reference: Path,
    dropped_path: Path | None,
    threads: int,
    multiplier: float,
    window_size: int,
    percentile: int,
    tmp_dir_path: None | Path,
) -> None:
    multiplier = float(multiplier)
    if multiplier < 1.0:
        raise ValueError("multiplier must be >= 1.0")
    if multiplier < 2.0:
        log.warning(
            "multiplier is lower than 2.0. This will result in a lot of dropped rafs."
        )
    keep_dropped = dropped_path is not None
    # get chr names
    chr_names = get_chr_names(reference)
    # run filter_rafs_by_excessive_indel_counts_per_chr in parallel
    # create tmp output paths for each chr
    tmp_output_paths = [
        tempfile.NamedTemporaryFile(
            dir=tmp_dir_path,
            delete=False if tmp_dir_path else True,
            suffix=f".tmp_out.{chr}.tsv",
        )
        for chr in chr_names
    ]
    if keep_dropped:
        tmp_dropped_paths = [
            tempfile.NamedTemporaryFile(
                dir=tmp_dir_path,
                delete=False if tmp_dir_path else True,
                suffix=f".tmp_dropped.{chr}.tsv",
            )
            for chr in chr_names
        ]

    jobs_args = [
        {
            "input": input,
            "output": tmp_output_paths[i].name,
            "chr": chr,
            "dropped": tmp_dropped_paths[i].name if keep_dropped else None,
            "multiplier": multiplier,
            "window_size": window_size,
            "percentile": percentile,
        }
        for i, chr in enumerate(chr_names)
    ]

    with mp.Pool(min(threads, len(jobs_args))) as pool:
        pool.map(mp_process_filter_rafs_by_excessive_indel_counts_per_chr, jobs_args)

    # prepare genome file
    tmp_bedtools_genome = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path,
        delete=False if tmp_dir_path else True,
        suffix=".tmp_bedtools_genome",
    )
    util.genome_file_for_bedtools(reference=reference, output=tmp_bedtools_genome.name)
    # cat all tmp outputs to one output
    # and sort, compress and index the output
    log.info(f"Concatenating and compressing output to {output}")
    tmp_concatenated_output = tempfile.NamedTemporaryFile(
        delete=False if tmp_dir_path else True,
        dir=tmp_dir_path,
        suffix=".tmp_concatenated_output.tsv",
    )
    cmd_cat = f"cat {' '.join([tmp_output.name for tmp_output in tmp_output_paths])}"
    with open(tmp_concatenated_output.name, "w") as f:
        subprocess.check_call(shlex.split(cmd_cat), stdout=f)
    alignments_to_rafs.compress_and_index_bedlike(
        input=tmp_concatenated_output.name,
        output=output,
        threads=threads,
        genome=tmp_bedtools_genome.name,
    )

    if keep_dropped:
        log.info(f"Concatenating and compressing dropped rafs to {dropped_path}")
        tmp_concatenated_dropped = tempfile.NamedTemporaryFile(
            delete=False if tmp_dir_path else True,
            dir=tmp_dir_path,
            suffix=".tmp_concatenated_dropped.tsv",
        )
        cmd_cat = (
            f"cat {' '.join([tmp_dropped.name for tmp_dropped in tmp_dropped_paths])}"
        )
        with open(tmp_concatenated_dropped.name, "w") as f:
            subprocess.check_call(shlex.split(cmd_cat), stdout=f)
        alignments_to_rafs.compress_and_index_bedlike(
            input=tmp_concatenated_dropped.name,
            output=dropped_path,
            threads=threads,
            genome=tmp_bedtools_genome.name,
        )


# %%


def run(args, **kwargs):
    filter_rafs_by_excessive_indel_counts(
        input=args.input,
        output=args.output,
        reference=args.reference,
        dropped_path=args.dropped,
        threads=args.threads,
        multiplier=args.multiplier,
        window_size=args.window_size,
        percentile=args.percentile,
        tmp_dir_path=args.tmpdir,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Filters read alignment fragments by excessive indel counts compared to rafs in the proximity."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="tsv.gz file of read alignment fragments with adjusted effective intervals.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="output tsv.gz file of filtered read alignment fragments.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="reference genome file (.fasta).",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=1,
        help="number of threads to use.",
    )
    parser.add_argument(
        "--tmpdir", required=False, default=None, help="temporary directory to use."
    )
    parser.add_argument(
        "--dropped",
        required=False,
        default=None,
        help="output tsv file of dropped read alignment fragments.",
    )
    parser.add_argument(
        "--multiplier",
        type=float,
        required=False,
        default=4,
        help="multiplier for the percentile of indel counts to filter excessive indel counts. The higher, the more indels per raf are tolerated.",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        required=False,
        default=30,
        help="sliding window size for percentile calculation.",
    )
    parser.add_argument(
        "--percentile",
        type=int,
        required=False,
        default=90,
        help="percentile value (0-100) for indel count threshold calculation.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
