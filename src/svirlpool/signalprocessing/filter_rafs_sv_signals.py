# this script iterates over a rafs file and creates counts on the sv signals of all rafs.
# %%
import argparse
import csv
import json
import tempfile
from pathlib import Path

from logzero import logger as log

from ..signalprocessing import alignments_to_rafs
from ..util import util
from . import rafs_indel_histograms

# %%


def filter_rafs_sv_signals(
    input: Path,
    reference: Path,
    output: Path,
    tmp_dir_path: Path = None,
    threads: int = 1,
    filter_threshold: float = 1.0,
):
    log.info(f"counting deletions and insertions in {input} ...")
    counts, n_rafs, sum_bp = rafs_indel_histograms.get_indel_counts(input=input)

    log.info(f"counted deletions and insertions in {n_rafs} RAFs")
    # iterate the sorted keys of counts["deletions"] and counts["insertions"]
    # if counts["deletions"][size] / n_rafs <= 0.03 or counts["insertions"][size] / n_rafs <= 0.03
    # break the loop and set min_deletion_size and min_insertion_size to size
    min_deletion_size = 999999
    min_insertion_size = 999999
    # deletions first
    for size in sorted(counts["deletions"].keys()):
        if counts["deletions"][size] / (sum_bp / 1e6) <= filter_threshold:
            min_deletion_size = max(size, 6)
            break
    # insertions second
    for size in sorted(counts["insertions"].keys()):
        if counts["insertions"][size] / (sum_bp / 1e6) <= filter_threshold:
            min_insertion_size = max(size, 6)
            break

    # iterate the input rafs again and filter out all deletions and insertions with size < min_deletion_size and min_insertion_size
    # create a new rafs file with the filtered rafs
    tmp_out = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=False if tmp_dir_path else True, suffix=".tsv"
    )

    log.info(
        f"Filtering deletions with size < {min_deletion_size} and insertions with size < {min_insertion_size} ..."
    )
    with open(tmp_out.name, "w") as out:
        writer = csv.writer(out, delimiter="\t")
        for raf in util.yield_from_raf(input=input):
            raf.SV_signals = [
                sv
                for sv in raf.SV_signals
                if sv.sv_type > 2
                or sv.sv_type == 0
                and sv.size >= min_insertion_size
                or sv.sv_type <= 2
                and sv.size >= min_deletion_size
            ]
            writer.writerow([
                raf.reference_name,
                raf.reference_alignment_start,
                raf.reference_alignment_end,
                json.dumps(raf.unstructure()),
            ])

    # compress and index bedlike
    tmp_genome = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path, delete=False if tmp_dir_path else True, suffix=".genome"
    )
    util.genome_file_for_bedtools(reference=reference, output=tmp_genome.name)

    alignments_to_rafs.compress_and_index_bedlike(
        genome=tmp_genome.name,
        input=Path(tmp_out.name),
        output=Path(output),
        threads=threads,
    )


def run(args, **kwargs):
    filter_rafs_sv_signals(
        input=args.input,
        reference=args.reference,
        output=args.output,
        tmp_dir_path=args.tmp_dir,
        threads=args.threads,
        filter_threshold=args.threshold,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Filter read alignment fragments (rafs) for gaps with a size that occurs in gaps more often than once per megabase."
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
        "--tmp-dir",
        required=False,
        default=None,
        help="Path to the temporary directory.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=1,
        help="Number of threads to use.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        required=False,
        default=1.0,
        help="Filter threshold for deletions and insertions. Default is 1.0. A higher value means more candidate regions and possible false positives. A lower value can make the tool blind towards smaller SVs or ripped SV signals.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
