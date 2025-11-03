# %%
import argparse
from pathlib import Path

import pysam
from logzero import logger as log
from tqdm import tqdm

from . import util

# %%

# path_alignments = Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d01/alignments/HG003.bam")
# path_out = Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d01/alignments/HG003.filtered.bam")


def filter_bam_for_incomplete_sequences(path_alignments: Path, path_out: Path):
    """
    Filter a bam file for reads that have incomplete sequences.
    A read is considered incomplete if there is no complete DNA sequence in any of its alignments or all alignments have hard clipped tails.
    The filtered bam file will contain only reads that have full sequences.
    :param path_alignments: Path to input bam file.
    :param path_out: Path to output bam file.
    """
    set_no_seq = set()
    # set_no_seq is the set of readIDs that might have no full sequence in the bam file

    set_with_seq = set()
    # set_with_seq is the set of readIDs that have full sequence in the bam file

    log.info("Reading the alignments and finding reads with incomplete sequences...")
    # at first, fill the sets
    with pysam.AlignmentFile(path_alignments, "rb") as input:
        for read in tqdm(input):
            if read.is_unmapped:
                continue
            if read.query_name in set_with_seq:
                continue
            if read.cigartuples[0][0] == 5 or read.cigartuples[-1][0] == 5:
                set_no_seq.add(read.query_name)
            else:
                if not read.query_sequence:
                    set_no_seq.add(read.query_name)
                    continue
                if util.query_total_length(aln=read) != len(read.query_sequence):
                    log.warning(
                        f"Read {read.query_name} has softclipping but total length does not match sequence length"
                    )
                set_with_seq.add(read.query_name)
                if read.query_name in set_no_seq:
                    set_no_seq.remove(read.query_name)

    log.info(f"Number of reads with no full sequence: {len(set_no_seq)}")
    log.info(f"Number of reads with full sequence: {len(set_with_seq)}")

    with pysam.AlignmentFile(path_out, "wb") as out:
        with pysam.AlignmentFile(path_alignments, "rb") as input:
            for read in tqdm(input):
                if read.is_unmapped:
                    continue
                if read.query_name in set_no_seq:
                    continue
                out.write(read)

    log.info(f"Filtered bam file written to {path_out}")


def run(args, **kwargs):
    filter_bam_for_incomplete_sequences(
        path_alignments=args.input, path_out=args.output
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Filter a bam file for reads that have incomplete sequences."
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="path to input bam file"
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="path to output bam file"
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
