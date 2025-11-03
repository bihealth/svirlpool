# this script subsamples a fastq file.
# inputs is a fastq (gzipped or not) and a probability of keeping each read

import argparse
import gzip
import random
import sys
from pathlib import Path

from Bio import SeqIO
from logzero import logger as log
from tqdm import tqdm


def subsample_fastq(input: Path, prob: float, output: Path = Path("")):
    # test if prob is between 0 and 1
    if not 0 < prob <= 1:
        raise ValueError("Probability x must be 0 < x <= 1")
    counter = 0
    filtered = 0
    if input.suffix == ".gz":
        f = gzip.open(input, "rt")
        log.info(f"Reading gzipped file {input}")
    else:
        f = open(input, "r")
        log.info(f"Reading file {input}")
    with f:
        # if no output is give, print to stdout
        if output == Path(""):
            out = sys.stdout
            log.info("Writing to stdout")
        else:
            if output.suffix == ".gz":
                out = gzip.open(output, "wt")
                log.info(f"Writing to gzipped file {output}")
            else:
                out = open(output, "w")
                log.info(f"Writing to file {output}")
        with out:
            for record in tqdm(SeqIO.parse(f, "fastq")):
                counter += 1
                if random.random() < prob:
                    SeqIO.write(record, out, "fastq")
                else:
                    filtered += 1
    log.info(
        f"Kept {counter-filtered} reads out of {counter} ({(counter-filtered)/counter*100:.2f}%)"
    )


def run(args, **kwargs):
    subsample_fastq(input=args.input, prob=args.prob, output=args.output)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Subsample a fastq file by a given probability"
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="path to (gzipped) fastq file"
    )
    parser.add_argument(
        "-p",
        "--prob",
        type=float,
        required=True,
        help="probability of keeping each read",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(""),
        help="output file (default: stdout)",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
