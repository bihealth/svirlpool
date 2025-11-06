# script to filter a signaldepth file
# in a first pass, the statistics of the depth for each sample on each chromosome is calculated
# in a second pass, all positions with a depth higher than N*mean are filtered out
# all not filtered positions are written to a new file
# %%
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from logzero import logger as log
from tqdm import tqdm

from . import util


# %%
def filter_signaldepths(
    path_signaldepths, path_sampledicts, max_median_coverage_factor, output_path
):
    sampledicts = util.load_sampledicts(path_sampledicts)
    csv.field_size_limit(sys.maxsize)
    # new csv reader
    reader = csv.reader(path_signaldepths.open(), delimiter="\t")
    # first pass: calculate statistics
    # create a dictionary per sampleID, per chromosme with entries: mean, median, Q1, Q3, std, min, max
    stats_dict = {int(sampleID): {} for sampleID in sampledicts[2].keys()}
    # create a data_dict that contains all data for all samples
    data_dict = {int(sampleID): {} for sampleID in sampledicts[2].keys()}
    # line in signaldepth file: 17      133177  133178  BNDR    39475   39476   109     a704b90c-00e0-4062-9be2-ab490892e83c    0       1       18      35
    # line in signaldepth file: chrom   pos     pos+1   type    start   end     length  readID                                 sampleID forward chr     depth
    for line in tqdm(reader):
        try:
            sampleID = int(line[8])
            chrom = line[10]
            depth = int(line[11])
            if chrom not in data_dict[sampleID]:
                data_dict[sampleID][chrom] = []
            data_dict[sampleID][chrom].append(depth)
        except:
            raise ValueError(f"Error in line: {line}")
    # now calc the stats from the data_dict
    for sampleID in data_dict.keys():
        for chrom in data_dict[sampleID].keys():
            data = data_dict[sampleID][chrom]
            stats_dict[sampleID][chrom] = {
                "mean": np.mean(data),
                "median": np.median(data),
                "Q1": np.quantile(data, 0.25),
                "Q3": np.quantile(data, 0.75),
                "std": np.std(data),
                "min": np.min(data),
                "max": np.max(data),
            }
    # second pass: filter out positions with depth > N*mean
    reader = csv.reader(path_signaldepths.open(), delimiter="\t")
    # new csv writer
    counter_filtered = 0
    counter_total = 0
    with open(output_path, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for line in tqdm(reader):
            sampleID = int(line[8])
            chrom = line[10]
            depth = int(line[11])
            counter_total += 1
            if (
                depth
                < max_median_coverage_factor * stats_dict[sampleID][chrom]["median"]
            ):
                writer.writerow(line)
            else:
                counter_filtered += 1
    if counter_total == 0:
        raise ("No signals found in signaldepth file.")
    log.info(
        f"Filtered {counter_filtered} signals from {counter_total} total signals. Ratio: {counter_filtered / counter_total}"
    )


def run(args, **kwargs):
    filter_signaldepths(
        path_signaldepths=args.path_signaldepths,
        path_sampledicts=args.path_sampledicts,
        max_median_coverage_factor=args.max_median_coverage_factor,
        output_path=args.output_path,
        **kwargs,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Filters a signaldepth file for positions with a depth higher than N*mean."
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="path_signaldepths",
        type=Path,
        required=True,
        help="Path to the signaldepth file.",
    )
    parser.add_argument(
        "-s",
        "--sampledicts",
        dest="path_sampledicts",
        type=Path,
        required=True,
        help="Path to the sampledicts.",
    )
    parser.add_argument(
        "-f",
        "--factor",
        dest="max_median_coverage_factor",
        type=float,
        required=True,
        help="Maximum median coverage factor.",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_path",
        type=Path,
        required=True,
        help="Path to the output file.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
