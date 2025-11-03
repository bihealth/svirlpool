# analysis script to investigate running times in log files
# grep all log files in the directory
# %%
import argparse
import datetime
import glob
import logging as log
import re
from pathlib import Path

import pandas as pd
from tqdm import tqdm

COLNAMES = ["load", "run", "save", "total", "num_clusters"]


def parse_file(logpath: str | Path) -> list[int]:
    try:
        lines = [line for line in open(logpath).readlines() if line.startswith("[I")]
        # get the start time
        time_start = datetime.datetime.strptime(lines[0].split(" ")[2], "%H:%M:%S")
        # find line with "loaded data"
        line_loaded = next(
            (line for line in lines if line.endswith("loaded data\n")), None
        )
        time_loaded = datetime.datetime.strptime(line_loaded.split(" ")[2], "%H:%M:%S")
        # num clusters = iterations, so find the line with "cluster:" searching from the end
        line_cluster = next(
            (line for line in reversed(lines) if re.search(r"cluster:\d+", line)), None
        )
        # get the number of clusters
        num_clusters = int(re.search(r"cluster:(\d+)", line_cluster).group(1))
        # find the line with ".txt" searching from the end
        line_save = next(
            (line for line in reversed(lines) if line.endswith(".txt.\n")), None
        )
        time_save = datetime.datetime.strptime(line_save.split(" ")[2], "%H:%M:%S")
        # get the time fromt he last line
        time_end = datetime.datetime.strptime(lines[-1].split(" ")[2], "%H:%M:%S")

        delta_load = time_loaded - time_start
        delta_run = time_save - time_loaded
        delta_save = time_end - time_save
        delta_total = time_end - time_start
    except:
        return [-1, -1, -1, -1, -1]

    return [
        delta_load.seconds,
        delta_run.seconds,
        delta_save.seconds,
        delta_total.seconds,
        num_clusters,
    ]


# %%


def run(args):
    base_path = Path(args.input_dir)
    log.info("globbing the log files..")
    logpaths = glob.glob(str(base_path / "*/*.log"))
    log.info("found %d log files", len(logpaths))
    log.info("parsing the log files..")
    df: pd.DataFrame = pd.DataFrame(
        [parse_file(fp) for fp in tqdm(logpaths)], columns=COLNAMES
    )
    # save the df to output
    df.to_csv(args.output, sep="\t", index=False)
    log.info("Done. Saved the df to %s", args.output)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Parses all times from the consensus logs."
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        type=str,
        required=True,
        help="Path to the consensus directory of a processed sample.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to the output tsv file that contains the dataframe data.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
