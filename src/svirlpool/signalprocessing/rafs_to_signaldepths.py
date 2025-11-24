# this script extracts the effective signals from a given RAF file and adds depth of coverage to the signals.
# a signal line is of the form: [chrID,ref_start,ref_end,type,read_start,read_end,size,readname,sampleID,forward,chr,deph_of_coverage]
# where depth of coverage counts all alignments of the corresponding sample in the effective interval.

import argparse
import logging
import subprocess
import tempfile
from collections.abc import Iterator
from pathlib import Path
from shlex import split

from tqdm import tqdm

from ..signalprocessing import alignments_to_rafs
from ..util import datatypes, util

log = logging.getLogger(__name__)

# 0     1         2       3    4          5        6    7        8          9       10
# [chrID,ref_start,ref_end,type,read_start,read_end,size,readname,samplename,forward,chr] # +deph_of_coverage (later)


def get_depths(rafs: Iterator[datatypes.ReadAlignmentFragment], output: Path) -> None:
    """Get the maximum depth of coverage over the given intervals, excluding the specified readname.

    Processes RAFs in a streaming fashion to avoid loading all data into memory.
    """
    with tempfile.TemporaryDirectory() as tdir:
        tmp_bed_intervals = Path(tdir) / "intervals.bed"
        tmp_bed_signals = Path(tdir) / "signals.bed"

        log.debug(f"Writing temporary intervals file: {tmp_bed_intervals}")
        # Stream RAFs and write both intervals and signals in one pass
        with (
            open(tmp_bed_intervals, "w") as intervals_f,
            open(tmp_bed_signals, "w") as signals_f,
        ):
            for raf in rafs:
                x: datatypes.ReadAlignmentFragment = raf
                # Write interval for this RAF
                print(
                    x.referenceID,
                    x.reference_alignment_start,
                    x.reference_alignment_end,
                    sep="\t",
                    file=intervals_f,
                )

                # Write all signal lines for this RAF
                for signal in x.SV_signals:
                    y: datatypes.SVsignal = signal
                    line = [
                        x.referenceID,
                        y.ref_start,
                        y.ref_end,
                        y.sv_type,
                        y.read_start,
                        y.read_end,
                        y.size,
                        raf.read_name,
                        x.samplename,
                        int(x.alignment_forward),
                        x.reference_name,
                    ]
                    print(*line, sep="\t", file=signals_f)

        # Sort the intervals file by position
        log.debug(f"Sorting temporary intervals file: {tmp_bed_intervals}")
        sorted_intervals = Path(tdir) / "intervals_sorted.bed"
        cmd_sort_intervals = f"sort -k1,1 -k2,2n {tmp_bed_intervals}"
        with open(sorted_intervals, "w") as out_f:
            subprocess.check_call(split(cmd_sort_intervals), stdout=out_f)

        # Sort the signal bed file
        log.debug(f"Sorting temporary signals file: {tmp_bed_signals}")
        sorted_signals = Path(tdir) / "signals_sorted.bed"
        cmd_sort = f"sort -k1,1 -k2,2n {tmp_bed_signals}"
        with open(sorted_signals, "w") as out_f:
            subprocess.check_call(split(cmd_sort), stdout=out_f)

        # now use bedtools coverage to get the coverage of each signal line by the intervals in the tmp bed file
        log.debug(
            f"Running bedtools coverage with -a {sorted_signals} and -b {sorted_intervals}"
        )
        cmd = f"bedtools coverage -a {sorted_signals} -b {sorted_intervals} -counts"
        with open(output, "w") as out_f:
            subprocess.check_call(split(cmd), stdout=out_f)


def signaldepths_from_rafs(input: Path, output: Path, threads: int = 8) -> None:
    """Creates a signaldepths file from a RAF file.

    Coverage is only counted for the unique region spanning intervals of all alignments (RAFs).
    The output file is a tsv with columns:
    [chrID,ref_start,ref_end,type,read_start,read_end,size,readname,sampleID,forward,chr,depth_of_coverage]

    Streams all RAFs from the input file without loading them into memory.
    """
    # Stream all RAFs and process them directly
    get_depths(tqdm(util.yield_from_raf(input)), output)


def create_combined_files(
    input: Path, output: Path, threads: int, tmp_dir_path: Path | None = None
) -> None:
    """creates effective intervals and effective signals files from a RAF file."""
    with tempfile.TemporaryDirectory(dir=tmp_dir_path) as tdir:
        tmp_output_txt = tempfile.NamedTemporaryFile(
            dir=tdir,
            delete=False if tmp_dir_path else True,
            prefix="output_effective_intervals",
            suffix=".bed",
        )
        signaldepths_from_rafs(
            input=input, output=Path(tmp_output_txt.name), threads=threads
        )
        alignments_to_rafs.compress_and_index_bedlike(
            input=Path(tmp_output_txt.name),
            output=Path(output),
            threads=threads,
            sort_numerically=True,
        )


def setup_logging(args):
    loglevel = logging.INFO
    if args.debug:
        loglevel = logging.DEBUG
    if args.silent:
        loglevel = logging.ERROR
    # basic config to ensure handlers/format is set
    logging.basicConfig(
        level=loglevel, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s"
    )
    if args.log_file:
        file_handler = logging.FileHandler(Path(args.log_file))
        file_handler.setLevel(loglevel)
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s")
        )
        logging.getLogger().addHandler(file_handler)


def run(args, **kwargs):
    setup_logging(args)
    create_combined_files(
        input=args.input,
        output=args.output,
        threads=args.threads,
        tmp_dir_path=args.tmp,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Generate a tsv of signals from a RAFs file. The output contains as the last column the effective depth of coverage."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the RAF file (cols: chr,start,end,raf).",
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Path to the output file."
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=8,
        help="Number of threads to use. Default is 8.",
    )
    parser.add_argument(
        "--tmp",
        default=None,
        help="Path to a temporary directory to store intermediate results. Default uses the system temporary directory.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Set logging level to DEBUG.",
    )
    parser.add_argument(
        "--silent",
        action="store_true",
        help="Set logging level to ERROR.",
    )
    parser.add_argument(
        "--log-file", required=False, default=None, help="Path to an optional log file."
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()


# %%
