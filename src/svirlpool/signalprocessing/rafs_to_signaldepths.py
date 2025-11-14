# this script extracts the effective signals from a given RAF file and adds depth of coverage to the signals.
# a signal line is of the form: [chrID,ref_start,ref_end,type,read_start,read_end,size,readname,sampleID,forward,chr,deph_of_coverage]
# where depth of coverage counts all alignments of the corresponding sample in the effective interval.

import argparse
import multiprocessing as mp
import subprocess
import tempfile
from pathlib import Path
from shlex import split

import logzero

from ..signalprocessing import alignments_to_rafs
from ..util import datatypes, util

# 0     1         2       3    4          5        6    7        8          9       10
# [chrID,ref_start,ref_end,type,read_start,read_end,size,readname,samplename,forward,chr] # +deph_of_coverage (later)


def get_depths(rafs: list[datatypes.ReadAlignmentFragment], output: Path) -> None:
    """Get the maximum depth of coverage over the given intervals, excluding the specified readname."""
    with tempfile.TemporaryDirectory() as tdir:
        tmp_bed_intervals = Path(tdir) / "intervals.bed"
        # write intervals to bed file for depth computation
        with open(tmp_bed_intervals, "w") as f:
            for x in sorted(rafs, key=lambda x: x.reference_alignment_start):
                print(
                    x.referenceID,
                    x.reference_alignment_start,
                    x.reference_alignment_end,
                    sep="\t",
                    file=f,
                )

        # now write the signal lines to a tmp file
        tmp_bed_signals = Path(tdir) / "signals.bed"
        with open(tmp_bed_signals, "w") as f:
            for raf in rafs:
                x: datatypes.ReadAlignmentFragment = raf
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
                    print(*line, sep="\t", file=f)

        # sort the signal bed file
        sorted_signals = Path(tdir) / "signals_sorted.bed"
        cmd_sort = f"sort -k1,1 -k2,2n {tmp_bed_signals}"
        with open(sorted_signals, "w") as out_f:
            subprocess.check_call(split(cmd_sort), stdout=out_f)

        # now use bedtools coverage to get the coverage of each signal line by the intervals in the tmp bed file
        cmd = f"bedtools coverage -a {sorted_signals} -b {tmp_bed_intervals} -counts"
        with open(output, "w") as out_f:
            subprocess.check_call(split(cmd), stdout=out_f)


def _process_chromosome(args_tuple):
    """Helper function for parallel processing of chromosomes."""
    chrom_rafs_serialized, temp_output = args_tuple
    # Deserialize RAFs
    chrom_rafs = [
        datatypes.ReadAlignmentFragment.from_unstructured(raf_data)
        for raf_data in chrom_rafs_serialized
    ]
    get_depths(chrom_rafs, temp_output)
    return temp_output


def signaldepths_from_rafs(input: Path, output: Path, threads: int = 8) -> None:
    """creates a signaldepths file from a RAF file. Coverage is only counted for the unique region spanning \
intervals of all alignments (RAFs). The output file is a tsv with olumns: \
[chrID,ref_start,ref_end,type,read_start,read_end,size,readname,sampleID,forward,chr,deph_of_coverage]"""

    rafs: dict[str : list[datatypes.ReadAlignmentFragment]] = {}  # chr:list of rafs
    for raf in util.yield_from_raf(input):
        if raf.reference_name not in rafs:
            rafs[raf.reference_name] = []
        rafs[raf.reference_name].append(raf)

    # procedure:
    # create jobs for multiprocessing. each job processes the ReadalignmentFragments of one chromosome and writes to a temporary file.
    # create a temporary directory with target files for each process
    # execute the jobs in parallel
    # combine the temporary files into the output file
    # sort the output file by chrID and ref_start

    with tempfile.TemporaryDirectory() as tdir:
        temp_dir = Path(tdir)

        # Prepare jobs for parallel processing - serialize RAFs
        jobs = []
        for chrom, chrom_rafs in rafs.items():
            if not chrom_rafs:
                continue
            temp_output = temp_dir / f"{chrom}_signals.bed"
            # Serialize RAFs for multiprocessing
            chrom_rafs_serialized = [raf.unstructure() for raf in chrom_rafs]
            jobs.append((chrom_rafs_serialized, temp_output))

        # Process chromosomes in parallel
        with mp.Pool(processes=min(threads, len(jobs))) as pool:
            temp_files = pool.map(_process_chromosome, jobs)

        # Combine all temporary files into the output file
        with open(output, "w") as out_f:
            for temp_file in temp_files:
                if temp_file.exists():
                    with open(temp_file, "r") as in_f:
                        out_f.write(in_f.read())


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
    loglevel = logzero.INFO
    if args.debug:
        loglevel = logzero.DEBUG
    if args.silent:
        loglevel = logzero.ERROR
    if args.log_file:
        logzero.logfile(Path(args.log_file), loglevel=loglevel)


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
        help="Set logzero logging level to logzero.DEBUG.",
    )
    parser.add_argument(
        "--silent",
        action="store_true",
        help="Set logzero logging level to logzero.ERROR.",
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
