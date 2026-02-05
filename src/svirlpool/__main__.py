#!/usr/bin/env python

import argparse
import os

from .scripts import cut_reads_from_alns, get_consensus_sequences, run_wf
from .svcalling import multisample_sv_calling
from .version import get_versions

#: The executables required for running svirlpool
# REQUIRED_EXECUTABLES = ("bedtools", "sort")


def get_parser():
    """Return argparse command line parser."""
    parser = argparse.ArgumentParser(
        description="SvirlPool generates SV calls and a database with all necessary data to perform high precision multi-samples SV calling.",
    )
    parser.add_argument(
        "--version", action="version", version=get_versions()["version"]
    )
    subparsers = parser.add_subparsers()

    # =========================================================================
    #  run_wf
    # =========================================================================

    parser_run_wf = subparsers.add_parser(
        "run", description="Run the analysis of a single sample."
    )

    parser_run_wf.add_argument(
        "--samplename",
        help="Name of the sample. Should be unique in your study.",
        required=True,
        type=str,
    )
    parser_run_wf.add_argument(
        "--workdir",
        help="base directory for the analysis, where all files are stored. Sub dirs are created in the process.",
        required=True,
        type=os.path.abspath,
    )
    parser_run_wf.add_argument(
        "--output",
        help="Output database filename (relative to workdir) or absolute path. Default: svirltile.db",
        required=False,
        type=str,
        default="svirltile.db",
    )
    parser_run_wf.add_argument(
        "--alignments",
        help="Path to the bam file",
        required=True,
        type=os.path.abspath,
    )
    # parser_run_wf.add_argument(
    #     "--coverage", help="average depth of read coverage on the genome. Can be computed e.g. with mosdepth", required=True, type=float,
    # )
    parser_run_wf.add_argument(
        "--reference",
        help="reference genome in fasta format",
        required=True,
        type=os.path.abspath,
    )
    parser_run_wf.add_argument(
        "--regions",
        help="regions of interest in bed format. If not provided, will be generated from reference fasta index.",
        required=False,
        type=os.path.abspath,
        default=None,
    )
    parser_run_wf.add_argument(
        "--trf",
        help="tandem repeat regions in bed format",
        required=True,
        type=os.path.abspath,
    )
    parser_run_wf.add_argument(
        "--mononucleotides",
        help="mononucleotide regions in bed format",
        required=True,
        type=os.path.abspath,
    )
    # parser_run_wf.add_argument(
    #     "--unique-regions",
    #     help="unique regions in bed format",
    #     required=True,
    #     type=os.path.abspath,
    # )
    parser_run_wf.add_argument(
        "--lamassemble-mat",
        help="lamassamble matrix file used for the final consensus assembly.",
        required=True,
        type=os.path.abspath,
    )
    parser_run_wf.add_argument(
        "--threads", help="number of threads to use", required=True, type=int
    )
    parser_run_wf.add_argument(
        "--cores-per-consensus",
        help="number of cores per consensus clustering and assembly [1]",
        required=False,
        type=int,
        default=1,
    )
    parser_run_wf.add_argument(
        "--max-coverage-per-region",
        type=int,
        default=400,
        help="Maximum coverage threshold. Alignments in regions exceeding this coverage will be skipped. Default is 400.",
    )
    parser_run_wf.add_argument(
        "--N-files-per-dir",
        help="number of files per directory",
        required=False,
        type=int,
        default=200,
    )
    parser_run_wf.add_argument(
        "--min-cr-size",
        help="minimum candidate number region size",
        required=False,
        type=int,
        default=500,
    )
    parser_run_wf.add_argument(
        "--cr-merge-buffer",
        help="candidate region merge buffer",
        required=False,
        type=int,
        default=1200,
    )
    parser_run_wf.add_argument(
        "--filter-absolute",
        help="filter absolute in signal collection. decrease with less noisy data",
        required=False,
        type=float,
        default=1.0,
    )
    parser_run_wf.add_argument(
        "--filter-normalized",
        help="filter normalized in signal collection. decrease with less noisy data",
        required=False,
        type=float,
        default=0.06,
    )
    parser_run_wf.add_argument(
        "--min-mapq",
        help="minimum mapping quality",
        required=False,
        type=int,
        default=1,
    )
    parser_run_wf.add_argument(
        "--min-signal-size",
        help="minimum signal size",
        required=False,
        type=int,
        default=12,
    )
    parser_run_wf.add_argument(
        "--cn-dispersion",
        help="Copy number dispersion parameter for Negative Binomial emission model (default: 0.1). "
        "Higher values = more tolerance for coverage variance. Range: 0.05-0.2. "
        "Use ~0.05-0.1 for clean data, ~0.15-0.2 for noisy data.",
        required=False,
        type=float,
        default=0.1,
    )
    parser_run_wf.add_argument(
        "--min-alt-reads",
        help="minimum alternative reads",
        required=False,
        type=int,
        default=3,
    )
    parser_run_wf.add_argument(
        "--min-sv-size",
        help="minimum structural variant size",
        required=False,
        type=int,
        default=30,
    )
    parser_run_wf.add_argument(
        "--snakemake-unlock",
        help="unlock snakemake",
        required=False,
        action="store_true",
    )
    # add executor slurm
    parser_run_wf.add_argument(
        "--executor-slurm-jobs",
        help="use slurm executor is used with N jobs",
        required=False,
        type=int,
        default=0,
    )
    parser_run_wf.add_argument(
        "--dont-merge-horizontally",
        help="do not merge horizontally, i.e. do not merge SVs that are close to each other or in tandem repeats",
        required=False,
        action="store_true",
    )

    parser_run_wf.set_defaults(fast=False, func=run_wf.run_wf)

    # =========================================================================
    #  multisample sv calling
    # =========================================================================

    parser_run_sv_calling = subparsers.add_parser(
        "sv-calling",
        description="Call SVs from the final databases (svirltiles) and write them to a vcf file.",
    )
    multisample_sv_calling.add_arguments(parser_run_sv_calling)

    parser_run_sv_calling.set_defaults(fast=False, func=multisample_sv_calling.run)

    # =========================================================================
    #  read cutting
    # =========================================================================

    parser_cut_reads = subparsers.add_parser(
        "cut-reads", description="Cut reads from alignments given a region of interest."
    )
    cut_reads_from_alns.add_arguments(parser_cut_reads)
    parser_cut_reads.set_defaults(fast=False, func=cut_reads_from_alns.run)

    # =========================================================================
    #  get consensus sequences
    # =========================================================================

    parser_get_consensus = subparsers.add_parser(
        "get-consensus",
        description="Extract consensus sequences from a SVIRLPOOL database to FASTA file.",
    )
    get_consensus_sequences.add_arguments(parser_get_consensus)
    parser_get_consensus.set_defaults(fast=False, func=get_consensus_sequences.run)

    return parser


# def execs_available(execs: list[str]) -> None:
#     import subprocess
#     ok = True
#     for exec in execs:
#         res = subprocess.run(["which", exec], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         if res.returncode != 0:
#             logger.warn("Could not find executable %s via `which %s`", exec, exec)
#             ok = False
#     return ok


def main():
    # if not execs_available(REQUIRED_EXECUTABLES):
    #     logger.error("Missing some executables. The program will eventually fail.")
    # else:
    #     logger.debug("External executables present: %s", ", ".join(REQUIRED_EXECUTABLES))
    parser = get_parser()
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
    else:
        args.func(args)


if __name__ == "__main__":
    main()
