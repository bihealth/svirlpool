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
    parser_run_wf.add_argument(
        "--lamassemble-mat",
        help="lamassamble matrix file used for the final consensus assembly. Required when --consensus-method is 'lamassemble'.",
        required=False,
        default=None,
        type=os.path.abspath,
    )
    parser_run_wf.add_argument(
        "--consensus-method",
        help="Method for consensus assembly: 'lamassemble' (default) or 'racon'.",
        required=False,
        type=str,
        choices=["lamassemble", "racon"],
        default="lamassemble",
    )
    parser_run_wf.add_argument(
        "--threads", help="number of threads to use", required=True, type=int
    )
    parser_run_wf.add_argument(
        "--consensus-escalation",
        help="THREADS:SECONDS levels at which a candidate-region container is processed "
        "in the consensus stage. Containers are processed at the first level; one in "
        "which the all-vs-all alignment or the assembly timed out is processed again at "
        "the next level, up to the last one (the hard ceiling; then its result is kept "
        "degraded). Threads are capped at --threads. Add levels on large machines, e.g. "
        "1:20,4:60,12:120,32:300 [1:20,4:60,12:120]",
        required=False,
        type=str,
        default="1:20,4:60,12:120",
    )
    parser_run_wf.add_argument(
        "--consensus-max-mem-mb",
        help="memory ceiling (MB) of a consensus batch job. A failed job is retried with "
        "double the memory, from 2048 MB up to this [16384]",
        required=False,
        type=int,
        default=16384,
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
        "--cutoff-median-readcount-per-region",
        help="cutoff multiplier for median read count per candidate region. regions with excessive read counts are filtered out",
        required=False,
        type=float,
        default=6.0,
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
        "--log-level",
        help="Set the logging level for all workflow modules that support it.",
        required=False,
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
    )
    parser_run_wf.add_argument(
        "--single-evidence-gt",
        help="Pass this flag through the workflow configuration for experimental single-read evidence genotyping.",
        required=False,
        action="store_true",
        default=False,
    )
    parser_run_wf.add_argument(
        "--max-padding-size",
        help="Maximum number of bases to use for padding flanks (default: 100000).",
        required=False,
        type=int,
        default=100000,
    )
    parser_run_wf.add_argument(
        "--max-consensus-copy-number",
        help="Maximum estimated copy number of a candidate-region container for which a "
        "consensus is still attempted (default: 4). Containers exceeding this threshold are "
        "skipped as too complex and produce no consensus (and therefore no SV calls) for that "
        "region. Raise (e.g. 6-8) to recover SVs in higher-copy / complex tandem-repeat regions "
        "at the cost of runtime and potential noise.",
        required=False,
        type=int,
        default=4,
    )
    parser_run_wf.add_argument(
        "--consensus-clustering-mode",
        help="Experimental: how the reads of a candidate-region container are split into "
        "alleles before assembly. 'legacy' (default): KMeans on summed indels, else spectral "
        "clustering with k = local copy number. 'phased': phase the reads by the SNVs and SVs "
        "in their all-vs-all alignments and build one consensus per allele found.",
        required=False,
        choices=("legacy", "phased"),
        default="legacy",
    )
    parser_run_wf.add_argument(
        "--phasing-flank",
        help="Experimental: flank (bp) around the candidate regions to which reads are cut "
        "for read phasing with --consensus-clustering-mode phased (default: 10000).",
        required=False,
        type=int,
        default=10000,
    )
    parser_run_wf.add_argument(
        "--phasing-fallback",
        help="Experimental: with --consensus-clustering-mode phased, what to do when the "
        "phasing finds fewer than two alleles: 'single' (default) one consensus from all "
        "reads, 'legacy' the legacy clustering.",
        required=False,
        choices=("single", "legacy"),
        default="single",
    )
    parser_run_wf.add_argument(
        "--rerun-triggers",
        help="Snakemake rerun triggers. Comma-separated list of triggers that cause a rule to be rerun. "
        "Allowed values: mtime, params, input, software-env, code. "
        "Default is the Snakemake default (all triggers). "
        "Use 'mtime' to skip reruns caused by code or parameter changes (useful when resuming after source-code edits).",
        required=False,
        type=str,
        default=None,
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
