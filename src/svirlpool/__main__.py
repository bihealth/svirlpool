#!/usr/bin/env python

import argparse
import os

from .scripts import (
    cut_reads_from_alns,
    get_consensus_sequences,
    multisample_sv_calling,
    run_wf,
)
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
        "--unique-regions",
        help="unique regions in bed format",
        required=True,
        type=os.path.abspath,
    )
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
    # parser_run_wf.add_argument(
    #     "--consensus-mode",
    #     help="consensus mode [greedy, ava, rapidfuzz] determines the method of consensus clustering and assembly",
    #     required=False,
    #     type=str,
    #     default="greedy",
    # )
    parser_run_wf.add_argument(
        "--N-files-per-dir",
        help="number of files per directory",
        required=False,
        type=int,
        default=200,
    )
    parser_run_wf.add_argument(
        "--optimal-cr-size",
        help="optimal candidate region size",
        required=False,
        type=int,
        default=600,
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
    parser_run_sv_calling.add_argument(
        "--input",
        help="Paths to the per sample svirltiles.",
        nargs="+",
        required=True,
        type=os.path.abspath,
    )
    parser_run_sv_calling.add_argument(
        "--reference",
        help="Path to the reference genome.",
        required=True,
        type=os.path.abspath,
    )
    parser_run_sv_calling.add_argument(
        "--output",
        help="Path to the output vcf file. An additional gzipped copy is created automatically",
        required=True,
        type=os.path.abspath,
    )
    parser_run_sv_calling.add_argument(
        "--sv-types",
        type=str,
        nargs="+",
        required=False,
        default=["DEL", "INS"],
        help=f"List of structural variant types to call. Default: DEL INS. Allowed: {str(multisample_sv_calling.SUPPORTED_SV_TYPES)}.",
    )
    parser_run_sv_calling.add_argument(
        "--min-sv-size",
        type=int,
        required=False,
        default=30,
        help="Minimum size of structural variants to call. Default: 30.",
    )
    parser_run_sv_calling.add_argument(
        "--threads",
        help="Number of threads to use.",
        required=False,
        default=8,
        type=int,
    )
    parser_run_sv_calling.add_argument(
        "--max-cohens-d",
        help="Maximum Cohen's d value for merging SVs (default: 2.0).",
        required=False,
        type=float,
        default=2.0,
    )
    parser_run_sv_calling.add_argument(
        "--min-kmer-overlap",
        help="Minimum k-mer overlap for merging SVs (default: 0.7).",
        required=False,
        type=float,
        default=0.7,
    )
    parser_run_sv_calling.add_argument(
        "--near",
        help="Maximum distance for merging SVs (default: 500).",
        required=False,
        type=int,
        default=500,
    )
    parser_run_sv_calling.add_argument(
        "-n",
        "--apriori-size-difference-fraction-tolerance",
        help="Size difference fraction tolerance for merging SVs (default: 0.1). Decrease for stronger separation of haplotypes",
        type=float,
        default=0.1,
    )

    parser_run_sv_calling.add_argument(
        "--find-leftmost-reference-position",
        help="When determining the reference position of an SVcomposite, use the leftmost position of all underlying SVpatterns instead of the one with most supporting reads * size. This might be better aligned with the giab SV benchmark, but should be discussed in the paper.",
        action="store_true",
        default=False,
    )

    parser_run_sv_calling.add_argument(
        "--candidate-regions-file",
        help="Optional TSV file with samplename and comma-separated candidate region IDs (crID) to filter SVpatterns. Format: samplename<TAB>crID1,crID2,crID3",
        type=os.path.abspath,
        default=None,
    )
    parser_run_sv_calling.add_argument(
        "--symbolic-threshold",
        help="Sequence length threshold for using symbolic alleles in VCF (default: 100000). Sequences longer than this will be written to a companion FASTA file.",
        type=int,
        default=100000,
    )
    parser_run_sv_calling.add_argument(
        "--tmp-dir-path",
        help="Path to temporary directory (default: system temp dir).",
        type=os.path.abspath,
        default=None,
    )
    parser_run_sv_calling.add_argument(
        "--skip-covtrees",
        help="Skip computing coverage trees from sample data. Instead, use uniform coverage of 30 across all chromosomes. This significantly speeds up execution when genotype coverage information is not critical.",
        action="store_true",
        default=False,
    )
    parser_run_sv_calling.add_argument(
        "--verbose", help="Enable verbose output.", action="store_true", default=False
    )

    parser_run_sv_calling.set_defaults(fast=False, func=multisample_sv_calling.run)

    # =========================================================================
    #  read cutting
    # =========================================================================

    parser_cut_reads = subparsers.add_parser(
        "cut-reads", description="Cut reads from alignments given a region of interest."
    )
    parser_cut_reads.add_argument(
        "--input",
        help="Path to the alignments file (bam/sam).",
        required=True,
        type=os.path.abspath,
    )
    parser_cut_reads.add_argument(
        "--region",
        help="Region to cut reads from. Format: chr:start-end",
        required=True,
        type=str,
    )
    parser_cut_reads.add_argument(
        "--output",
        help="Path to the output fasta/fastq (.gz) file. If the filename extension is .gz, the output file will be written as a gzip compressed file.",
        required=True,
        type=os.path.abspath,
    )
    parser_cut_reads.add_argument(
        "--buffer-clipped-length",
        help="The maximum number of bases that are included in the cut sequences, if they have been hard or soft clipped within the region bounds. Defaults to 1000.",
        required=False,
        type=int,
        default=1000,
    )
    parser_cut_reads.set_defaults(fast=False, func=cut_reads_from_alns.run)

    # =========================================================================
    #  get consensus sequences
    # =========================================================================

    parser_get_consensus = subparsers.add_parser(
        "get-consensus",
        description="Extract consensus sequences from a SVIRLPOOL database to FASTA file.",
    )
    parser_get_consensus.add_argument(
        "-i",
        "--input",
        help="Path to the SVIRLPOOL database file.",
        required=True,
        type=os.path.abspath,
    )
    parser_get_consensus.add_argument(
        "-o",
        "--output",
        help="Path to the output FASTA file for consensus sequences. Use .gz extension for gzip compression (will be indexed with samtools faidx).",
        required=True,
        type=os.path.abspath,
    )
    parser_get_consensus.add_argument(
        "--batch-size",
        help="Number of sequences to batch before writing to file (default: 1000).",
        required=False,
        type=int,
        default=1000,
    )
    parser_get_consensus.add_argument(
        "--consensus-ids",
        help="Optional list of consensus IDs to extract. If not provided, all sequences will be extracted.",
        required=False,
        type=str,
        nargs="+",
        default=None,
    )
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
