# # this module is for debugging. It reads one or more files that were generated with the consensus module
# # read the scripts.consensus.py module and study the format of its output. These outputs are also parsed in consensus_containers_to_db
# # the input parameters are: one or more consensus result files, a reference genome, the output path of the alignments
# # procedure: the consensus result files are parsed and written to a tmp fasta file, then aligned with minimap2 to the reference genome  (using scripts.util.align_reads_with_minimap)
# # keep the style and CLI to the other modules

# # %%
# import argparse
# import json
# import sys
# import tempfile
# from pathlib import Path

# import cattrs
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from logzero import logfile
# from logzero import logger as log

# from ..localassembly import consensus_class

# # Add the parent directory to the path to allow imports
# sys.path.insert(0, str(Path(__file__).parent.parent))

# from ..scripts import util

# # %%


# def parse_crs_container_results(path: Path):
#     """
#     Parse CrsContainerResult objects from a consensus results file.
#     Supports both plain text and gzipped files.

#     Args:
#         path: Path to the consensus results file (plain or .gz)

#     Yields:
#         CrsContainerResult objects
#     """
#     # check if file is gzipped
#     if path.suffix == ".gz":
#         import gzip

#         open_file = gzip.open
#     else:
#         open_file = open

#     with open_file(path, "rt") as f:
#         for line in f:
#             crs_container_result = json.loads(line)
#             yield cattrs.structure(
#                 crs_container_result, consensus_class.CrsContainerResult
#             )


# def extract_consensus_sequences(
#     crs_container_results: list[consensus_class.CrsContainerResult],
# ) -> list[SeqRecord]:
#     """
#     Extract consensus sequences from CrsContainerResult objects.

#     Args:
#         crs_container_results: List of CrsContainerResult objects

#     Returns:
#         List of SeqRecord objects containing consensus sequences
#     """
#     consensus_sequences = []

#     for crs_result in crs_container_results:
#         for consensus_id, consensus_obj in crs_result.consensus_dicts.items():
#             if consensus_obj is None:
#                 log.warning(f"Skipping None consensus for ID {consensus_id}")
#                 continue

#             # Create a SeqRecord from the consensus object
#             # Note: consensus_sequence is already a string in the Consensus class
#             seq_record = SeqRecord(
#                 Seq(consensus_obj.consensus_sequence),
#                 id=consensus_id,
#                 name=consensus_id,
#                 description=f"Consensus sequence for {consensus_id}",
#             )
#             consensus_sequences.append(seq_record)

#     return consensus_sequences


# def align_consensus_results(
#     input_files: list[Path],
#     reference: Path,
#     output: Path,
#     threads: int = 4,
#     tech: str = "map-ont",
#     aln_args: str = "",
#     tmp_dir_path: Path | None = None,
# ) -> None:
#     """
#     Main function to align consensus results to a reference genome.

#     Args:
#         input_files: List of paths to consensus result files
#         reference: Path to reference genome (fasta)
#         output: Path to output BAM file
#         threads: Number of threads to use for alignment
#         tech: Minimap2 technology preset (default: map-ont)
#         aln_args: Additional alignment arguments for minimap2
#         tmp_dir_path: Optional temporary directory path
#     """
#     log.info(f"Processing {len(input_files)} consensus result file(s)")

#     # Parse all consensus results from all input files
#     all_crs_results = []
#     for input_file in input_files:
#         log.info(f"Reading consensus results from {input_file}")
#         crs_results = list(parse_crs_container_results(input_file))
#         all_crs_results.extend(crs_results)
#         log.info(f"  Loaded {len(crs_results)} CrsContainerResult objects")

#     # Extract consensus sequences
#     log.info("Extracting consensus sequences")
#     consensus_sequences = extract_consensus_sequences(all_crs_results)

#     if len(consensus_sequences) == 0:
#         log.warning("No consensus sequences found in input files. Nothing to align.")
#         return

#     log.info(f"Extracted {len(consensus_sequences)} consensus sequences")

#     # Create temporary directory for intermediate files
#     with tempfile.TemporaryDirectory(dir=tmp_dir_path) as tmp_dir:
#         # Write consensus sequences to temporary fasta file
#         tmp_fasta = Path(tmp_dir) / "consensus_sequences.fasta"
#         log.info(f"Writing consensus sequences to {tmp_fasta}")
#         with open(tmp_fasta, "w") as f:
#             SeqIO.write(consensus_sequences, f, "fasta")

#         # Align consensus sequences to reference genome
#         log.info(f"Aligning consensus sequences to reference genome {reference}")
#         util.align_reads_with_minimap(
#             reference=reference,
#             reads=tmp_fasta,
#             bamout=output,
#             tech=tech,
#             threads=threads,
#             aln_args=aln_args,
#         )

#     log.info(f"Alignment complete. Output written to {output}")
#     log.info(f"You can view the alignments with: samtools view {output}")


# # %%


# def run(args):
#     """Run the alignment from command line arguments."""
#     align_consensus_results(
#         input_files=args.input,
#         reference=args.reference,
#         output=args.output,
#         threads=args.threads,
#         tech=args.tech,
#         aln_args=args.aln_args,
#         tmp_dir_path=args.tmp_dir_path,
#     )


# def get_parser():
#     """Create argument parser."""
#     parser = argparse.ArgumentParser(
#         description="Align consensus sequences from consensus result files to a reference genome. "
#         "This is a debugging tool to visualize consensus sequences."
#     )
#     parser.add_argument(
#         "-i",
#         "--input",
#         type=Path,
#         required=True,
#         nargs="+",
#         help="Path(s) to consensus result file(s) (plain text or .gz). Can specify multiple files.",
#     )
#     parser.add_argument(
#         "-r",
#         "--reference",
#         type=Path,
#         required=True,
#         help="Path to reference genome (fasta format).",
#     )
#     parser.add_argument(
#         "-o",
#         "--output",
#         type=Path,
#         required=True,
#         help="Path to output BAM file with alignments.",
#     )
#     parser.add_argument(
#         "-t",
#         "--threads",
#         type=int,
#         default=4,
#         help="Number of threads to use for alignment (default: 4).",
#     )
#     parser.add_argument(
#         "--tech",
#         type=str,
#         default="map-ont",
#         help="Minimap2 technology preset (default: map-ont). Options: map-ont, map-pb, etc.",
#     )
#     parser.add_argument(
#         "--aln-args",
#         type=str,
#         default="",
#         help="Additional arguments to pass to minimap2 (default: empty string).",
#     )
#     parser.add_argument(
#         "--tmp-dir-path",
#         type=Path,
#         default=None,
#         help="Path to temporary directory for intermediate files (default: system temp dir).",
#     )
#     parser.add_argument(
#         "--logfile",
#         type=Path,
#         default=None,
#         help="Path to log file (default: log to console).",
#     )
#     return parser


# def main():
#     """Main entry point for the script."""
#     parser = get_parser()
#     args = parser.parse_args()

#     if args.logfile:
#         logfile(str(args.logfile))

#     run(args)


# if __name__ == "__main__":
#     main()
