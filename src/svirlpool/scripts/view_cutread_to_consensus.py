import argparse
from pathlib import Path

from . import consensus_class, util


def view_cut_reads_from_consensus(
    consensusID: str,
    path_db: Path,
    terminal_width: int = 120,
    print_signals: bool = False,
):
    consensus: consensus_class.Consensus = next(
        util.yield_consensus_objects(consensusID=consensusID, path_database=path_db)
    )
    read_alignments = [
        recr.alignment.to_pysam() for recr in consensus.reconstructible_reads
    ]
    util.display_ascii_alignments(
        alignments=read_alignments, terminal_width=terminal_width
    )
    # then print the cut_read_alignment_signals sorted by position
    if print_signals:
        print("start pos", "size", "sv_type", sep="\t")
        flat_signals = sorted(
            [
                signal
                for ras in consensus.cut_read_alignment_signals
                for signal in ras.SV_signals
            ],
            key=lambda x: x.ref_start,
        )
        for signal in flat_signals:
            print(signal.ref_start, signal.size, signal.sv_type, sep="\t")


def run(args, **kwargs):
    view_cut_reads_from_consensus(
        consensusID=args.consensusID,
        path_db=args.database,
        terminal_width=args.terminal_width,
        print_signals=args.print_signals,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Get the cut reads for a given crID from the database and write them to a bam file."
    )
    parser.add_argument(
        "-c",
        "--consensusID",
        type=str,
        required=True,
        help="The crID to get the cut reads for.",
    )
    parser.add_argument(
        "-d", "--database", type=Path, required=True, help="The path to the database."
    )
    parser.add_argument(
        "-w",
        "--terminal-width",
        type=int,
        required=False,
        default=120,
        help="The width of the terminal to display the alignments. [120]",
    )
    parser.add_argument(
        "-s",
        "--print-signals",
        action="store_true",
        help="Print the signals of the cut reads. [False]",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
