"""Plot the copy-number tracks of copynumber_tracks.py from its database.

A separate step from generating the tracks, so that redrawing the plot never
rewrites the tracks (and with them everything downstream in the workflow).
"""

import argparse
import logging
from pathlib import Path

from .copynumber_tracks import load_copynumber_trees_from_db, plot_copynumber_tracks

log = logging.getLogger(__name__)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot the copy-number tracks from the database of copynumber_tracks."
    )
    parser.add_argument(
        "--db", type=Path, required=True, help="Copy-number tracks database (--output-db)"
    )
    parser.add_argument(
        "--reference-fai",
        type=Path,
        required=True,
        help="Reference FASTA index (.fai): chromosome order and sizes",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output figure")
    parser.add_argument("--title", default="Copy Number Tracks", help="Plot title")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    return parser


def main():
    args = get_parser().parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot_copynumber_tracks(
        cn_trees=load_copynumber_trees_from_db(args.db),
        output_figure=args.output,
        reference_fai=args.reference_fai,
        title=args.title,
    )


if __name__ == "__main__":
    main()
