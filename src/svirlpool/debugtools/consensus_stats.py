#!/usr/bin/env python3
"""Generate statistics about consensus sequences from FASTA files.

Reads consensus FASTA files where core sequences are in uppercase letters
and flanking padding is in lowercase letters. For each sequence, extracts
the left padding length, core length, and right padding length, then
computes min/median/max statistics across all sequences per file.
"""

import argparse
import logging
import re
import statistics
import sys
from pathlib import Path
from typing import Generator

from tqdm import tqdm

log = logging.getLogger(__name__)


def _parse_consensus_sequence(sequence: str) -> dict[str, int]:
    """Parse a single consensus sequence into padding and core lengths.

    Core consensus sequences are uppercase; flanking padding is lowercase.
    Pattern: lowercase (left pad) + uppercase (core) + lowercase (right pad).

    Returns a dict with keys: left_padding, core, right_padding.
    """
    # Split into alternating runs of upper/lower case
    # Each block is a maximal run of either upper or lower case letters
    blocks = re.findall(r"[A-Z]+|[a-z]+", sequence)

    left_padding = 0
    core = 0
    right_padding = 0

    if not blocks:
        return {"left_padding": 0, "core": 0, "right_padding": 0}

    # Determine the pattern based on the first block
    first_is_upper = blocks[0].isupper()

    if first_is_upper:
        # Pattern: UPPER (core) or UPPER (core) + lower (right pad) + ...
        if len(blocks) == 1:
            core = len(blocks[0])
        elif len(blocks) == 2:
            core = len(blocks[0])
            right_padding = len(blocks[1])
        else:
            # Multiple blocks: first upper = core, last lower = right pad,
            # middle blocks summed into core
            core = len(blocks[0])
            for b in blocks[1:-1]:
                core += len(b)
            right_padding = len(blocks[-1])
    else:
        # Pattern: lower (left pad) + upper (core) + optional lower (right pad)
        if len(blocks) == 1:
            # All lowercase — no core
            left_padding = len(blocks[0])
        elif len(blocks) == 2:
            left_padding = len(blocks[0])
            core = len(blocks[1])
        else:
            # Multiple blocks: first lower = left pad, last lower = right pad,
            # middle blocks summed into core
            left_padding = len(blocks[0])
            core = len(blocks[1])
            for b in blocks[2:-1]:
                core += len(b)
            right_padding = len(blocks[-1])

    return {"left_padding": left_padding, "core": core, "right_padding": right_padding}


def _count_fasta(fasta_path: Path) -> int:
    """Quick count of sequences in a FASTA file."""
    count = 0
    with open(fasta_path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def _parse_fasta(fasta_path: Path) -> Generator[dict[str, int], None, None]:
    """Yield consensus sequence stats one at a time from a FASTA file.

    Yields dicts with keys: left_padding, core, right_padding.
    """
    current_seq_parts: list[str] = []

    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # Flush previous sequence
                if current_seq_parts:
                    seq = "".join(current_seq_parts)
                    yield _parse_consensus_sequence(seq)
                    current_seq_parts = []
            elif line.strip():
                current_seq_parts.append(line)

        # Flush last sequence
        if current_seq_parts:
            seq = "".join(current_seq_parts)
            yield _parse_consensus_sequence(seq)


def _median_or_zero(values: list[int]) -> int:
    """Return the median of *values*, or 0 if the list is empty."""
    return int(statistics.median(values)) if values else 0


def _min_or_zero(values: list[int]) -> int:
    return min(values) if values else 0


def _max_or_zero(values: list[int]) -> int:
    return max(values) if values else 0


def compute_stats(fasta_path: Path) -> dict[str, int | str]:
    """Compute consensus statistics for a single FASTA file.

    Returns a dict with keys:
        file, core_min, core_median, core_max,
        pad_left_min, pad_left_median, pad_left_max,
        pad_right_min, pad_right_median, pad_right_max, count
    """
    core_lengths: list[int] = []
    pad_left_sizes: list[int] = []
    pad_right_sizes: list[int] = []

    total = _count_fasta(fasta_path)

    for seq_stats in tqdm(_parse_fasta(fasta_path), total=total, desc=str(fasta_path)):
        core_lengths.append(seq_stats["core"])
        pad_left_sizes.append(seq_stats["left_padding"])
        pad_right_sizes.append(seq_stats["right_padding"])

    count = len(core_lengths)

    return {
        "file": str(fasta_path),
        "core_min": _min_or_zero(core_lengths),
        "core_median": _median_or_zero(core_lengths),
        "core_max": _max_or_zero(core_lengths),
        "pad_left_min": _min_or_zero(pad_left_sizes),
        "pad_left_median": _median_or_zero(pad_left_sizes),
        "pad_left_max": _max_or_zero(pad_left_sizes),
        "pad_right_min": _min_or_zero(pad_right_sizes),
        "pad_right_median": _median_or_zero(pad_right_sizes),
        "pad_right_max": _max_or_zero(pad_right_sizes),
        "count": count,
    }


def consensus_stats(input_fasta: list[Path], output: Path | None = None) -> None:
    """Print consensus statistics as TSV to *output* (or stdout)."""
    headers = [
        "file",
        "core_min",
        "core_median",
        "core_max",
        "pad_left_min",
        "pad_left_median",
        "pad_left_max",
        "pad_right_min",
        "pad_right_median",
        "pad_right_max",
        "count",
    ]

    if output is not None:
        out = open(output, "w", encoding="utf-8")
    else:
        out = sys.stdout

    try:
        print("\t".join(headers), file=out)

        for fasta_path in input_fasta:
            stats = compute_stats(fasta_path)
            values = [str(stats[h]) for h in headers]
            print("\t".join(values), file=out)
            log.info("Processed %s: %d consensus sequences", fasta_path, stats["count"])
    finally:
        if out is not sys.stdout:
            out.close()


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate statistics about consensus sequences from FASTA files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        nargs="+",
        help="One or more consensus FASTA files.",
        metavar="FASTA",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output TSV file (default: stdout).",
        metavar="TSV",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    else:
        logging.basicConfig(level=logging.WARNING, format="%(name)s: %(message)s")

    consensus_stats(args.input, args.output)


if __name__ == "__main__":
    main()
