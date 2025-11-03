# this script analysis input alignment files
# by counting all indel signals (extracted from the cigar strings)
# and recording alignment lengths
# all indels are counted stratified by length : 5-30, 31-100, 100-350, 350-1000, 1000-10000, >10000

# A histogram of alignment lengths is also produced. To do so, all alignment lengths are recorded
# then the histogram is computed.

import argparse
import multiprocessing as mp
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

# Define size bins for indel analysis
SIZE_BINS = [
    (5, 30),
    (31, 100),
    (101, 350),
    (351, 1000),
    (1001, 10000),
    (10001, float("inf")),
]
SIZE_BIN_LABELS = ["5-30bp", "31-100bp", "101-350bp", "351bp-1kb", "1kb-10kb", ">10kb"]


def parse_regions_from_fai(fai_path: Path) -> list[tuple[str, int, int]]:
    """Parse chromosomes/contigs from a fasta index file."""
    regions = []
    with open(fai_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom = parts[0]
            length = int(parts[1])
            regions.append((chrom, 0, length))
    return regions


def process_region(bam_path: Path, region: tuple[str, int, int], min_mapq: int) -> dict:
    """Process a single region and count indels and alignment lengths.
    Large regions are divided into 5Mb sub-regions to avoid memory issues."""
    chrom, start, end = region

    # Counters for indels by size bin
    insertions_by_bin = defaultdict(int)
    deletions_by_bin = defaultdict(int)
    alignment_lengths = []

    # Divide region into 2Mb sub-regions
    sub_region_size = 2_000_000
    region_length = end - start

    # Process sub-regions sequentially
    for sub_start in range(start, end, sub_region_size):
        sub_end = min(sub_start + sub_region_size, end)

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for aln in bam.fetch(chrom, sub_start, sub_end):
                # Skip unmapped, secondary, supplementary, and low quality alignments
                if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                    continue
                if aln.mapping_quality < min_mapq:
                    continue

                # Skip alignments that were already counted in previous sub-region
                # (only count alignments that start in this sub-region)
                if aln.reference_start < sub_start:
                    continue

                # Record alignment length
                alignment_lengths.append(aln.query_length)

                # Extract indels directly from cigar tuples
                for cigar_op, length in aln.cigartuples:
                    # cigar_op: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
                    if cigar_op == 1 and length >= 5:  # Insertion
                        for i, (lower, upper) in enumerate(SIZE_BINS):
                            if lower <= length <= upper:
                                insertions_by_bin[i] += 1
                                break
                    elif cigar_op == 2 and length >= 5:  # Deletion
                        for i, (lower, upper) in enumerate(SIZE_BINS):
                            if lower <= length <= upper:
                                deletions_by_bin[i] += 1
                                break

    return {
        "insertions": dict(insertions_by_bin),
        "deletions": dict(deletions_by_bin),
        "alignment_lengths": alignment_lengths,
    }


def mp_process_region(args):
    """Multiprocessing wrapper."""
    return process_region(*args)


def process_sample(
    bam_path: Path, sample_name: str, regions: list, min_mapq: int, threads: int
) -> dict:
    """Process a single sample's BAM file."""
    print(f"\n=== Processing sample: {sample_name} ===")

    # Prepare jobs for multiprocessing
    jobs = [(bam_path, region, min_mapq) for region in regions]

    # Process regions in parallel
    print(f"Processing alignments with {threads} threads...")
    if threads > 1:
        with mp.Pool(threads) as pool:
            results = list(
                tqdm(
                    pool.imap(mp_process_region, jobs),
                    total=len(jobs),
                    desc=sample_name,
                )
            )
    else:
        results = [
            process_region(bam_path, region, min_mapq)
            for region in tqdm(jobs, desc=sample_name)
        ]

    # Aggregate results
    print("Aggregating results...")
    total_insertions = defaultdict(int)
    total_deletions = defaultdict(int)
    all_alignment_lengths = []

    for result in results:
        for bin_idx, count in result["insertions"].items():
            total_insertions[bin_idx] += count
        for bin_idx, count in result["deletions"].items():
            total_deletions[bin_idx] += count
        all_alignment_lengths.extend(result["alignment_lengths"])

    total_aligned_bases = sum(all_alignment_lengths)

    return {
        "sample_name": sample_name,
        "total_insertions": dict(total_insertions),
        "total_deletions": dict(total_deletions),
        "alignment_lengths": all_alignment_lengths,
        "total_aligned_bases": total_aligned_bases,
    }


def run(args):
    """Main analysis function."""
    # check if histogram has extensions vsg, png, pdf
    valid_exts = [".png", ".svg", ".pdf"]
    if not any(str(args.histogram).lower().endswith(ext) for ext in valid_exts):
        raise ValueError(
            f"Histogram output file must have one of the following extensions: {', '.join(valid_exts)}"
        )
    # Parse input BAM files and sample names
    bam_files = [Path(p.strip()) for p in args.alignments]
    sample_names = [s.strip() for s in args.sample_names]

    if len(bam_files) != len(sample_names):
        raise ValueError(
            f"Number of BAM files ({len(bam_files)}) must match number of sample names ({len(sample_names)})"
        )

    # Parse regions from fai file
    print(f"Parsing regions from {args.fai}")
    regions = parse_regions_from_fai(args.fai)
    print(f"Found {len(regions)} chromosomes/contigs")

    # Process each sample
    sample_results = []
    for bam_path, sample_name in zip(bam_files, sample_names):
        result = process_sample(
            bam_path, sample_name, regions, args.min_mapq, args.threads
        )
        sample_results.append(result)

    # Create combined statistics dataframe
    all_stats = []
    for result in sample_results:
        total_aligned_bases = result["total_aligned_bases"]
        for i, label in enumerate(SIZE_BIN_LABELS):
            ins_count = result["total_insertions"].get(i, 0)
            del_count = result["total_deletions"].get(i, 0)
            all_stats.append(
                {
                    "sample": result["sample_name"],
                    "size_bin": label,
                    "insertions": ins_count,
                    "deletions": del_count,
                    "insertions_per_mb": (
                        (ins_count / total_aligned_bases * 1_000_000)
                        if total_aligned_bases > 0
                        else 0
                    ),
                    "deletions_per_mb": (
                        (del_count / total_aligned_bases * 1_000_000)
                        if total_aligned_bases > 0
                        else 0
                    ),
                }
            )

    stats_df = pd.DataFrame(all_stats)
    stats_df.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved statistics to {args.output}")

    # Create histograms
    print("\nCreating histograms...")
    n_samples = len(sample_results)
    fig, axes = plt.subplots(n_samples, 2, figsize=(14, 5 * n_samples))

    # Ensure axes is always 2D
    if n_samples == 1:
        axes = axes.reshape(1, -1)

    # Color palette for samples
    colors = plt.cm.tab10(np.linspace(0, 1, n_samples))

    for idx, result in enumerate(sample_results):
        sample_name = result["sample_name"]
        all_alignment_lengths = result["alignment_lengths"]
        total_insertions = result["total_insertions"]
        total_deletions = result["total_deletions"]
        total_aligned_bases = result["total_aligned_bases"]

        # Plot 1: Alignment length distribution (capped at 95th percentile)
        ax1 = axes[idx, 0]
        if len(all_alignment_lengths) > 0:
            aln_lengths_array = np.array(all_alignment_lengths)
            percentile_95 = np.percentile(aln_lengths_array, 95)
            capped_lengths = aln_lengths_array[aln_lengths_array <= percentile_95]

            ax1.hist(
                capped_lengths, bins=50, color=colors[idx], edgecolor="black", alpha=0.7
            )
            ax1.set_xlabel("Alignment Length (bp)", fontsize=11)
            ax1.set_ylabel("Count", fontsize=11)
            ax1.set_title(
                f"{sample_name}: Alignment Length Distribution\n(capped at 95th percentile)",
                fontsize=12,
                fontweight="bold",
            )
            ax1.grid(True, alpha=0.3)

            # Add statistics text
            median_len = np.median(aln_lengths_array)
            mean_len = np.mean(aln_lengths_array)
            ax1.text(
                0.98,
                0.97,
                f"N = {len(all_alignment_lengths):,}\nMedian = {median_len:.0f} bp\nMean = {mean_len:.0f} bp\n95th %ile = {percentile_95:.0f} bp",
                transform=ax1.transAxes,
                fontsize=9,
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
            )

        # Plot 2: Indel counts by size bin
        ax2 = axes[idx, 1]
        x_pos = np.arange(len(SIZE_BIN_LABELS))
        width = 0.35

        ins_counts = [total_insertions.get(i, 0) for i in range(len(SIZE_BIN_LABELS))]
        del_counts = [total_deletions.get(i, 0) for i in range(len(SIZE_BIN_LABELS))]

        ax2.bar(
            x_pos - width / 2,
            ins_counts,
            width,
            label="Insertions",
            color="goldenrod",
            edgecolor="black",
            alpha=0.8,
        )
        ax2.bar(
            x_pos + width / 2,
            del_counts,
            width,
            label="Deletions",
            color="slategray",
            edgecolor="black",
            alpha=0.8,
        )

        ax2.set_xlabel("SV Size Bin", fontsize=11)
        ax2.set_ylabel("Count (log scale)", fontsize=11)
        ax2.set_title(
            f"{sample_name}: Indel Counts by Size Bin", fontsize=12, fontweight="bold"
        )
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(SIZE_BIN_LABELS, rotation=45, ha="right")
        ax2.set_yscale("log")
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(args.histogram, dpi=200, bbox_inches="tight")
    print(f"Saved histogram to {args.histogram}")

    # Print summary for all samples
    print("\n=== Summary ===")
    for result in sample_results:
        sample_name = result["sample_name"]
        total_aligned_bases = result["total_aligned_bases"]
        ins_counts = [
            result["total_insertions"].get(i, 0) for i in range(len(SIZE_BIN_LABELS))
        ]
        del_counts = [
            result["total_deletions"].get(i, 0) for i in range(len(SIZE_BIN_LABELS))
        ]

        print(f"\n--- {sample_name} ---")
        print(f"Total alignments analyzed: {len(result['alignment_lengths']):,}")
        print(f"Total aligned bases: {total_aligned_bases:,}")
        print(f"Total insertions: {sum(ins_counts):,}")
        print(f"Total deletions: {sum(del_counts):,}")
        print(
            f"Insertions per Mb: {sum(ins_counts) / total_aligned_bases * 1_000_000:.2f}"
            if total_aligned_bases > 0
            else "N/A"
        )
        print(
            f"Deletions per Mb: {sum(del_counts) / total_aligned_bases * 1_000_000:.2f}"
            if total_aligned_bases > 0
            else "N/A"
        )

    print("\n\nDetailed statistics saved to:", args.output)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Extract alignment statistics from one or more BAM files."
    )
    parser.add_argument(
        "-a",
        "--alignments",
        type=str,
        required=True,
        nargs="+",
        help="list of paths to BAM files (e.g., 'sample1.bam sample2.bam').",
    )
    parser.add_argument(
        "-s",
        "--sample-names",
        type=str,
        required=True,
        nargs="+",
        help="list of sample names corresponding to BAM files (e.g., 'Sample1 Sample2').",
    )
    parser.add_argument(
        "--fai",
        type=Path,
        required=True,
        help="Path to the reference fasta index (.fai) file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output stats file (TSV).",
    )
    parser.add_argument(
        "--histogram",
        type=Path,
        required=True,
        help="Path to the output histogram file (PNG, SVG, PDF).",
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=0,
        help="Minimum mapping quality to consider an alignment.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use per sample.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
