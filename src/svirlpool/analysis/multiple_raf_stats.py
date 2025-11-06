# %%
import argparse
import multiprocessing as mp
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from intervaltree import IntervalTree
from tqdm import tqdm

from ..scripts.util import yield_from_raf


def load_mononucleotide_stretches(bed_file: Path) -> dict[str, IntervalTree]:
    """Load mononucleotide stretches from BED file into intervaltree structure.

    Args:
        bed_file: Path to BED file with mononucleotide stretches

    Returns:
        Dictionary mapping chromosome to IntervalTree of mononucleotide stretches
    """
    intervaltrees = {}

    with open(bed_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])

            if chrom not in intervaltrees:
                intervaltrees[chrom] = IntervalTree()

            # Store the size of the mononucleotide stretch as data
            intervaltrees[chrom][start:end] = end - start

    return intervaltrees


def sort_bed_file(input_bed: Path, output_bed: Path) -> None:
    """Sort a BED file using Unix sort.

    Args:
        input_bed: Path to input BED file
        output_bed: Path to output sorted BED file
    """
    cmd = ["sort", "-k1,1", "-k2,2n", str(input_bed)]
    with open(output_bed, "w") as f:
        subprocess.run(cmd, stdout=f, check=True)


def process_single_raf_file(args: tuple[Path, str]) -> pd.DataFrame:
    """Process a single RAF file and extract SV statistics.

    Args:
        args: tuple of (raf_path, sample_name)

    Returns:
        DataFrame with columns: sample_name, chr, ref_start, ref_end, sv_type, sv_size
    """
    raf_path, sample_name = args

    results = []

    for raf in tqdm(yield_from_raf(raf_path), desc=f"Processing {sample_name}"):
        for sv in raf.SV_signals:
            # Only process insertions (sv_type=0) and deletions (sv_type=1)
            if sv.sv_type >= 2:
                continue

            # Determine SV type
            sv_type = "ins" if sv.sv_type == 0 else "del"

            results.append({
                "sample_name": sample_name,
                "chr": raf.reference_name,
                "ref_start": sv.ref_start,
                "ref_end": sv.ref_end,
                "sv_type": sv_type,
                "sv_size": sv.size,
            })

    return pd.DataFrame(results)


def process_chromosome_bedtools(args: tuple[str, pd.DataFrame, Path, Path]) -> set[int]:
    """Process bedtools intersect for a single chromosome.

    Args:
        args: tuple of (chromosome, sv_df_chr, mononucleotide_bed, tmpdir_path)

    Returns:
        Set of SV indices that overlap with mononucleotide repeats
    """
    chrom, sv_df_chr, mononucleotide_bed, tmpdir_path = args

    # Write SV intervals for this chromosome
    sv_bed_unsorted = tmpdir_path / f"sv_intervals.{chrom}.unsorted.bed"
    sv_bed = tmpdir_path / f"sv_intervals.{chrom}.bed"

    with open(sv_bed_unsorted, "w") as f:
        for idx, row in sv_df_chr.iterrows():
            f.write(f"{row['chr']}\t{row['ref_start']}\t{row['ref_end']}\t{idx}\n")

    # Sort SV BED file
    sort_bed_file(sv_bed_unsorted, sv_bed)

    # Filter mononucleotide BED for this chromosome and sort
    mono_bed_chr = tmpdir_path / f"mononucleotide.{chrom}.bed"
    with open(mononucleotide_bed, "r") as f_in, open(mono_bed_chr, "w") as f_out:
        for line in f_in:
            if line.startswith("#"):
                continue
            if line.startswith(chrom + "\t"):
                f_out.write(line)

    # Sort the filtered mononucleotide file
    mono_bed_chr_sorted = tmpdir_path / f"mononucleotide.{chrom}.sorted.bed"
    sort_bed_file(mono_bed_chr, mono_bed_chr_sorted)

    # Run bedtools intersect
    intersect_output = tmpdir_path / f"intersect_result.{chrom}.bed"
    cmd = [
        "bedtools",
        "intersect",
        "-f",
        "0.5",
        "-a",
        str(sv_bed),
        "-b",
        str(mono_bed_chr_sorted),
        "-wa",
        "-wb",
    ]

    with open(intersect_output, "w") as f:
        subprocess.run(cmd, stdout=f, check=True)

    # Parse intersect results for this chromosome
    # With -f 0.5, bedtools already filters for sufficient overlap
    # We just need to collect the indices that have any overlap
    overlapping_indices = set()

    with open(intersect_output, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue

            # Parse SV entry (column 3 is the index)
            sv_idx_str = parts[3]
            if sv_idx_str != ".":
                sv_idx = int(sv_idx_str)
                overlapping_indices.add(sv_idx)

    return overlapping_indices


def annotate_with_bedtools_intersect(
    sv_df: pd.DataFrame, mononucleotide_bed: Path, threads: int = 4
) -> pd.DataFrame:
    """Annotate SV dataframe with mononucleotide overlap using bedtools intersect.

    Parallelizes by chromosome for faster processing.

    Args:
        sv_df: DataFrame with SV information (chr, ref_start, ref_end, sv_size)
        mononucleotide_bed: Path to BED file with mononucleotide stretches
        threads: Number of threads for parallel processing

    Returns:
        DataFrame with additional column 'in_mononucleotide_repeat'
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        # Group SVs by chromosome
        print("Grouping SVs by chromosome...")
        chromosomes = sv_df["chr"].unique()
        print(
            f"Processing {len(chromosomes)} chromosomes in parallel with {threads} threads..."
        )

        # Prepare arguments for parallel processing
        process_args = []
        for chrom in chromosomes:
            sv_df_chr = sv_df.query("chr == @chrom")
            process_args.append((chrom, sv_df_chr, mononucleotide_bed, tmpdir_path))

        # Process chromosomes in parallel
        if threads > 1:
            with mp.Pool(processes=threads) as pool:
                results = list(
                    tqdm(
                        pool.imap(process_chromosome_bedtools, process_args),
                        total=len(process_args),
                        desc="Processing chromosomes with bedtools",
                    )
                )
        else:
            results = [
                process_chromosome_bedtools(arg)
                for arg in tqdm(process_args, desc="Processing chromosomes")
            ]

        # Combine results from all chromosomes
        print("Combining results from all chromosomes...")
        all_overlapping_indices = set()
        for overlapping_indices in results:
            all_overlapping_indices.update(overlapping_indices)

        # Add annotation column
        sv_df["in_mononucleotide_repeat"] = False
        sv_df.loc[list(all_overlapping_indices), "in_mononucleotide_repeat"] = True

        print(
            f"Found {len(all_overlapping_indices)} SVs overlapping with mononucleotide repeats"
        )

    return sv_df


def process_multiple_rafs(
    raf_paths: list[Path],
    sample_names: list[str],
    mononucleotide_bed: Path,
    output_tsv: Path | None = None,
    threads: int = 4,
) -> pd.DataFrame:
    """Process multiple RAF files in parallel and combine results.

    Args:
        raf_paths: list of paths to RAF files
        sample_names: list of sample names corresponding to RAF files
        mononucleotide_bed: Path to BED file with mononucleotide stretches
        output_tsv: Path to output TSV file
        threads: Number of threads for multiprocessing

    Returns:
        Combined DataFrame with all SV statistics
    """
    if len(raf_paths) != len(sample_names):
        raise ValueError("Number of RAF paths must match number of sample names")

    print(f"Preparing to process {len(raf_paths)} RAF files...")
    # Prepare arguments for multiprocessing
    process_args = [
        (raf_path, sample_name)
        for raf_path, sample_name in zip(raf_paths, sample_names)
    ]

    # Process files in parallel
    print(f"Processing {len(raf_paths)} RAF files with {threads} threads...")
    if threads > 1:
        with mp.Pool(processes=threads) as pool:
            dfs = list(
                tqdm(
                    pool.imap(process_single_raf_file, process_args),
                    total=len(process_args),
                    desc="Processing RAF files",
                )
            )
    else:
        dfs = [
            process_single_raf_file(arg)
            for arg in tqdm(process_args, desc="Processing RAF files")
        ]

    # Combine all DataFrames
    print("Combining results...")
    combined_df = pd.concat(dfs, ignore_index=True)

    # Annotate with mononucleotide overlaps using bedtools (parallelized by chromosome)
    print("Annotating with mononucleotide repeat overlaps...")
    combined_df = annotate_with_bedtools_intersect(
        combined_df, mononucleotide_bed, threads=threads
    )

    # Save to TSV
    if output_tsv is not None:
        print(f"Saving results to {output_tsv}...")
        combined_df.to_csv(output_tsv, sep="	", index=False)

    print(
        f"Done! Processed {len(combined_df)} SV signals from {len(raf_paths)} samples"
    )
    return combined_df


def run(args, **kwargs):
    """Run the multiple RAF stats analysis."""
    process_multiple_rafs(
        raf_paths=args.rafs,
        sample_names=args.samples,
        mononucleotide_bed=args.mononucleotide_bed,
        output_tsv=args.output,
        threads=args.threads,
    )


def get_parser():
    """Create argument parser for multiple RAF statistics."""
    parser = argparse.ArgumentParser(
        description="Extract SV statistics from multiple RAF files and check overlap with mononucleotide repeats."
    )
    parser.add_argument(
        "-r",
        "--rafs",
        type=Path,
        nargs="+",
        required=True,
        help="List of paths to RAF files",
    )
    parser.add_argument(
        "-s",
        "--samples",
        type=str,
        nargs="+",
        required=True,
        help="List of sample names corresponding to RAF files",
    )
    parser.add_argument(
        "-m",
        "--mononucleotide-bed",
        type=Path,
        required=True,
        help="Path to BED file with mononucleotide stretches",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output path for TSV file with SV statistics",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of threads for multiprocessing (default: 4)",
    )
    return parser


def main():
    """Main entry point."""
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
