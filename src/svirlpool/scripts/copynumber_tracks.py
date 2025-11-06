"""
Copy Number Track Generation Module

This module generates copy number tracks from coverage data using HMM/Markov models
and stores them in both BED format (for random access via tabix) and SQLite database
(for integration with svirltile).
"""

# %%
import argparse
import gzip
import logging
import pickle
import sqlite3
import subprocess
import time
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pysam
from intervaltree import IntervalTree
from scipy import stats
from tqdm import tqdm

from . import covtree, rafs_to_coverage

log = logging.getLogger(__name__)

# %%


def create_copynumber_database(output_db_path: Path) -> None:
    """
    Create a database to store copy number interval trees.

    Args:
        output_db_path: Path to the output database file
    """
    with sqlite3.connect(str(output_db_path)) as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS copynumber_tracks (
                chromosome VARCHAR(50) PRIMARY KEY,
                intervaltree BLOB NOT NULL
            )
        """
        )
        conn.commit()
    log.debug(f"Created copy number database: {output_db_path}")


def write_copynumber_trees_to_db(
    cn_trees: dict[str, IntervalTree], db_path: Path
) -> None:
    """
    Write copy number interval trees to database.

    Args:
        cn_trees: Dictionary mapping chromosome names to IntervalTrees
        db_path: Path to the database file
    """
    # Create database if it doesn't exist
    create_copynumber_database(db_path)

    with sqlite3.connect(str(db_path)) as conn:
        for chromosome, tree in cn_trees.items():
            # Serialize the IntervalTree using pickle
            serialized_tree = pickle.dumps(tree)

            # Insert or replace the tree in the database
            conn.execute(
                """
                INSERT OR REPLACE INTO copynumber_tracks
                (chromosome, intervaltree)
                VALUES (?, ?)
            """,
                (chromosome, serialized_tree),
            )

        # Create index for faster lookups
        conn.execute(
            """
            CREATE INDEX IF NOT EXISTS idx_copynumber_chromosome
            ON copynumber_tracks (chromosome)
        """
        )
        conn.commit()

    log.info(f"Wrote {len(cn_trees)} chromosome copy number tracks to {db_path}")


def load_copynumber_trees_from_db(
    db_path: Path, chromosomes: list[str] | None = None
) -> dict[str, IntervalTree]:
    """
    Load copy number interval trees from database.

    Args:
        db_path: Path to the database file
        chromosomes: Optional list of chromosomes to load (loads all if None)

    Returns:
        Dictionary mapping chromosome names to IntervalTrees
    """
    cn_trees = {}

    with sqlite3.connect(f"file:{str(db_path)}?mode=ro", uri=True) as conn:
        cursor = conn.cursor()

        if chromosomes is None:
            # Load all chromosomes
            cursor.execute("SELECT chromosome, intervaltree FROM copynumber_tracks")
        else:
            # Load specific chromosomes
            placeholders = ",".join("?" * len(chromosomes))
            cursor.execute(
                f"SELECT chromosome, intervaltree FROM copynumber_tracks WHERE chromosome IN ({placeholders})",
                chromosomes,
            )

        for chromosome, serialized_tree in cursor.fetchall():
            cn_trees[chromosome] = pickle.loads(serialized_tree)

        cursor.close()

    log.info(f"Loaded {len(cn_trees)} chromosome copy number tracks from {db_path}")
    return cn_trees


# ============================================================================
# Copy Number Computation Functions
# ============================================================================


def _parse_bed_file(bed_path: Path | str) -> dict[str, list[tuple[int, int]]]:
    """Parse BED file and return dictionary of chromosome -> list of (start, end) tuples."""
    regions = {}
    with open(bed_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chr_name = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if chr_name not in regions:
                regions[chr_name] = []
            regions[chr_name].append((start, end))
    return regions


def _to_numpy_per_base_coverage(
    data: np.ndarray, chr_filterlist: list[str] = None
) -> dict[str, np.ndarray]:
    """Converts coverage data to per-base numpy arrays for each chromosome."""
    if chr_filterlist is None:
        chr_filterlist = []
    chrs = set(data[:, 0]) - set(chr_filterlist)
    cov_arrays = {}
    for chr in tqdm(chrs, desc="Converting to per-base coverage arrays"):
        chr = str(chr)
        data_chr_intervals = data[data[:, 0] == chr][:, 1:3].astype(int)
        _data_chr_readnames = data[data[:, 0] == chr][:, 3].astype(
            str
        )  # FIXME: unused?
        if len(data_chr_intervals) == 0:
            continue
        max_pos = np.max(data_chr_intervals)
        coverage_array = np.zeros(max_pos, dtype=np.int32)
        for start, end in tqdm(
            data_chr_intervals, desc=f"Processing {chr}", leave=False
        ):
            coverage_array[start:end] += 1
        cov_arrays[chr] = coverage_array
    return cov_arrays


def _to_numpy_per_base_coverage_with_regions(
    data: np.ndarray,
    regions: dict[str, list[tuple[int, int]]],
    chr_filterlist: list[str] = None,
) -> dict[str, np.ndarray]:
    """Converts coverage data to per-base numpy arrays for specified regions only."""
    if chr_filterlist is None:
        chr_filterlist = []
    chrs = set(data[:, 0]) - set(chr_filterlist)
    cov_arrays = {}

    for chr in tqdm(chrs, desc="Converting to per-base coverage arrays (regions)"):
        chr = str(chr)
        if chr not in regions:
            continue

        data_chr_intervals = data[data[:, 0] == chr][:, 1:3].astype(int)
        if len(data_chr_intervals) == 0:
            continue

        # Create coverage array for each region
        region_coverages = []
        for region_start, region_end in regions[chr]:
            region_length = region_end - region_start
            region_coverage = np.zeros(region_length, dtype=np.int32)

            # Find overlapping intervals
            for start, end in data_chr_intervals:
                # Calculate overlap with region
                overlap_start = max(start, region_start)
                overlap_end = min(end, region_end)

                if overlap_start < overlap_end:
                    # Convert to region-relative coordinates
                    rel_start = overlap_start - region_start
                    rel_end = overlap_end - region_start
                    region_coverage[rel_start:rel_end] += 1

            region_coverages.append(region_coverage)

        # Concatenate all regions for this chromosome
        if region_coverages:
            cov_arrays[chr] = np.concatenate(region_coverages)

    return cov_arrays


def _bins_from_per_base_coverage(
    cov_arrays: dict[str, np.ndarray], bin_size: int
) -> list[int]:
    """Bins per-base coverage arrays and returns a list of binned coverages."""
    all_bin_coverages = []
    for _chr, coverage_array in cov_arrays.items():
        max_pos = len(coverage_array)
        num_bins = (max_pos + bin_size - 1) // bin_size  # Ceiling division
        for i in range(num_bins):
            bin_start = i * bin_size
            bin_end = min((i + 1) * bin_size, max_pos)
            bin_coverage = coverage_array[bin_start:bin_end]
            if len(bin_coverage) > 0:
                max_bin_cov = np.max(bin_coverage)
                all_bin_coverages.append(max_bin_cov)
    return all_bin_coverages


def _compute_emission_probability(
    observed_cov: int, expected_cov: float, dispersion: float = 0.1
) -> float:
    """
    Computes emission probability using Negative Binomial distribution.

    The Negative Binomial better models sequencing coverage which shows overdispersion
    (variance > mean) due to PCR biases, GC content, mappability, etc.

    Args:
        observed_cov: Observed coverage at a bin
        expected_cov: Expected coverage for this CN state
        dispersion: Dispersion parameter (higher = more variance).
                   Typical values: 0.05-0.2 for well-behaved samples

    Returns:
        Emission probability P(observed | expected)
    """
    if expected_cov <= 0.1:
        # For CN=0, we expect very low coverage
        # Use Poisson with small lambda
        return stats.poisson.pmf(observed_cov, mu=0.1)

    # Negative Binomial parameterization:
    # mean = expected_cov
    # variance = expected_cov + dispersion * expected_cov^2

    # Convert to (n, p) parameterization used by scipy
    # n = mean^2 / (variance - mean)
    # p = mean / variance

    var = expected_cov + dispersion * (expected_cov**2)

    if var <= expected_cov:
        # Fallback to Poisson if variance is too small
        return stats.poisson.pmf(observed_cov, mu=expected_cov)

    n = (expected_cov**2) / (var - expected_cov)
    p = expected_cov / var

    # Ensure valid parameters
    n = max(n, 1e-6)
    p = min(max(p, 1e-6), 1.0 - 1e-6)

    return stats.nbinom.pmf(observed_cov, n=n, p=p)


def _viterbi_copy_number(
    bin_coverages: list[int],
    median_cov: float,
    transition_probs: np.ndarray,
    dispersion: float = 0.1,
) -> list[int]:
    """
    Viterbi algorithm to find most likely copy number sequence.

    Uses Negative Binomial emission probabilities with strong bias toward CN=2.

    Args:
        bin_coverages: List of coverage values for each bin
        median_cov: Median coverage of the genome (corresponds to CN=2)
        transition_probs: Transition probability matrix (CN states 0-4)
        dispersion: Overdispersion parameter for Negative Binomial (default: 0.1)

    Returns:
        List of copy number assignments for each bin
    """
    # Define copy number states (0, 1, 2, 3, 4)
    states = [0, 1, 2, 3, 4]
    n_states = len(states)
    n_obs = len(bin_coverages)

    # Expected coverage for each CN state
    expected_cov = {
        0: 0.01,  # Near-zero coverage
        1: median_cov * 0.6,
        2: median_cov * 1.0,
        3: median_cov * 1.6,
        4: median_cov * 2.2,
    }

    # Initialize Viterbi tables
    viterbi = np.zeros((n_obs, n_states))
    backpointer = np.zeros((n_obs, n_states), dtype=int)

    initial_probs = np.array([0.2, 0.2, 0.2, 0.2, 0.2])

    # First observation
    for s in range(n_states):
        cn = states[s]
        emission_prob = _compute_emission_probability(
            bin_coverages[0], expected_cov[cn], dispersion
        )

        # Apply prior penalty in log space
        viterbi[0, s] = np.log(initial_probs[s] + 1e-10) + np.log(emission_prob + 1e-10)

    # Forward pass
    for t in range(1, n_obs):
        for s in range(n_states):
            cn = states[s]

            # Compute transition probabilities from all previous states
            trans_probs = viterbi[t - 1, :] + np.log(transition_probs[:, s] + 1e-10)
            backpointer[t, s] = np.argmax(trans_probs)
            max_trans_prob = np.max(trans_probs)

            emission_prob = _compute_emission_probability(
                bin_coverages[t], expected_cov[cn], dispersion
            )

            # Apply prior penalty in log space
            viterbi[t, s] = max_trans_prob + np.log(emission_prob + 1e-10)

    # Backtrack to find best path
    best_path = np.zeros(n_obs, dtype=int)
    best_path[-1] = np.argmax(viterbi[-1, :])

    for t in range(n_obs - 2, -1, -1):
        best_path[t] = backpointer[t + 1, best_path[t + 1]]

    # Convert state indices to copy numbers
    return [states[s] for s in best_path]


def _get_default_transition_matrix(stay_prob: float) -> np.ndarray:
    """
    Returns default transition probability matrix for copy number states.
    States: CN=0, CN=1, CN=2, CN=3, CN=4

    The matrix is computed from two probabilities:
    - stay_prob: Probability of staying in the current state
    - change_prob: Probability of changing to any other state (distributed equally)

    Args:
        stay_prob: Probability of remaining in current state (default: 0.96)
                   Higher values make the HMM more conservative about CN changes.

    Returns:
        5x5 transition probability matrix
    """
    n_states = 5
    change_prob = 1.0 - stay_prob
    # Distribute change probability equally among all other states
    off_diagonal_prob = change_prob / (n_states - 1)

    # Create matrix with off-diagonal probabilities
    trans_matrix = np.full((n_states, n_states), off_diagonal_prob)
    # Set diagonal to stay probability
    np.fill_diagonal(trans_matrix, stay_prob)

    return trans_matrix


def median_binned_total_coverage(
    path_db: Path | str,
    bin_size: int,
    chr_filterlist: list[str] = None,
    regions: Path | str | None = None,
) -> int:
    """Computes median coverage of all bins across all chromosomes, except those in chr_filterlist."""
    if chr_filterlist is None:
        chr_filterlist = []
    data = rafs_to_coverage.load_cov_from_db(path_db)
    data = np.array(data)

    if regions is not None:
        bed_regions = _parse_bed_file(regions)
        np_data = _to_numpy_per_base_coverage_with_regions(
            data, bed_regions, chr_filterlist
        )
    else:
        np_data = _to_numpy_per_base_coverage(data, chr_filterlist)

    bins = _bins_from_per_base_coverage(np_data, bin_size)
    return int(np.median(bins))


def load_copynumber_track_as_intervaltree(
    bgzip_bed: Path | str, chromosomes: list[str] | None = None
) -> dict[str, IntervalTree]:
    """
    Loads a bgzipped copy number BED file into interval trees for random access.

    Args:
        bgzip_bed: Path to bgzipped and tabix-indexed BED file
        chromosomes: Optional list of chromosomes to load (loads all if None)

    Returns:
        Dictionary mapping chromosome names to IntervalTrees with CN data
    """
    cn_trees = {}

    with gzip.open(bgzip_bed, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue

            chr_name = parts[0]
            if chromosomes is not None and chr_name not in chromosomes:
                continue

            start = int(parts[1])
            end = int(parts[2])
            cn = int(parts[3])

            if chr_name not in cn_trees:
                cn_trees[chr_name] = IntervalTree()

            cn_trees[chr_name].addi(start, end, cn)

    return cn_trees


def query_copynumber_from_regions(
    bgzip_bed: Path | str, regions: list[tuple[str, int, int]]
) -> int:
    """
    Query copy numbers from a bgzipped BED file for given regions and return the maximum CN.

    Uses tabix for efficient random access to query specific regions.

    Args:
        bgzip_bed: Path to bgzipped and tabix-indexed copy number BED file
        regions: List of tuples (chromosome, start, end) defining regions to query

    Returns:
        Maximum copy number found across all queried regions
    """
    if len(regions) == 0:
        log.warning("No regions provided for copy number query, returning default CN=2")
        return 2

    max_cn = 0

    # Use context manager for automatic file closure
    with pysam.TabixFile(str(bgzip_bed)) as tabix_file:
        for chr_name, start, end in regions:
            try:
                # Query the region using tabix
                # fetch() returns an iterator over tab-delimited rows
                for row in tabix_file.fetch(chr_name, start, end):
                    parts = row.split("\t")
                    if len(parts) >= 4:
                        cn = int(parts[3])
                        max_cn = max(max_cn, cn)
            except (ValueError, KeyError) as e:
                # ValueError: chromosome not in index
                # KeyError: chromosome not found
                log.warning(f"Could not query region {chr_name}:{start}-{end}: {e}")
                continue

    # Ensure at least CN=2 (diploid baseline)
    if max_cn == 0:
        log.warning("No copy number data found for regions, returning default CN=2")
        return 2

    log.info(f"Maximum copy number across {len(regions)} region(s): {max_cn}")
    return max_cn


# ============================================================================
# Main Copy Number Track Generation
# ============================================================================


def generate_copynumber_tracks(
    fasta_index: Path,
    coverage_db: Path,
    output_bed: Path,
    output_db: Path,
    threads: int,
    stay_prob: float,
    dispersion: float,
    chr_filterlist: list[str] = None,
    regions_path: Path | None = None,
    transition_matrix=None,
) -> tuple[Path, Path]:
    """
    Generate copy number tracks from coverage data.

    Creates two outputs:
    1. A bgzipped and tabix-indexed BED file for random access
    2. A SQLite database with IntervalTree representation for fast querying

    Args:
        coverage_db: Path to coverage database (from rafs_to_coverage)
        output_bed: Path to output BED file (will be bgzipped)
        output_db: Path to output database file
        chr_filterlist: Chromosomes to exclude from analysis
        regions: Optional BED file with regions to restrict analysis
        transition_matrix: Optional custom transition probability matrix
        dispersion: Overdispersion parameter for Negative Binomial emission probabilities
                   (default: 0.1). Higher = more tolerance for variance. Range: 0.05-0.2

    Returns:
        Tuple of (bgzipped_bed_path, database_path)
    """
    if chr_filterlist is None:
        chr_filterlist = []
    log.info("Generating copy number tracks...")

    # step 1: use covtree to compute non-overlapping intervals of coverage
    log.info(f"loading data from {coverage_db}")
    data = covtree.load_reference_data(path_db=coverage_db)
    log.info("constructing interval trees")
    intervall_trees_reads, all_positions = covtree.construct_interval_trees(data=data)
    log.info("computing coverages")
    start_time = time.time()
    cov_trees = covtree.parallel_coverage_computation(
        all_positions=all_positions,
        intervall_trees=intervall_trees_reads,
        num_workers=threads,
    )
    log.info(f"Execution time: {time.time() - start_time:.2f} seconds")

    # step 2a: if no regions file is provided, generate regions from the fasta index
    if regions_path is None:
        log.info(
            f"No regions file provided, generating regions from fasta index: {fasta_index}"
        )
        chromosomes = _parse_fai_file(fasta_index)
        regions = {}
        for chr_name, chr_length in chromosomes:
            if chr_name in chr_filterlist:
                continue
            regions[chr_name] = [(0, chr_length)]
        log.info(f"Generated regions for {len(regions)} chromosomes from fasta index")
    else:  # parse bed file to regions
        log.info(f"Parsing regions from BED file: {regions_path}")
        regions = _parse_bed_file(regions_path)
        log.info(f"Parsed regions for {len(regions)} chromosomes from BED file")

    # step 3: compute the genome-wide median coverage from all bins (intervals) across all chromosomes. An interval is of the form (start, end, coverage)
    all_coverages = []
    for _chr, tree in cov_trees.items():
        all_coverages.extend([interval.data for interval in tree])
    median_cov = int(np.median(all_coverages))

    # step 4: compute the median read length from intervall_trees_reads
    all_read_lengths = []
    for _chr, tree in intervall_trees_reads.items():
        all_read_lengths.extend([interval.end - interval.begin for interval in tree])

    # step 5: new bin size is 10 * median read length
    median_read_length = int(np.median(all_read_lengths))
    new_bin_size = 10 * median_read_length

    # regions might be separated by small gaps (smaller than median_read_length)
    # if that is the case, regions should be merged.
    _merge_close_regions_inplace(regions=regions, median_read_length=median_read_length)

    # step 6: re-bin the coverage intervals to the new bin size, taking the median coverage in each bin
    # 1) for each region, get the coverage intervals from cov_trees
    # 2) generate a new chr:intervalltree datas tructure that holds all new big bins with median coverage
    binned_cov_trees: dict[str, IntervalTree] = {}  # nested interval trees of sub trees
    for chr in tqdm(sorted(cov_trees.keys()), desc="Re-binning coverage intervals"):
        if chr not in regions:
            continue
        binned_cov_trees[chr] = IntervalTree()
        chr_tree = cov_trees[chr]
        for region_start, region_end in regions[chr]:
            # Generate bins within this region
            sub_it: list[tuple[int, int, int]] = []
            for bin_start in range(region_start, region_end, new_bin_size):
                bin_end = min(bin_start + new_bin_size, region_end)
                # Get all coverage intervals overlapping this bin
                overlapping_intervals = chr_tree[bin_start:bin_end]
                if len(overlapping_intervals) == 0:
                    median_bin_cov = 0
                else:
                    median_bin_cov = int(
                        np.median([interval.data for interval in overlapping_intervals])
                    )
                sub_it.append((bin_start, bin_end, median_bin_cov))
            sub_tree = IntervalTree.from_tuples(sub_it)
            binned_cov_trees[chr].addi(region_start, region_end, sub_tree)

    # now we have a data structure of binned_cov_trees[chr] = IntervalTree of (start, end, IntervalTree of (bin_start, bin_end, median_cov))
    # this allows viterbi to be executed on each region separately.
    # The assumption is that each region consists of consecutive intervals,
    # while the regions themselves may be separated by gaps. close-by regions have been merged already.

    # step 7: Get transition matrix
    if transition_matrix is None:
        transition_matrix = _get_default_transition_matrix(stay_prob=stay_prob)

    # step 8: Apply Viterbi algorithm per region on each chromosome
    log.info(
        f"Applying Viterbi algorithm to determine copy numbers (dispersion={dispersion})"
    )
    cn_tracks = {}

    for chr in tqdm(sorted(binned_cov_trees.keys()), desc="Processing chromosomes"):
        cn_tracks[chr] = []

        # Iterate over each region in this chromosome
        for region_interval in sorted(binned_cov_trees[chr], key=lambda x: x.begin):
            # Extract the sub-tree for this region
            sub_tree = region_interval.data

            # Convert sub-tree to sorted list of (start, end, coverage) tuples
            bins = sorted(
                [(iv.begin, iv.end, iv.data) for iv in sub_tree], key=lambda x: x[0]
            )

            if len(bins) == 0:
                continue

            # Extract coverage values for Viterbi
            bin_coverages = [cov for _, _, cov in bins]

            # Apply Viterbi to this region
            copy_numbers = _viterbi_copy_number(
                bin_coverages=bin_coverages,
                median_cov=median_cov,
                transition_probs=transition_matrix,
                dispersion=dispersion,
            )

            # Combine positions with copy numbers
            for (start, end, _), cn in zip(bins, copy_numbers, strict=True):
                cn_tracks[chr].append((start, end, cn))

    # Write to BED file
    log.info(f"Writing copy number track to {output_bed}")
    output_path = Path(output_bed)
    with open(output_path, "w") as f:
        for chr in sorted(cn_tracks.keys()):
            for start, end, cn in cn_tracks[chr]:
                f.write(f"{chr}\t{start}\t{end}\t{cn}\n")

    # Bgzip compress
    log.info("Compressing with bgzip")
    bgzip_path = Path(str(output_path) + ".gz")
    subprocess.run(["bgzip", "-f", str(output_path)], check=True)

    # Tabix index
    log.info("Indexing with tabix")
    subprocess.run(["tabix", "-f", "-p", "bed", str(bgzip_path)], check=True)

    # Step 2: Load the generated BED file as interval trees
    log.info("Loading copy number tracks into interval trees...")
    cn_trees = load_copynumber_track_as_intervaltree(bgzip_path)

    # Step 3: Write interval trees to database
    log.info("Writing copy number interval trees to database...")
    write_copynumber_trees_to_db(cn_trees, output_db)

    log.info("Copy number tracks generated successfully:")
    log.info(f"  - BED file: {bgzip_path}")
    log.info(f"  - Database: {output_db}")

    return bgzip_path, output_db


def _merge_close_regions_inplace(
    regions: dict[str, list[tuple[int, int]]], median_read_length: int
) -> None:
    """Merge regions that are closer than median_read_length."""
    for chr in regions:
        if chr not in regions:
            continue
        merged_regions = []
        sorted_regions = sorted(regions[chr], key=lambda x: x[0])
        current_start, current_end = sorted_regions[0]
        for start, end in sorted_regions[1:]:
            if start <= current_end + median_read_length:
                # Merge regions
                current_end = max(current_end, end)
            else:
                merged_regions.append((current_start, current_end))
                current_start, current_end = start, end
        merged_regions.append((current_start, current_end))
        regions[chr] = merged_regions


def _parse_fai_file(fai_path: Path) -> list[tuple[str, int]]:
    chromosomes = []
    with open(fai_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chr_name = parts[0]
                chr_length = int(parts[1])
                chromosomes.append((chr_name, chr_length))
    return chromosomes


def plot_copynumber_tracks(
    cn_trees: dict[str, IntervalTree],
    output_figure: Path,
    reference_fai: Path,
    title: str = "Copy Number Tracks",
) -> None:
    """
    Plot copy number tracks as horizontal colored lines per chromosome.

    Args:
        cn_trees: Dictionary mapping chromosome names to IntervalTrees with CN data
        output_figure: Path to save the output figure
        reference_fai: Path to reference FASTA index (.fai) for chromosome order/sizes
        title: Title for the plot
    """
    # Color scheme for copy numbers
    CN_COLORS = {
        0: "#DC143C",  # Crimson red (deletion)
        1: "#FF8C00",  # Dark orange (hemizygous)
        2: "#000000",  # Black (diploid normal)
        3: "#1E90FF",  # Dodger blue (triplication)
        4: "#006400",  # Dark green (quadruplication)
        5: "#004d00",  # Even darker green (5+)
    }
    NO_DATA_COLOR = "#D3D3D3"  # Light gray

    # Parse chromosome order and sizes from .fai
    chromosomes = _parse_fai_file(reference_fai)

    if len(chromosomes) == 0:
        log.error(f"No chromosomes found in {reference_fai}")
        return

    # Set up the figure
    n_chrom = len(chromosomes)
    fig_height = max(
        0.3 * n_chrom, 6
    )  # At least 6 inches, scale with number of chromosomes
    fig, ax = plt.subplots(figsize=(16, fig_height))

    # Track maximum chromosome length for x-axis
    max_length = max(length for _, length in chromosomes)

    # Plot each chromosome
    for idx, (chr_name, chr_length) in enumerate(chromosomes):
        y_pos = n_chrom - idx - 1  # Draw from top to bottom

        if chr_name not in cn_trees or len(cn_trees[chr_name]) == 0:
            # No data - draw gray line for entire chromosome
            ax.plot(
                [0, chr_length],
                [y_pos, y_pos],
                color=NO_DATA_COLOR,
                linewidth=4,
                solid_capstyle="butt",
            )
        else:
            # Get all intervals for this chromosome and sort by start position
            intervals = sorted(cn_trees[chr_name], key=lambda x: x.begin)

            if len(intervals) == 0:
                # Empty tree - draw gray
                ax.plot(
                    [0, chr_length],
                    [y_pos, y_pos],
                    color=NO_DATA_COLOR,
                    linewidth=4,
                    solid_capstyle="butt",
                )
            else:
                # Group consecutive intervals with same CN
                segments = []
                current_cn = intervals[0].data
                current_start = intervals[0].begin
                current_end = intervals[0].end

                for interval in intervals[1:]:
                    if interval.data == current_cn and interval.begin <= current_end:
                        # Same CN and consecutive/overlapping - extend segment
                        current_end = max(current_end, interval.end)
                    else:
                        # Different CN or gap - save current segment and start new one
                        segments.append((current_start, current_end, current_cn))
                        current_cn = interval.data
                        current_start = interval.begin
                        current_end = interval.end

                # Add the last segment
                segments.append((current_start, current_end, current_cn))

                # Draw segments
                for start, end, cn in segments:
                    # Map CN to color (CN > 4 gets same color as CN=4+)
                    color = CN_COLORS.get(min(cn, 5), CN_COLORS[2])
                    ax.plot(
                        [start, end],
                        [y_pos, y_pos],
                        color=color,
                        linewidth=4,
                        solid_capstyle="butt",
                    )

                # Fill gaps with gray (regions not covered by CN calls)
                if segments[0][0] > 0:
                    # Gap at the beginning
                    ax.plot(
                        [0, segments[0][0]],
                        [y_pos, y_pos],
                        color=NO_DATA_COLOR,
                        linewidth=4,
                        solid_capstyle="butt",
                    )

                for i in range(len(segments) - 1):
                    gap_start = segments[i][1]
                    gap_end = segments[i + 1][0]
                    if gap_end > gap_start:
                        ax.plot(
                            [gap_start, gap_end],
                            [y_pos, y_pos],
                            color=NO_DATA_COLOR,
                            linewidth=4,
                            solid_capstyle="butt",
                        )

                if segments[-1][1] < chr_length:
                    # Gap at the end
                    ax.plot(
                        [segments[-1][1], chr_length],
                        [y_pos, y_pos],
                        color=NO_DATA_COLOR,
                        linewidth=4,
                        solid_capstyle="butt",
                    )

        # Add chromosome label
        ax.text(
            -max_length * 0.02, y_pos, chr_name, ha="right", va="center", fontsize=9
        )

        # Add tick marks every 25 Mb
        tick_interval = 25_000_000  # 25 Mb
        for tick_pos in range(0, chr_length, tick_interval):
            ax.plot(
                [tick_pos, tick_pos],
                [y_pos - 0.15, y_pos + 0.15],
                color="black",
                linewidth=0.5,
                zorder=10,
            )

    # Formatting
    ax.set_xlim(-max_length * 0.05, max_length * 1.02)
    ax.set_ylim(-0.5, n_chrom - 0.5)
    ax.set_xlabel("Position (bp)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")

    # Format x-axis to show Mb
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{int(x / 1e6)}"))
    ax.set_xlabel("Position (Mb)", fontsize=12)

    # Remove y-axis ticks (chromosome labels are on the left)
    ax.set_yticks([])

    # Add legend
    legend_elements = [
        mpatches.Patch(color=CN_COLORS[0], label="CN=0 (Deletion)"),
        mpatches.Patch(color=CN_COLORS[1], label="CN=1 (Hemizygous)"),
        mpatches.Patch(color=CN_COLORS[2], label="CN=2 (Diploid)"),
        mpatches.Patch(color=CN_COLORS[3], label="CN=3 (Triplication)"),
        mpatches.Patch(color=CN_COLORS[4], label="CN=4+"),
        mpatches.Patch(color=NO_DATA_COLOR, label="No data"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=9)

    # Grid for easier reading
    ax.grid(axis="x", alpha=0.3, linestyle="--", linewidth=0.5)

    # Tight layout
    plt.tight_layout()

    # Save figure
    plt.savefig(output_figure, dpi=150, bbox_inches="tight")
    plt.close()

    log.info(f"Copy number track plot saved to {output_figure}")


def get_parser():
    parser = argparse.ArgumentParser(
        description="Generate copy number tracks from coverage data"
    )
    parser.add_argument(
        "--coverage-db",
        type=Path,
        required=True,
        help="Path to coverage database (from rafs_to_coverage)",
    )
    parser.add_argument(
        "--output-bed",
        type=Path,
        required=True,
        help="Path to output BED file (will be bgzipped and indexed)",
    )
    parser.add_argument(
        "--output-db", type=Path, required=True, help="Path to output database file"
    )
    parser.add_argument(
        "--threads", type=int, default=4, help="Number of threads to use (default: 4)"
    )
    parser.add_argument(
        "--dispersion",
        type=float,
        default=0.1,
        help="Overdispersion parameter for Negative Binomial emission model (default: 0.1). "
        "Higher values = more tolerance for coverage variance. Typical range: 0.05-0.2. "
        "Use ~0.05-0.1 for clean data, ~0.15-0.2 for noisy data.",
    )
    parser.add_argument(
        "--stay-prob",
        type=float,
        default=0.96,
        help="Stay probability for the Hidden Markov Model (default: 0.98). "
        "Higher values = more likely to stay in the same state.",
    )
    parser.add_argument(
        "--chr-filterlist",
        nargs="*",
        default=[],
        help="Chromosomes to exclude from analysis",
    )
    parser.add_argument(
        "--regions",
        type=Path,
        default=None,
        help="Optional BED file with regions to restrict analysis",
    )
    parser.add_argument(
        "--plot-output",
        type=Path,
        default=None,
        help="Optional path to save a copy number track visualization plot",
    )
    parser.add_argument(
        "--reference-fai", type=Path, help="Path to reference FASTA index (.fai)d"
    )
    parser.add_argument(
        "--plot-title",
        type=str,
        default="Copy Number Tracks",
        help="Title for the plot (if --plot-output is specified)",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Validate dispersion parameter
    if args.dispersion <= 0:
        log.error(f"--dispersion must be positive, got {args.dispersion}")
        return
    if args.dispersion > 0.5:
        log.warning(
            f"--dispersion={args.dispersion} is unusually high (typical: 0.05-0.2). "
            "This may lead to overly permissive CN calls."
        )

    # Generate copy number tracks
    generate_copynumber_tracks(
        fasta_index=args.reference_fai,
        coverage_db=args.coverage_db,
        output_bed=args.output_bed,
        output_db=args.output_db,
        chr_filterlist=args.chr_filterlist,
        regions_path=args.regions,
        threads=args.threads,
        dispersion=args.dispersion,
        stay_prob=args.stay_prob,
    )

    # Generate plot if requested
    if args.plot_output is not None:
        log.info("Generating copy number track plot...")

        # Load the CN trees from the database we just created
        cn_trees = load_copynumber_trees_from_db(args.output_db)

        # Create the plot
        plot_copynumber_tracks(
            cn_trees=cn_trees,
            output_figure=args.plot_output,
            reference_fai=args.reference_fai,
            title=args.plot_title,
        )


if __name__ == "__main__":
    main()
