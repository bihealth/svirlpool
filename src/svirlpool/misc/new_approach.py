# %%
import logging as log
import math
from collections import defaultdict
from pathlib import Path

# %%
import matplotlib.pyplot as plt

# 1) coverage trees
# 2) identify SV regions
# 3) select representative reads based on N clusters (N estimated based on the local coverage)
# 4) create assembly sequences
# 5) re-align reads to asssemblies
import numpy as np
from intervaltree import IntervalTree
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
from tqdm import tqdm

from ..scripts import datatypes, util


def fast_plot_coverage_intervals(
    M, output_path, pixels_per_bp=1, max_width=10000, dpi=300
):
    """
    Fast pixel-based version of the coverage plot using imshow.
    Now with chromosome labels and black background.
    """

    # Organize by chromosome
    chr_data = {}
    chr_order = []
    for row in M:
        chr_, start, end, cov = row
        if chr_ not in chr_data:
            chr_data[chr_] = []
            chr_order.append(chr_)
        chr_data[chr_].append((int(start), int(end), int(cov)))

    # Compute average coverage and sizes
    chr_avg_cov = {}
    chr_lengths = {}
    for chr_ in chr_order:
        intervals = chr_data[chr_]
        total_cov = sum((end - start) * cov for start, end, cov in intervals)
        total_len = intervals[-1][1] - intervals[0][0] if intervals else 0
        avg_cov = total_cov / total_len if total_len > 0 else 0
        chr_avg_cov[chr_] = avg_cov
        chr_lengths[chr_] = max(end for _, end, _ in intervals)

    # Compute total image height
    chr_y_offsets = {}
    y_offset = 0
    chr_heights = {}
    for chr_ in chr_order:
        chr_y_offsets[chr_] = y_offset
        height = int(4 * chr_avg_cov[chr_]) + 1
        chr_heights[chr_] = height
        y_offset += height

    total_height = y_offset

    # Determine scaling
    max_chr_len = max(chr_lengths.values())
    scale = max(1, max_chr_len // max_width)  # number of bp per pixel in x
    width = max_chr_len // scale + 1

    # Create image canvas (float32 array for log2FC values)
    img = np.zeros((total_height, width), dtype=np.float32)

    # Fill in image
    for chr_ in chr_order:
        y0 = chr_y_offsets[chr_]
        avg_cov = chr_avg_cov[chr_]
        for start, end, cov in tqdm(chr_data[chr_]):
            if cov == 0 or avg_cov == 0:
                continue
            log2_fc = math.log2(cov / avg_cov)
            log2_fc = np.clip(log2_fc, 0, 4)
            x_start = start // scale
            x_end = end // scale + 1
            height = min(cov, chr_heights[chr_])  # cap box height

            # Fill pixels: from y0 to y0 + height
            img[y0 : y0 + height, x_start:x_end] = log2_fc

    # Plot image
    fig, ax = plt.subplots(figsize=(width / 100, total_height / 100), dpi=dpi)
    cmap = get_cmap("viridis")
    norm = Normalize(vmin=0, vmax=4)
    ax.imshow(
        img,
        aspect="auto",
        cmap=cmap,
        origin="lower",
        interpolation="nearest",
        norm=norm,
    )

    # Background and styling
    ax.set_facecolor("black")
    fig.patch.set_facecolor("black")
    ax.axis("off")

    # Add chromosome labels
    for chr_ in chr_order:
        y0 = chr_y_offsets[chr_]
        h = chr_heights[chr_]
        ax.text(
            x=5,
            y=y0 + h // 2,
            s=chr_,
            color="white",
            va="center",
            ha="left",
            fontsize=8,
            fontweight="bold",
            backgroundcolor="black",
        )

    plt.tight_layout()
    plt.savefig(
        output_path, dpi=dpi, facecolor="black", bbox_inches="tight", pad_inches=0
    )
    plt.close()


def rafs_to_coverage(input: Path) -> np.ndarray:
    """Iterates all rafs and returns a numpy array of intervals."""
    intervals = []
    for raf in util.yield_from_raf(input):
        intervals.append(
            [
                raf.reference_name,
                raf.reference_alignment_start,
                raf.reference_alignment_end,
            ]
        )
    return np.array(intervals, dtype=object)


def compute_coverage_intervals(N: np.ndarray) -> np.ndarray:
    """
    Compute coverage matrix M from interval matrix N.
    M has columns (chr, start, end, cov), where cov is constant across each interval.
    """
    log.info("Computing coverage intervals from N")
    # Step 1: Generate coverage events
    events = defaultdict(list)  # chr -> list of (pos, delta)

    for row in N:
        chr_, start, end = row
        events[chr_].append((start, +1))
        events[chr_].append((end, -1))

    # Step 2: Process events and build M
    M_rows = []

    for chr_ in tqdm(events):
        chr_events = sorted(events[chr_])
        coverage = 0
        prev_pos = None

        for pos, delta in chr_events:
            if prev_pos is not None and pos != prev_pos and coverage > 0:
                M_rows.append((chr_, prev_pos, pos, coverage))
            coverage += delta
            prev_pos = pos

    M = np.array(M_rows, dtype=object)
    return M


def write_M_to_bed(M: np.ndarray, output_bed_path: Path | str):
    """
    Write coverage matrix M to a BED file.
    """
    log.info(f"Writing coverage matrix M to BED file: {output_bed_path}")
    with open(output_bed_path, "w") as bed_file:
        for row in tqdm(M):
            chr_, start, end, cov = row
            bed_file.write(f"{chr_}\t{start}\t{end}\t{cov}\n")


def load_coverage_bed(bed_path: Path | str) -> np.ndarray:
    """
    Load coverage BED file into a numpy array.
    """
    log.info(f"Loading coverage BED file: {bed_path}")
    intervals = []
    with open(bed_path, "r") as bed_file:
        for line in tqdm(bed_file):
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            chr_ = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            cov = int(parts[3])
            intervals.append([chr_, start, end, cov])
    return np.array(intervals, dtype=object)


def construct_interval_trees(data: np.ndarray) -> dict[str, IntervalTree]:
    """Constructs interval trees and position sets from the input data."""
    chrs = set(data[:, 0])
    intervall_trees = {}
    for chr in tqdm(chrs, desc="Building interval trees"):
        data_chr_intervals = data[data[:, 0] == chr][:, 1:3].astype(int)
        tree = IntervalTree()
        tree = IntervalTree.from_tuples(
            zip(data_chr_intervals[:, 0], data_chr_intervals[:, 1])
        )
        intervall_trees[chr] = tree
    return intervall_trees


def signals_matrix_from_rafs(rafs: list[datatypes.ReadAlignmentFragment]) -> np.ndarray:
    L = []
    for raf in tqdm(rafs, desc="Converting RAFs to signals matrix"):
        readhash = hash(raf.read_name)
        for signal in raf.SV_signals:
            L.append(
                (
                    raf.reference_name,
                    signal.ref_start + (signal.ref_start + signal.ref_end) / 2.0,
                    signal.sv_type,
                    signal.size,
                    readhash,
                )
            )
    return np.array(L, dtype=object)


# %%

wdir = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/test/test_new"
)
rafs_path = wdir / "rafs.filtered_by_excessive_indels.tsv.gz"
output = wdir / "coverage.bed"

# %%
N: np.ndarray = rafs_to_coverage(rafs_path)
M: np.ndarray = compute_coverage_intervals(N)
write_M_to_bed(M, output)

# %% COVERAGES
M: np.ndarray = load_coverage_bed(output)
fast_plot_coverage_intervals(M, output.with_suffix(".png"), max_width=10000, dpi=300)

# %% INTERVAL TREES
covtrees = construct_interval_trees(data=M)
# %%
rafs = list(util.yield_from_raf(rafs_path))

# %%
signals_matrix = signals_matrix_from_rafs(rafs)
signals_trees = signaltrees(signals_matrix)

# %%


def signaltrees(signals_matrix: np.ndarray) -> dict[str, IntervalTree]:
    """
    Constructs interval trees for SV signals from the signals matrix.
    The signals matrix is expected to have columns: chr, pos, sv_type, size, readhash.
    """
    log.info("Constructing signal trees from signals matrix")
    chrs = set(signals_matrix[:, 0])
    signal_trees = {}
    for chr in tqdm(chrs, desc="Building signal trees"):
        data_chr_signals = signals_matrix[signals_matrix[:, 0] == chr][:, 1:3].astype(
            float
        )
        tree = IntervalTree()
        tree = IntervalTree.from_tuples(
            zip(data_chr_signals[:, 0], data_chr_signals[:, 2])
        )
        signal_trees[chr] = tree
    return signal_trees


# %%
# identify regions
# iterate all signals with a window size of W=20_000 bp
# get the median coverage in the window
# if more then 25% of readnames have at least one sv signal in the window, then this is a candidate region
def get_candidate_regions(
    signals_matrix: np.ndarray,
    covtrees: dict[str, IntervalTree],
    window_size: int = 20_000,
    min_signal_fraction: float = 0.25,
) -> list[tuple[str, int, int]]:
    signals_matrix.sort(axis=0, order=[0, 1])  # Sort by chromosome and position
    candidate_regions = []
    for signal in tqdm(signals_matrix, desc="Identifying candidate regions"):
        chr_, pos, sv_type, size, readhash = signal
        # Get the coverage in the window
        start = pos - window_size // 2
        end = pos + window_size // 2
        coverage = covtrees[chr_].count_overlapping(start, end)
        if coverage > 0 and coverage / len(rafs) > min_signal_fraction:
            candidate_regions.append((chr_, start, end))
    return candidate_regions
