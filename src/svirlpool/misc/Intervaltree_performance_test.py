# EXPERIMENT: given a number of intervals A and a number of intervals B, for each b in B, we want to know if there is at least one overlap with any a in A.
# We want to determine when to use an intervaltree (O(log n)) or when to check the overlaps by checking each a in A against each b in B with (in O(n**2)).
# Build a test suite that tests the range N(A)=[x**2 for x in range(1, 25)], N(B)=[x**2 for x in range(1, 25)].
# plot a heatmap of the time it takes to check the overlaps for each combination of N(A) and N(B) and print a table with the results in seconds.
# %%
import random
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from intervaltree import Interval, IntervalTree


def generate_random_intervals(
    n: int, max_pos: int = 10000, max_length: int = 100
) -> list[tuple[int, int]]:
    """Generate n random intervals with random start positions and lengths."""
    intervals = []
    for _ in range(n):
        start = random.randint(0, max_pos)
        length = random.randint(1, max_length)
        intervals.append((start, start + length))
    return intervals


def brute_force_overlap_check(
    intervals_A: list[tuple[int, int]], intervals_B: list[tuple[int, int]]
) -> list[bool]:
    """For each interval in B, check if it overlaps with any interval in A using brute force."""
    results = []
    _intervals_A = [Interval(start, end) for start, end in intervals_A]
    _intervals_B = [Interval(start, end) for start, end in intervals_B]
    for b in _intervals_B:
        has_overlap = any(a.overlaps(b) for a in _intervals_A)
        results.append(has_overlap)
    return results


def intervaltree_overlap_check(
    intervals_A: list[tuple[int, int]], intervals_B: list[tuple[int, int]]
) -> list[bool]:
    """For each interval in B, check if it overlaps with any interval in A using IntervalTree."""
    # Build IntervalTree from intervals A
    tree = IntervalTree()
    for start, end in intervals_A:
        tree[start:end] = True

    # Check overlaps for intervals B
    results = []
    for b_start, b_end in intervals_B:
        overlapping = tree[b_start:b_end]
        results.append(len(overlapping) > 0)
    return results


def benchmark_overlap_methods(experiments: int = 25):
    """Benchmark both methods across different dataset sizes."""
    # Test ranges: N(A) and N(B) = [x^2 for x in range(1, 25)]
    test_sizes = [x**2 for x in range(1, experiments)]  # 1, 4, 9, 16, ..., 576

    brute_force_times = np.zeros((len(test_sizes), len(test_sizes)))
    intervaltree_times = np.zeros((len(test_sizes), len(test_sizes)))

    print("Running interval overlap benchmark...")
    print(f"Testing {len(test_sizes)} x {len(test_sizes)} combinations")

    for i, n_a in enumerate(test_sizes):
        for j, n_b in enumerate(test_sizes):
            print(
                f"Testing N(A)={n_a}, N(B)={n_b} ({i*len(test_sizes) + j + 1}/{len(test_sizes)**2})"
            )

            # Generate random test data
            intervals_A = generate_random_intervals(n_a)
            intervals_B = generate_random_intervals(n_b)

            # Benchmark brute force method
            start_time = time.perf_counter()
            brute_results = brute_force_overlap_check(intervals_A, intervals_B)
            brute_force_times[i, j] = time.perf_counter() - start_time

            # Benchmark IntervalTree method
            start_time = time.perf_counter()
            tree_results = intervaltree_overlap_check(intervals_A, intervals_B)
            intervaltree_times[i, j] = time.perf_counter() - start_time

            # Verify results are identical
            assert (
                brute_results == tree_results
            ), f"Results differ for N(A)={n_a}, N(B)={n_b}"

    return test_sizes, brute_force_times, intervaltree_times


def plot_benchmark_results(test_sizes, brute_force_times, intervaltree_times):
    """Create heatmaps and comparison plots of the benchmark results."""

    # Create DataFrames for easier plotting
    bf_df = pd.DataFrame(brute_force_times, index=test_sizes, columns=test_sizes)
    it_df = pd.DataFrame(intervaltree_times, index=test_sizes, columns=test_sizes)

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))

    # Brute Force heatmap
    sns.heatmap(
        bf_df, ax=axes[0, 0], cmap="viridis", cbar_kws={"label": "Time (seconds)"}
    )
    axes[0, 0].set_title("Brute Force Method")
    axes[0, 0].set_xlabel("N(B) - Number of intervals in set B")
    axes[0, 0].set_ylabel("N(A) - Number of intervals in set A")

    # IntervalTree heatmap
    sns.heatmap(
        it_df, ax=axes[0, 1], cmap="viridis", cbar_kws={"label": "Time (seconds)"}
    )
    axes[0, 1].set_title("IntervalTree Method")
    axes[0, 1].set_xlabel("N(B) - Number of intervals in set B")
    axes[0, 1].set_ylabel("N(A) - Number of intervals in set A")

    # Ratio heatmap (brute_force / intervaltree)
    ratio_df = bf_df / it_df
    sns.heatmap(
        ratio_df, ax=axes[1, 0], cmap="RdYlBu_r", cbar_kws={"label": "Speedup Factor"}
    )
    axes[1, 0].set_title("Speedup Factor (BruteForce/IntervalTree)")
    axes[1, 0].set_xlabel("N(B) - Number of intervals in set B")
    axes[1, 0].set_ylabel("N(A) - Number of intervals in set A")

    # Line plot comparing methods for diagonal (N(A) = N(B))
    diagonal_bf = np.diag(brute_force_times)
    diagonal_it = np.diag(intervaltree_times)
    axes[1, 1].loglog(test_sizes, diagonal_bf, "o-", label="Brute Force", linewidth=2)
    axes[1, 1].loglog(test_sizes, diagonal_it, "s-", label="IntervalTree", linewidth=2)
    axes[1, 1].set_xlabel("Dataset Size N (N(A) = N(B))")
    axes[1, 1].set_ylabel("Time (seconds)")
    axes[1, 1].set_title("Performance Comparison (N(A) = N(B))")
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(
        "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/figures/interval_tree_tests.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()

    # Print summary table
    print("\n" + "=" * 80)
    print("BENCHMARK SUMMARY")
    print("=" * 80)

    # Find crossover point where IntervalTree becomes faster
    crossover_found = False
    for i in range(len(test_sizes)):
        for j in range(len(test_sizes)):
            if (
                brute_force_times[i, j] > intervaltree_times[i, j]
                and not crossover_found
            ):
                print(
                    f"IntervalTree becomes faster at approximately N(A)={test_sizes[i]}, N(B)={test_sizes[j]}"
                )
                crossover_found = True
                break
        if crossover_found:
            break

    # Print some key statistics
    max_speedup = np.max(ratio_df.values)
    max_idx = np.unravel_index(np.argmax(ratio_df.values), ratio_df.shape)
    print(
        f"Maximum speedup: {max_speedup:.2f}x at N(A)={test_sizes[max_idx[0]]}, N(B)={test_sizes[max_idx[1]]}"
    )

    print("\nDiagonal comparison (N(A) = N(B)):")
    print(f"{'N':<6} {'BruteForce':<12} {'IntervalTree':<12} {'Speedup':<8}")
    print("-" * 40)
    for i, n in enumerate(test_sizes[::3]):  # Show every 3rd point to avoid clutter
        idx = i * 3
        if idx < len(test_sizes):
            speedup = diagonal_bf[idx] / diagonal_it[idx]
            print(
                f"{n:<6} {diagonal_bf[idx]:<12.6f} {diagonal_it[idx]:<12.6f} {speedup:<8.2f}x"
            )


# Run the experiment

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

# Run benchmark
test_sizes, brute_force_times, intervaltree_times = benchmark_overlap_methods(30)

# Plot results
plot_benchmark_results(test_sizes, brute_force_times, intervaltree_times)
# %%
