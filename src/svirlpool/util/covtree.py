# %%
import logging as log
import multiprocessing as mp
from pathlib import Path

import numpy as np
from intervaltree import IntervalTree
from tqdm import tqdm

from ..signalprocessing import rafs_to_coverage

log.basicConfig(level=log.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def load_reference_data(path_db: Path | str) -> np.ndarray:
    """
    Loads reference dictionary and coverage data from the database.
    data is a container that holds all parsed RAF intervals as numpy array. chr, start, end, readname
    """
    data = np.array(list(rafs_to_coverage.load_cov_from_db(path_db)))
    return data


def construct_interval_trees(
    data: np.ndarray,
) -> tuple[dict[str, IntervalTree], dict[str, list[int]]]:
    """Constructs interval trees and position sets from the input data."""
    if data.shape[0] == 0:
        log.warning("No data found to construct interval trees.")
        return {}, {}
    if data.shape[1] != 4:
        raise ValueError(
            f"Data shape is incorrect. Expected 4 columns, got {data.shape[1]} columns."
        )

    chrs = set(data[:, 0])
    intervall_trees = {}
    all_positions: dict[str, list[int]] = {}
    for chr in tqdm(chrs, desc="Building interval trees"):
        chr = str(chr)
        data_chr_intervals = data[data[:, 0] == chr][:, 1:3].astype(int)
        data_chr_readnames = data[data[:, 0] == chr][:, 3].astype(
            np.uint64
        )  # 64 bit hashed readnames
        tree = IntervalTree()
        tree = IntervalTree.from_tuples(
            zip(
                data_chr_intervals[:, 0],
                data_chr_intervals[:, 1],
                data_chr_readnames,
                strict=True,
            )
        )
        intervall_trees[chr] = tree
        all_positions[chr] = sorted(
            set(data_chr_intervals[:, 0]) | set(data_chr_intervals[:, 1])
        )
    return intervall_trees, all_positions


def process_chromosome(
    chr: str, positions: list[int], intervall_tree: IntervalTree
) -> tuple[str, IntervalTree]:
    """Processes a single chromosome to compute coverage intervals."""
    log.info(f"Processing chromosome {chr} for coverage computation")
    local_cov_tree = IntervalTree()
    for i in range(len(positions) - 1):
        start, end = positions[i], positions[i + 1]
        intervals = intervall_tree[start:end]
        reads = {i.data for i in intervals}
        coverage = len(reads)
        local_cov_tree.addi(start, end, coverage)
    return chr, local_cov_tree


def parallel_coverage_computation(
    all_positions: dict[str, list[int]],
    intervall_trees: dict[str, IntervalTree],
    num_workers: int = 4,
) -> dict[str, IntervalTree]:
    """Computes coverage in parallel for all chromosomes."""
    log.info(f"Using {num_workers} workers for parallel coverage computation")
    cov_trees = {}
    with mp.Pool(processes=num_workers) as pool:
        results = list(
            tqdm(
                pool.starmap(
                    process_chromosome,
                    [
                        (chr, all_positions[chr], intervall_trees[chr])
                        for chr in intervall_trees.keys()
                    ],
                ),
                total=len(intervall_trees),
                desc="Computing coverage",
            )
        )
    for chr, tree in results:
        cov_trees[chr] = tree
    return cov_trees


def covtree(path_db: Path | str) -> dict[str, IntervalTree]:
    """Returns a dict of interval trees. Each tree consists of intervals (start, end, readname) of the 'read alignment fragments (RAFs)'."""
    log.info(f"loading data from {path_db}")
    data = load_reference_data(path_db)
    # print the memory usage of data
    log.info(f"Data loaded. Memory usage: {data.nbytes / (1024**2):.2f} MB")
    log.info("constructing interval trees")
    intervall_trees, all_positions = construct_interval_trees(data)
    return intervall_trees
    # log.info(f"computing coverages")
    # start_time = time.time()
    # cov_trees = parallel_coverage_computation(all_positions, intervall_trees, num_workers)
    # log.info(f"Execution time: {time.time() - start_time:.2f} seconds")
    # return cov_trees


# def compute_coverage_trees(path_db:Path|str, num_workers:int) -> dict[str, IntervalTree]:
#     """Returns a dict of interval trees. Each tree consists of intervals (start, end, coverage) representing coverage. The intervals are non-overlapping."""
#     log.info(f"loading data from {path_db}")
#     data = load_reference_data(path_db)
#     log.info(f"constructing interval trees")
#     intervall_trees, all_positions = construct_interval_trees(data)
#     log.info(f"computing coverages")
#     start_time = time.time()
#     cov_trees = parallel_coverage_computation(all_positions, intervall_trees, num_workers)
#     log.info(f"Execution time: {time.time() - start_time:.2f} seconds")
#     return cov_trees
