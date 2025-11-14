"""
Minimizer-based approximate string matching for fast similarity search.

This module provides a fast alternative to brute-force string matching by using
minimizers (smallest k-mers in windows) to quickly identify potential matches
before performing detailed validation.
"""

from collections import defaultdict

import numpy as np
from Bio.Seq import Seq


class MinimizerIndex:
    """
    A minimizer-based index for fast approximate string matching.

    This class pre-processes a reference string into minimizers and provides
    fast lookup methods for finding similar regions.
    """

    def __init__(self, reference: str, k: int = 15, w: int = 10):
        """
        Initialize the minimizer index.

        Args:
            reference: The reference string to index
            k: k-mer size for minimizers (default 15)
            w: window size for minimizer selection (default 10)
        """
        # check if reference is a non-empty string
        if not isinstance(reference, str) or not reference:
            raise ValueError("Reference must be a non-empty string")
        self.reference = reference.upper()  # Ensure reference is uppercase
        self.k = k
        self.w = w
        self.ref_len = len(reference)

        # Pre-compute minimizers and their positions
        self.minimizers, self.positions = self._compute_minimizers()

        # Build position index: minimizer hash -> list of positions
        self.minimizer_index = defaultdict(list)
        for pos, minimizer_hash in enumerate(self.minimizers):
            if minimizer_hash is not None:
                self.minimizer_index[minimizer_hash].append(pos)

    def _hash_kmer(self, kmer: str) -> int:
        """Hash a k-mer to an integer using a fast hash function."""
        return hash(kmer) & 0x7FFFFFFF  # Keep positive

    def _compute_minimizers(self) -> tuple[np.ndarray, np.ndarray]:
        reference_hashes = np.array([
            self._hash_kmer(self.reference[i : i + self.k])
            for i in range(len(self.reference) - self.k + 1)
        ])
        reference_revcomp = str(Seq(self.reference).reverse_complement())
        reference_hashes_revcomp = np.array([
            self._hash_kmer(reference_revcomp[i : i + self.k])
            for i in range(len(reference_revcomp) - self.k + 1)
        ])
        reference_hashes = np.minimum(reference_hashes, reference_hashes_revcomp)

        full_minimizers = np.full(len(reference_hashes), -1, dtype=int)
        full_positions = np.full(len(reference_hashes), -1, dtype=int)
        # compute minimizers by rolling over all kmers with a window of size k+w-1
        # in the rolling loop: if the next kmer is smaller than the current minimizer, update the current minimizer
        # the current mimimizer i always set to the ith position in the full_minimizers array
        # first initialize the first window
        window = reference_hashes[: self.w]
        current_minimizer = np.min(window)
        current_minimizer_position = np.argmin(window)
        full_minimizers[: self.w] = current_minimizer
        full_positions[: self.w] = current_minimizer_position
        # now iterate the rest of the reference and update the current minimizer if the next kmer is smaller or current_minimizer_position is out of the window
        for i in range(self.w, len(reference_hashes)):
            if i - current_minimizer_position >= self.w:
                current_minimizer_position = np.argmin(
                    reference_hashes[i - self.w + 1 : i + 1]
                ) + (i - self.w + 1)
                current_minimizer = reference_hashes[current_minimizer_position]
            else:
                # check if the next kmer is smaller than the current minimizer
                if reference_hashes[i] < current_minimizer:
                    current_minimizer = reference_hashes[i]
                    current_minimizer_position = i
            full_minimizers[i] = current_minimizer
            full_positions[i] = current_minimizer_position
        # then compress the minimizers to eliminate all repetitions by building a mask that is 1 when two adjacent minimizers are not equal
        mask = full_minimizers[1:] != full_minimizers[:-1]
        minimizers = full_minimizers[:-1][mask != 0]
        positions = full_positions[
            :-1
        ][
            mask != 0
        ]  # the i-th minimizer maps to the j-th position on the reference_hashes. The seq letter pos is j:j+k

        return minimizers, positions
        # now I can build the minimizer index. maps a minimizer hash to a list of positions in all minimizers.

    def reference_positions_of_minimizers(self, _minimizers: np.ndarray) -> np.ndarray:
        ref_positions = np.full(len(_minimizers), -1, dtype=int)
        for minimizer_hash in _minimizers:
            if minimizer_hash in self.minimizer_index:
                for pos in self.minimizer_index[minimizer_hash]:
                    ref_positions[pos] = self.positions[pos]
        return ref_positions

    def minimizers_of_query(self, query: str) -> np.ndarray:
        query_revcom = str(Seq(query).reverse_complement())
        hashes = np.array([
            self._hash_kmer(query[i : i + self.k])
            for i in range(len(query) - self.k + 1)
        ])
        hashes_revcomp = np.array([
            self._hash_kmer(query_revcom[i : i + self.k])
            for i in range(len(query_revcom) - self.k + 1)
        ])
        hashes = np.minimum(hashes, hashes_revcomp)
        minimizers = np.full(len(hashes), -1, dtype=int)
        # compute minimizers by rolling over all kmers with a window of size k+w-
        window = hashes[: self.w]
        current_minimizer = np.min(window)
        current_minimizer_position = np.argmin(window)
        minimizers[: self.w] = current_minimizer
        # now iterate the rest of the query and update the current minimizer if the next k
        for i in range(self.w, len(hashes)):
            if i - current_minimizer_position >= self.w:
                current_minimizer_position = np.argmin(
                    hashes[i - self.w + 1 : i + 1]
                ) + (i - self.w + 1)
                current_minimizer = hashes[current_minimizer_position]
            else:
                # check if the next kmer is smaller than the current minimizer
                if hashes[i] < current_minimizer:
                    current_minimizer = hashes[i]
                    current_minimizer_position = i
            minimizers[i] = current_minimizer
        # then compress the minimizers to eliminate all repetitions by building a mask that is 1
        mask = minimizers[1:] != minimizers[:-1]
        compressed_minimizers = minimizers[:-1][mask != 0]

        # Ensure the first minimizer is always included
        if len(compressed_minimizers) == 0 or compressed_minimizers[0] != minimizers[0]:
            compressed_minimizers = np.insert(compressed_minimizers, 0, minimizers[0])

        return compressed_minimizers

    def find_similar_regions_of_query(
        self, query: str, hits_ratio: int = 0.75, minimizer_hit_window_size: int = 5
    ) -> list[tuple[int, int]]:
        # tranforms the query to minimizers, then checks in the minimizer index for all minimizers
        query = query.upper()
        if len(query) < self.k:
            raise ValueError(
                f"Query is too short to generate minimizers. It is {len(query)} characters long, but k is {self.k}."
            )
        query_minimizers: np.ndarray = self.minimizers_of_query(query)
        if len(query_minimizers) == 0:
            query_minimizers = self._hash_kmer(query[: self.k])

        if len(query_minimizers) < minimizer_hit_window_size:
            minimizer_hit_window_size = len(query_minimizers)
        # create bool array of size self.positions (query_hits_array)
        query_hits_array = np.zeros(len(self.positions), dtype=bool)
        # for each minimizer in the query, go to the positions in the minimizer index
        # and set the value in the array to true. If the minimizer is not in the index, continue
        for minimizer in query_minimizers:
            value: list[int] | None = self.minimizer_index.get(minimizer, None)
            if value:
                for pos in value:
                    query_hits_array[pos] = True
        # now roll a window over the query_hits_array and count the number of hits in each window.
        # the window has a start and an end position. The window cannot be larger than len(query) * a factor

        seeds: list[tuple[int, int]] = []
        # new approach of checking minimizer_hit_window_size minimizer hits and calculating the fraction
        for i in range(len(query_hits_array) - minimizer_hit_window_size):
            if (
                np.sum(query_hits_array[i : i + minimizer_hit_window_size])
                >= hits_ratio * minimizer_hit_window_size
            ):
                seeds.append((
                    self.positions[i],
                    self.positions[i + minimizer_hit_window_size - 1] + self.k - 1,
                ))

        if len(seeds) == 0:
            return []

        result: list[tuple[int, int]] = []
        last_start, last_end = seeds[0]
        for start, end in seeds:
            if start > last_end:
                result.append((last_start, last_end))
                last_start = start
                last_end = end
            else:
                last_end = max(last_end, end)
        result.append((last_start, last_end))
        return result
