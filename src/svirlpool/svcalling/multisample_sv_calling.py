# this script reads the svPrimitives of each provided sample's file and calls SVs from it
# %%

import argparse
import gzip
import json
import logging
import os
import pickle
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from itertools import groupby
from math import ceil, floor
from pathlib import Path
from shlex import split

import attrs
import numpy as np
from intervaltree import Interval, IntervalTree
from pandas import read_csv
from scipy.stats import binom
from tqdm import tqdm

from ..localassembly import SVpatterns, svirltile
from ..util.covtree import covtree
from ..util.datastructures import UnionFind
from ..util.util import kmer_similarity_of_groups
from .SVcomposite import SVcomposite

log = logging.getLogger(__name__)

SUPPORTED_SV_TYPES: frozenset[type[SVpatterns.SVpatternType]] = frozenset({
    SVpatterns.SVpatternInversionDeletion,
    SVpatterns.SVpatternInversionDuplication,
    SVpatterns.SVpatternInsertion,
    SVpatterns.SVpatternDeletion,
    SVpatterns.SVpatternInversion,
    SVpatterns.SVpatternSingleBreakend,
    SVpatterns.SVpatternAdjacency,
})

SUPPORTED_SV_TYPE_STRINGS: frozenset[str] = frozenset({
    pattern_type.get_sv_type() for pattern_type in SUPPORTED_SV_TYPES
})
# Add a sorted list constant for consistent ordering and help text
SUPPORTED_SV_TYPE_STRINGS_LIST: list[str] = sorted(SUPPORTED_SV_TYPE_STRINGS)

SUPPORTED_SV_TYPE_STRINGS_INVERSE: dict[str, type[SVpatterns.SVpatternType]] = {
    pattern_type.get_sv_type(): pattern_type for pattern_type in SUPPORTED_SV_TYPES
}

# %% FUNCTIONS


def parse_candidate_regions_file(filepath: Path) -> dict[str, set[int]]:
    """
    Parse a TSV file containing sample names and candidate region IDs.

    Args:
        filepath: Path to TSV file with format: samplename\tcrID1,crID2,crID3

    Returns:
        dict mapping samplename to set of candidate region IDs (crIDs)
    """
    result: dict[str, set[int]] = {}
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                log.warning(
                    f"Skipping malformed line in candidate regions file: {line}"
                )
                continue
            samplename = parts[0]
            crIDs = set(map(int, parts[1].split(",")))
            result[samplename] = crIDs
    log.info(f"Parsed candidate regions for {len(result)} samples")
    return result


def cohens_d(x: list | np.ndarray, y: list | np.ndarray) -> float:
    """
    Calculate Cohen's d effect size between two samples.
    Cohen's d = (mean1 - mean2) / pooled_standard_deviation
    Where pooled_standard_deviation = sqrt(((n1-1)*s1² + (n2-1)*s2²) / (n1+n2-2))
    """
    if len(x) == 0 or len(y) == 0:
        raise ValueError("Both samples must contain at least one value")

    nx = len(x)
    ny = len(y)

    x_arr = np.array(x) if not isinstance(x, np.ndarray) else x
    y_arr = np.array(y) if not isinstance(y, np.ndarray) else y

    # Calculate means
    mean_x = np.mean(x_arr)
    mean_y = np.mean(y_arr)

    # Handle edge cases
    if nx == 1 and ny == 1:
        # Can't calculate pooled standard deviation with only one observation each
        # Return a large effect size if means differ, 0 if they're the same
        return float("inf") if mean_x != mean_y else 0.0

    # Calculate sample standard deviations (with Bessel's correction, ddof=1)
    std_x = np.std(x_arr, ddof=1) if nx > 1 else 0.0
    std_y = np.std(y_arr, ddof=1) if ny > 1 else 0.0

    # Calculate pooled standard deviation
    pooled_var = ((nx - 1) * std_x**2 + (ny - 1) * std_y**2) / (nx + ny - 2)
    pooled_std = np.sqrt(pooled_var)

    # Handle case where pooled standard deviation is 0 (all values identical)
    if pooled_std == 0:
        return 0.0 if mean_x == mean_y else float("inf")

    # Calculate Cohen's d
    cohens_d_value = (mean_x - mean_y) / pooled_std

    return float(cohens_d_value)


def get_groups_by_reference_overlaps(
    svPatterns: list[SVpatterns.SVpatternType],
    idx_mask: tuple[np.ndarray, dict[str, int]],
    used_indices: list[int],
) -> list[list[int]]:
    """Groups SVpatterns by reference overlaps of the underlying SVprimitives' alignments. Returns lists of indices in the input list."""

    mask = idx_mask[0]  # binary mask of SVpattern indices

    # Create interval trees for each chromosome
    chr_intervals: dict[str, IntervalTree] = {}

    log.info(f"Building interval trees for {len(svPatterns)} SVpatterns.")
    for svpattern_idx, svPattern in tqdm(enumerate(svPatterns)):
        if not any(mask[svpattern_idx, idx] for idx in used_indices):
            continue
        # Get all SVprimitives from this pattern
        for svprimitive in svPattern.SVprimitives:
            # Get chromosome and reference positions
            chr_name = svprimitive.consensus_aln_interval[0]
            ref_start = min(
                svprimitive.consensus_aln_interval[1],
                svprimitive.consensus_aln_interval[2],
            )
            ref_end = max(
                svprimitive.consensus_aln_interval[1],
                svprimitive.consensus_aln_interval[2],
            )

            # Ensure we have an interval tree for this chromosome
            if chr_name not in chr_intervals:
                chr_intervals[chr_name] = IntervalTree()

            # Add this interval to the tree with the SVpattern index as data
            chr_intervals[chr_name].addi(ref_start, ref_end, [svpattern_idx])

    log.info("Merging overlapping intervals.")
    for tree in chr_intervals.values():
        tree.merge_overlaps(data_reducer=lambda x, y: x + y)

    # Map each SVpattern index to its corresponding merged intervals
    merged_groups = []
    for tree in chr_intervals.values():
        for interval in tree:
            merged_groups.append(list(set(interval.data)))

    return merged_groups


def binary_svpattern_index_mask(
    svPatterns: list[SVpatterns.SVpatternType],
) -> tuple[np.ndarray, dict[str, int]]:
    """Create a binary mask for the SVpattern indices. Each index is a row and for each SV type there is a column.
    Returns a tuple of (mask, sv_type_to_index) where mask is a 2D numpy array and sv_type_to_index is a dict mapping SV types to column indices.
    """

    sv_type_to_index = {}
    sv_types = set()

    for svPattern in svPatterns:
        sv_types.add(svPattern.get_sv_type())

    sv_types = sorted(sv_types)  # Sort to ensure consistent ordering
    sv_type_to_index = {sv_type: idx for idx, sv_type in enumerate(sv_types)}

    mask = np.zeros((len(svPatterns), len(sv_types)), dtype=bool)

    for i, svPattern in enumerate(svPatterns):
        mask[i, sv_type_to_index[svPattern.get_sv_type()]] = True

    return mask, sv_type_to_index


# def load_covtrees(input:list[str]) -> dict[str,dict[str,IntervalTree]]:
#     covtrees = {}
#     for path in input:
#         samplename:str = svirltile.get_metadata(path)["samplename"]
#         covtrees[samplename] = covtree(path_db=path)
#     return covtrees


# the horizontal merge allows to merge inter with intra alignment SVpatterns
# this really only concerns insertions and deletions, which can be split across multiple aligned fragments
def svPatterns_to_horizontally_merged_svComposites(
    svPatterns: list[SVpatterns.SVpatternType],
    sv_types: set[type[SVpatterns.SVpatternType]],
) -> list[SVcomposite]:
    
    result: list[SVcomposite] = []
    if len(svPatterns) == 0:
        log.warning(f"No SVpatterns found in {svPatterns}. Returning empty list.")
        return result
    # split the insertiosn and deletions. Add the rest to results immediately.
    # filter svPatterns for supported types
    _svPatterns = [svp for svp in svPatterns if type(svp) in sv_types]
    unsupported = [svp for svp in svPatterns if type(svp) not in sv_types]
    if unsupported:
        log.warning(f"Skipping {len(unsupported)} unsupported SVpatterns during horizontal merging: {[type(svp).__name__ for svp in unsupported]}")
        unsupported.clear()
    indels = [svp for svp in _svPatterns if type(svp) in (SVpatterns.SVpatternInsertion, SVpatterns.SVpatternDeletion)]
    others = [svp for svp in _svPatterns if type(svp) not in (SVpatterns.SVpatternInsertion, SVpatterns.SVpatternDeletion)]
    _svPatterns.clear()

    # convert other svPatterns directly to svComposites
    for svp in others:
        result.append(SVcomposite.from_SVpattern(svp))

    # sort svPatterns by consensusID and read_start
    indels = sorted(indels, key=lambda x: (x.consensusID, x.read_start))

    groups = groupby(
        indels, key=lambda x: x.consensusID
    )

    # loop svPatterns of each group and connect them if they share at least one repeatID
    for (_consensusID, group) in groups:
        group = list(group)
        if len(group) == 1:
            result.append(SVcomposite.from_SVpatterns(group))
            continue

        # Create a union-find structure for this group
        uf_group = UnionFind(range(len(group)))

        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                if group[i].repeatIDs.intersection(group[j].repeatIDs):
                    # TODO: Edge case, where indels in duplicated overlapping aligned fragments are concatenated horizontally
                    uf_group.union(i, j)

        # Collect connected components
        connected_components = uf_group.get_connected_components(allow_singletons=True)

        for component in connected_components:
            sv_patterns_in_component = [group[idx] for idx in component]
            if len(sv_patterns_in_component) > 0:
                result.append(SVcomposite.from_SVpatterns(sv_patterns_in_component))

    return result


def can_merge_svComposites_insertions(
    a: SVcomposite,
    b: SVcomposite,
    apriori_size_difference_fraction_tolerance: float,
    near: int,
    d: float = 2.0,
    min_kmer_overlap: float = 0.7,
    verbose: bool = False,
) -> bool:
    # throroughly check the input svComposites
    # 1) test if they have svPatterns
    if (
        not a.svPatterns
        or not b.svPatterns
        or len(a.svPatterns) == 0
        or len(b.svPatterns) == 0
    ):
        raise ValueError(
            "can_merge_svComposites_insertions: SVcomposites must contain at least one SVpattern to merge."
        )
    # 2) check every svPattern if it contains svPrimitives
    for svp in a.svPatterns + b.svPatterns:
        if not svp.SVprimitives or len(svp.SVprimitives) == 0:
            raise ValueError(
                "can_merge_svComposites_insertions: SVpattern must contain at least one SVprimitive to merge."
            )
    # 3) test every svPattern if it has a genotypeMeasurement with supporting_reads_start
    for svPattern in a.svPatterns + b.svPatterns:
        for svPrimitive in svPattern.SVprimitives:
            if (
                not svPrimitive.genotypeMeasurement
                or not svPrimitive.genotypeMeasurement.supporting_reads_start
            ):
                raise ValueError(
                    "can_merge_svComposites_insertions: SVprimitive must have a genotypeMeasurement with supporting_reads_start to merge."
                )

    are_near: bool = a.overlaps_any(
        b, tolerance_radius=max(near, max(a.get_size(), b.get_size()))
    )  # either spatially close or share at least one repeatID
    if not are_near:
        if verbose:
            print(
                f"Cannot merge insertions: SVcomposites are not near (tolerance_radius={near}) ({a.get_regions()} vs {b.get_regions()})"
            )
        return False

    # - have similar sizes (computes from a combined population of both composites), tested with Cohen's D
    # get population of sizes

    size_tolerance_a = sizetolerance_from_SVcomposite(a)
    size_a = abs(a.get_size())
    size_tolerance_b = sizetolerance_from_SVcomposite(b)
    size_b = abs(b.get_size())

    population_a: np.ndarray = (
        np.array(a.get_size_populations(), dtype=np.int32) + a.get_size()
    )
    if len(population_a) == 0:
        raise ValueError(
            f"can_merge_svComposites_insertions: SVcomposite 'a' has no size population to compare. SVcomposite: {a}"
        )
    population_b: np.ndarray = (
        np.array(b.get_size_populations(), dtype=np.int32) + b.get_size()
    )
    if len(population_b) == 0:
        raise ValueError(
            f"can_merge_svComposites_insertions: SVcomposite 'b' has no size population to compare. SVcomposite: {b}"
        )

    # Calculate original means
    mean_a = np.mean(population_a)
    mean_b = np.mean(population_b)
    mean_diff = abs(mean_a - mean_b)

    # Check if means are within tolerance range (trivial case)
    if mean_diff <= (size_tolerance_a + size_tolerance_b):
        # Means are close enough - consider them similar
        similar_size = True
        if verbose:
            print(
                f"Size comparison - Trivial case: mean_diff={mean_diff:.1f} <= tolerance_sum={size_tolerance_a + size_tolerance_b:.1f}"
            )
            print(f"  Mean A: {mean_a:.1f}, Mean B: {mean_b:.1f}")
            print(
                f"  Tolerance A: {size_tolerance_a:.1f}, Tolerance B: {size_tolerance_b:.1f}"
            )
    else:
        # Shift means towards each other and calculate Cohen's d
        if mean_a > mean_b:
            shifted_mean_a = mean_a - size_tolerance_a
            shifted_mean_b = mean_b + size_tolerance_b
        else:
            shifted_mean_a = mean_a + size_tolerance_a
            shifted_mean_b = mean_b - size_tolerance_b

        # Create shifted populations by adjusting all values by the shift amount
        shift_a = shifted_mean_a - mean_a
        shift_b = shifted_mean_b - mean_b
        shifted_population_a = population_a + shift_a
        shifted_population_b = population_b + shift_b

        # Calculate Cohen's d with shifted populations
        cohensD = cohens_d(shifted_population_a, shifted_population_b)

        similar_size = abs(cohensD) <= abs(d) or abs(
            size_a - size_b
        ) <= apriori_size_difference_fraction_tolerance * max(size_a, size_b)
        # or abs(size_a - size_b) <= int(ceil(np.log2((size_a + size_b)*0.5))) \

        if verbose:
            print(
                f"Size comparison - Shifted Cohen's D: {cohensD:.3f}, threshold: {d}, similar_size: {similar_size}"
            )
            print(
                f"  Original means - A: {mean_a:.1f}, B: {mean_b:.1f}, diff: {mean_diff:.1f}"
            )
            print(f"  Shifted means - A: {shifted_mean_a:.1f}, B: {shifted_mean_b:.1f}")
            print(
                f"  Tolerances - A: {size_tolerance_a:.1f}, B: {size_tolerance_b:.1f}"
            )
            print(f"  Population A sizes: {population_a.tolist()}")
            print(f"  Population B sizes: {population_b.tolist()}")

    # if the inserton is between two break points, then it is necessary to take the sequence from inserted_sequence
    if size_a < 100 and size_b < 100:
        similarity = 1.0  # don't compare k-mers for small deletions, they are too short to be meaningful
    else:
        similarity: float = kmer_similarity_of_groups(
            group_a=a.get_inserted_sequences(), group_b=b.get_inserted_sequences()
        )

    similar_insertion_kmers: bool = similarity >= min_kmer_overlap

    if verbose:
        print(
            f"K-mer similarity: {similarity:.3f}, threshold: {min_kmer_overlap}, similar_kmers: {similar_insertion_kmers}"
        )

    if not similar_size:
        if verbose:
            print("Cannot merge insertions: sizes are not similar enough")
        return False
    if not similar_insertion_kmers:
        if verbose:
            print("Cannot merge insertions: k-mer similarity is too low")
        return False

    if verbose:
        print("Can merge insertions: all criteria passed")
    return True


def sizetolerance_from_SVcomposite(a: SVcomposite) -> float:
    # if a is sv type insertion or inversion, get inserted complexity tracks
    mean_complexity: float = 1.0
    if issubclass(a.sv_type, SVpatterns.SVpatternInsertion) or issubclass(
        a.sv_type, SVpatterns.SVpatternInversion
    ):
        complexities: list[np.ndarray] = a.get_inserted_complexity_tracks()
        if (
            complexities is None
            or len(complexities) == 0
            or np.sum((np.sum(s) for s in complexities)) == 0
        ):
            mean_complexity = 0.0
        else:
            mean_complexity = float(
                np.average(
                    [np.mean(c) for c in complexities],
                    weights=[len(c) for c in complexities],
                )
            )
    elif issubclass(a.sv_type, SVpatterns.SVpatternDeletion):
        complexities: list[np.ndarray] = a.get_reference_complexity_tracks()
        if (
            complexities is None
            or len(complexities) == 0
            or np.sum((np.sum(s) for s in complexities)) == 0
        ):
            mean_complexity = 0.0
        else:
            mean_complexity = float(
                np.average(
                    [np.mean(c) for c in complexities],
                    weights=[len(c) for c in complexities],
                )
            )
    size_tolerance = (1.0 - mean_complexity) * float(abs(a.get_size()))
    return size_tolerance


def can_merge_svComposites_deletions(
    a: SVcomposite,
    b: SVcomposite,
    apriori_size_difference_fraction_tolerance: float,
    near: int,
    d: float = 2.0,
    min_kmer_overlap: float = 0.7,
    verbose: bool = False,
) -> bool:
    # throroughly check the input svComposites
    # 1) test if they have svPatterns
    if (
        not a.svPatterns
        or not b.svPatterns
        or len(a.svPatterns) == 0
        or len(b.svPatterns) == 0
    ):
        raise ValueError(
            f"can_merge_svComposites_deletions: SVcomposites must contain at least one SVpattern to merge. SVcomposites: {a}, {b}"
        )
    # 2) check every svPattern if it contains svPrimitives
    for svp in a.svPatterns + b.svPatterns:
        if not svp.SVprimitives or len(svp.SVprimitives) == 0:
            raise ValueError(
                f"can_merge_svComposites_deletions: SVpattern must contain at least one SVprimitive to merge. SVpattern: {svp}."
            )
    # 3) test every svPattern if it has a genotypeMeasurement with supporting_reads_start
    for svPattern in a.svPatterns + b.svPatterns:
        for svPrimitive in svPattern.SVprimitives:
            if (
                not svPrimitive.genotypeMeasurement
                or not svPrimitive.genotypeMeasurement.supporting_reads_start
            ):
                raise ValueError(
                    f"can_merge_svComposites_deletions: SVprimitive must have a genotypeMeasurement with supporting_reads_start to merge svPrimitive: {svPrimitive}\nsvPattern: {svPattern}"
                )

    are_near: bool = a.overlaps_any(
        b, tolerance_radius=near
    )  # eiter spatially close or share at least one repeatID
    if not are_near:
        if verbose:
            print(
                f"Cannot merge deletions: SVcomposites are not near (tolerance_radius={near}) ({a.get_regions()} vs {b.get_regions()})"
            )
        return False

    size_tolerance_a = sizetolerance_from_SVcomposite(a)
    size_a = abs(a.get_size())
    size_tolerance_b = sizetolerance_from_SVcomposite(b)
    size_b = abs(b.get_size())

    population_a: np.ndarray = (
        np.array(a.get_size_populations(), dtype=np.int32) + a.get_size()
    )
    if len(population_a) == 0:
        raise ValueError(
            f"can_merge_svComposites_deletions: SVcomposite 'a' has no size population to compare. SVcomposite: {a}"
        )
    population_b: np.ndarray = (
        np.array(b.get_size_populations(), dtype=np.int32) + b.get_size()
    )
    if len(population_b) == 0:
        raise ValueError(
            f"can_merge_svComposites_deletions: SVcomposite 'b' has no size population to compare. SVcomposite: {b}"
        )

    # Calculate original means
    mean_a = np.mean(population_a)
    mean_b = np.mean(population_b)
    mean_diff = abs(mean_a - mean_b)

    # Check if means are within tolerance range (trivial case)
    if mean_diff <= (size_tolerance_a + size_tolerance_b):
        # Means are close enough - consider them similar
        similar_size = True
        if verbose:
            print(
                f"Size comparison - Trivial case: mean_diff={mean_diff:.1f} <= tolerance_sum={size_tolerance_a + size_tolerance_b:.1f}"
            )
            print(f"  Mean A: {mean_a:.1f}, Mean B: {mean_b:.1f}")
            print(
                f"  Tolerance A: {size_tolerance_a:.1f}, Tolerance B: {size_tolerance_b:.1f}"
            )
    else:
        # Shift means towards each other and calculate Cohen's d
        if mean_a > mean_b:
            shifted_mean_a = mean_a - size_tolerance_a
            shifted_mean_b = mean_b + size_tolerance_b
        else:
            shifted_mean_a = mean_a + size_tolerance_a
            shifted_mean_b = mean_b - size_tolerance_b

        # Create shifted populations by adjusting all values by the shift amount
        shift_a = shifted_mean_a - mean_a
        shift_b = shifted_mean_b - mean_b
        shifted_population_a = population_a + shift_a
        shifted_population_b = population_b + shift_b

        # Calculate Cohen's d with shifted populations
        cohensD = cohens_d(shifted_population_a, shifted_population_b)

        similar_size = abs(cohensD) <= abs(d) or abs(
            size_a - size_b
        ) <= apriori_size_difference_fraction_tolerance * max(size_a, size_b)
        # or abs(size_a - size_b) <= int(ceil(np.log2((size_a + size_b)*0.5))) \

        if verbose:
            print(
                f"Size comparison - Shifted Cohen's D: {cohensD:.3f}, threshold: {d}, similar_size: {similar_size}"
            )
            print(
                f"  Original means - A: {mean_a:.1f}, B: {mean_b:.1f}, diff: {mean_diff:.1f}"
            )
            print(f"  Shifted means - A: {shifted_mean_a:.1f}, B: {shifted_mean_b:.1f}")
            print(
                f"  Tolerances - A: {size_tolerance_a:.1f}, B: {size_tolerance_b:.1f}"
            )
            print(f"  Population A sizes: {population_a.tolist()}")
            print(f"  Population B sizes: {population_b.tolist()}")

    if size_a < 100 and size_b < 100:
        similarity = 1.0  # don't compare k-mers for small deletions, they are too short to be meaningful
    else:
        ref_seqs_a = a.get_reference_sequences()
        ref_seqs_b = b.get_reference_sequences()
        similarity: float = kmer_similarity_of_groups(
            group_a=ref_seqs_a, group_b=ref_seqs_b
        )

    similar_deletion_kmers: bool = similarity >= min_kmer_overlap

    if verbose:
        print(
            f"K-mer similarity: {similarity:.3f}, threshold: {min_kmer_overlap}, similar_kmers: {similar_deletion_kmers}"
        )

    if not similar_size:
        if verbose:
            print("Cannot merge deletions: sizes are not similar enough")
        return False
    if not similar_deletion_kmers:
        if verbose:
            print("Cannot merge deletions: k-mer similarity is too low")

        return False

    if verbose:
        print("Can merge deletions: all criteria passed")
    return True


def vertically_merged_svComposites_from_group(
    group: list[SVcomposite],
    d: float,
    near: int,
    min_kmer_overlap: float,
    apriori_size_difference_fraction_tolerance: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    # The group is divided into subgroups with shared SV types
    insertions = {
        i: svc
        for i, svc in enumerate(group)
        if issubclass(svc.sv_type, SVpatterns.SVpatternInsertion)
    }
    deletions = {
        i: svc
        for i, svc in enumerate(group)
        if issubclass(svc.sv_type, SVpatterns.SVpatternDeletion)
    }
    inversions = {
        i: svc
        for i, svc in enumerate(group)
        if issubclass(svc.sv_type, SVpatterns.SVpatternInversion)
    }
    breakends = {
        i: svc
        for i, svc in enumerate(group)
        if issubclass(svc.sv_type, SVpatterns.SVpatternSingleBreakend)
    }
    cpx = {
        i: svc
        for i, svc in enumerate(group)
        if issubclass(svc.sv_type, SVpatterns.SVpatternComplex)
    }

    used_indices = set().union(
        insertions.keys(),
        deletions.keys(),
        inversions.keys(),
        breakends.keys(),
        cpx.keys(),
    )
    others = [
        svc for i, svc in enumerate(group) if i not in used_indices
    ]  # BND, DUP, etc.
    if verbose:
        # get the set of types of others and print ehir counts per type
        other_types = defaultdict(int)
        for svc in others:
            other_types[svc.sv_type.get_sv_type()] += 1
        print("Others svComposites counts by type:")
        for k, v in other_types.items():
            print(f"  {k}:\t{v}")

    merged_insertions = merge_insertions(
        insertions=list(insertions.values()),
        d=d,
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        verbose=verbose,
        apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
    )
    merged_deletions = merge_deletions(
        deletions=list(deletions.values()),
        d=d,
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        verbose=verbose,
        apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
    )
    merged_inversions = merge_inversions(
        inversions=list(inversions.values()),
        d=d,
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        verbose=verbose,
        apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
    )
    merged_breakends = merge_breakends(
        breakends=list(breakends.values()),
        d=d,
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        verbose=verbose,
        apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
    )
    merged_cpx = merge_complexes(
        complexes=list(cpx.values()),
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        verbose=verbose,
    )

    return (
        merged_insertions
        + merged_deletions
        + merged_inversions
        + merged_breakends
        + merged_cpx
        + others
    )


def merge_alt_count_dicts(
    dicts: list[dict[str, dict[str, int]]],
) -> dict[str, dict[str, int]]:
    merged: dict[str, dict[str, int]] = {}
    for data in dicts:
        for samplename, consensus_counts in data.items():
            target = merged.setdefault(samplename, {})
            for consensus_id, alt_count in consensus_counts.items():
                target[consensus_id] = max(target.get(consensus_id, 0), alt_count)
    return merged


def merge_insertions(
    insertions: list[SVcomposite],
    d: float,
    near: int,
    min_kmer_overlap: float,
    apriori_size_difference_fraction_tolerance: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    # 1) test if they have svPatterns
    if len(insertions) == 0:
        return []
    for svComposite in insertions:
        if not svComposite.svPatterns or len(svComposite.svPatterns) == 0:
            raise ValueError(
                "merge_insertions: SVcomposites must contain at least one SVpattern to merge."
            )
        # 2) check every svPattern if it contains svPrimitives
        for svp in svComposite.svPatterns:
            if not svp.SVprimitives or len(svp.SVprimitives) == 0:
                raise ValueError(
                    "merge_insertions: SVpattern must contain at least one SVprimitive to merge."
                )
        # 3) test every svPattern if it has a genotypeMeasurement with supporting_reads_start
        for svPattern in svComposite.svPatterns:
            for svPrimitive in svPattern.SVprimitives:
                if (
                    not svPrimitive.genotypeMeasurement
                    or not svPrimitive.genotypeMeasurement.supporting_reads_start
                ):
                    raise ValueError(
                        "merge_insertions: SVprimitive must have a genotypeMeasurement with supporting_reads_start to merge."
                    )
    if not all(
        issubclass(sv.sv_type, SVpatterns.SVpatternInsertion) for sv in insertions
    ):
        raise ValueError(
            "All SVcomposites must be of type SVpatternInsertion to merge them."
        )

    """Merge insertions that overlap on the reference and have similar sizes."""
    if len(insertions) == 0:
        return []
    uf = UnionFind(range(len(insertions)))
    for i in range(len(insertions)):
        for j in range(i + 1, len(insertions)):
            if verbose:
                print(
                    f"Checking insertions {i} and {j}: {insertions[i].svPatterns[0].samplename}:{insertions[i].svPatterns[0].consensusID} vs {insertions[j].svPatterns[0].samplename}:{insertions[j].svPatterns[0].consensusID}"
                )
                # print(f"with support: {insertions[i].get_alt_readnames_per_sample()} vs {insertions[j].get_alt_readnames_per_sample()}")
            if can_merge_svComposites_insertions(
                a=insertions[i],
                b=insertions[j],
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
                d=d,
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                verbose=verbose,
            ):
                uf.union(i, j)
            else:
                if verbose:
                    print(f"Not merging insertions {i} and {j}")
    result: list[SVcomposite] = []
    for cc in uf.get_connected_components(
        allow_singletons=True
    ):  # connected components of svComposites to merge
        result.append(
            SVcomposite.from_SVpatterns(
                svPatterns=[
                    svPattern for idx in cc for svPattern in insertions[idx].svPatterns
                ]
            )
        )
    return result


def merge_deletions(
    deletions: list[SVcomposite],
    d: float,
    near: int,
    min_kmer_overlap: float,
    apriori_size_difference_fraction_tolerance: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    """Merge deletions that overlap on the reference and have similar sizes."""
    if len(deletions) == 0:
        return []
    # check if all svPatterns are deletions
    if not all(
        issubclass(sv.sv_type, SVpatterns.SVpatternDeletion) for sv in deletions
    ):
        raise ValueError(
            "All SVcomposites must be of type SVpatternDeletion to merge them."
        )
    uf = UnionFind(range(len(deletions)))
    for i in range(len(deletions)):
        for j in range(i + 1, len(deletions)):
            if verbose:
                print(
                    f"Checking deletions {i} and {j}: {deletions[i].svPatterns[0].samplename}:{deletions[i].svPatterns[0].consensusID} vs {deletions[j].svPatterns[0].samplename}:{deletions[j].svPatterns[0].consensusID}"
                )
                # print(f"with support: {deletions[i].get_alt_readnames_per_sample()} vs {deletions[j].get_alt_readnames_per_sample()}")
            if can_merge_svComposites_deletions(
                a=deletions[i],
                b=deletions[j],
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
                d=d,
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                verbose=verbose,
            ):
                uf.union(i, j)
            else:
                if verbose:
                    print(f"Not merging deletions {i} and {j}")
    result: list[SVcomposite] = []
    for cc in uf.get_connected_components(allow_singletons=True):
        svPatterns_to_merge = [
            svPattern for idx in cc for svPattern in deletions[idx].svPatterns
        ]
        result.append(SVcomposite.from_SVpatterns(svPatterns=svPatterns_to_merge))
    return result


def can_merge_svComposites_inversions(
    a: SVcomposite,
    b: SVcomposite,
    apriori_size_difference_fraction_tolerance: float,
    near: int,
    d: float = 2.0,
    min_kmer_overlap: float = 0.7,
    verbose: bool = False,
) -> bool:
    """"""
    if (
        not a.svPatterns
        or not b.svPatterns
        or len(a.svPatterns) == 0
        or len(b.svPatterns) == 0
    ):
        raise ValueError(
            "can_merge_svComposites_inversions: SVcomposites must contain at least one SVpattern to merge."
        )
    # 2) check every svPattern if it contains svPrimitives
    for svp in a.svPatterns + b.svPatterns:
        if not svp.SVprimitives or len(svp.SVprimitives) == 0:
            raise ValueError(
                "can_merge_svComposites_inversions: SVpattern must contain at least one SVprimitive to merge."
            )
    # 3) test every svPattern if it has a genotypeMeasurement with supporting_reads_start
    for svPattern in a.svPatterns + b.svPatterns:
        for svPrimitive in svPattern.SVprimitives:
            if (
                not svPrimitive.genotypeMeasurement
                or not svPrimitive.genotypeMeasurement.supporting_reads_start
            ):
                raise ValueError(
                    "can_merge_svComposites_inversions: SVprimitive must have a genotypeMeasurement with supporting_reads_start to merge."
                )
    if issubclass(a.sv_type, SVpatterns.SVpatternInversionTranslocation) != issubclass(
        b.sv_type, SVpatterns.SVpatternInversionTranslocation
    ):
        log.debug(
            f"Cannot merge inversions: SVcomposites {a.svPatterns[0].consensusID},{b.svPatterns[0].consensusID} are of different subclasses (InversionTranslocation vs Inversion)."
        )
        return False

    are_near: bool = a.overlaps_any(
        b, tolerance_radius=near
    )  # eiter spatially close or share at least one repeatID
    if not are_near:
        log.debug(
            f"Cannot merge inversions: SVcomposites {a.svPatterns[0].consensusID},{b.svPatterns[0].consensusID} are not near (tolerance_radius={near}) ({a.get_regions()} vs {b.get_regions()})"
        )
        return False

    if (
        a.get_representative_SVpattern().SVprimitives[0].chr
        != b.get_representative_SVpattern().SVprimitives[0].chr
        or a.get_representative_SVpattern().SVprimitives[-1].chr
        != b.get_representative_SVpattern().SVprimitives[-1].chr
    ):
        log.debug(
            f"Cannot merge inversions: SVcomposites {a.svPatterns[0].consensusID},{b.svPatterns[0].consensusID} are on different chromosomes ({a.get_representative_SVpattern().SVprimitives[0].chr},{a.get_representative_SVpattern().SVprimitives[-1].chr} vs {b.get_representative_SVpattern().SVprimitives[0].chr},{b.get_representative_SVpattern().SVprimitives[-1].chr})"
        )
        return False

    # in the case that inversions of different subclasses are merged, it is necessary to measure their inner sizes to find a good size comparison,
    # since outer break ends can align very well while the inverted interval can be of very different sizes.
    # this allows to merge e.g. inverted dels with inverted dups given a stronger local noise.
    size_tolerance_a = sizetolerance_from_SVcomposite(a)
    size_a = abs(a.get_size(inner=True))
    size_tolerance_b = sizetolerance_from_SVcomposite(b)
    size_b = abs(b.get_size(inner=True))

    population_a: np.ndarray = (
        np.array(a.get_size_populations(), dtype=np.int32) + a.get_size()
    )
    if len(population_a) == 0:
        raise ValueError(
            f"can_merge_svComposites_inversions: SVcomposite 'a' has no size population to compare. SVcomposite: {a}"
        )
    population_b: np.ndarray = (
        np.array(b.get_size_populations(), dtype=np.int32) + b.get_size()
    )
    if len(population_b) == 0:
        raise ValueError(
            f"can_merge_svComposites_inversions: SVcomposite 'b' has no size population to compare. SVcomposite: {b}"
        )

    # Calculate original means
    mean_a = np.mean(population_a)
    mean_b = np.mean(population_b)
    mean_diff = abs(mean_a - mean_b)

    # Check if means are within tolerance range (trivial case)
    if mean_diff <= (size_tolerance_a + size_tolerance_b):
        # Means are close enough - consider them similar
        similar_size = True
        if verbose:
            print(
                f"Size comparison - Trivial case: mean_diff={mean_diff:.1f} <= tolerance_sum={size_tolerance_a + size_tolerance_b:.1f}"
            )
            print(f"  Mean A: {mean_a:.1f}, Mean B: {mean_b:.1f}")
            print(f"  Tolerance A: {size_tolerance_a:.1f}")
            print(f"  Tolerance B: {size_tolerance_b:.1f}")
    else:
        # Shift means towards each other and calculate Cohen's d
        if mean_a > mean_b:
            shifted_mean_a = mean_a - size_tolerance_a
            shifted_mean_b = mean_b + size_tolerance_b
        else:
            shifted_mean_a = mean_a + size_tolerance_a
            shifted_mean_b = mean_b - size_tolerance_b

        # Create shifted populations by adjusting all values by the shift amount
        shift_a = shifted_mean_a - mean_a
        shift_b = shifted_mean_b - mean_b
        shifted_population_a = population_a + shift_a
        shifted_population_b = population_b + shift_b

        # Calculate Cohen's d with shifted populations
        cohensD = cohens_d(shifted_population_a, shifted_population_b)

        similar_size = abs(cohensD) <= abs(d) or abs(
            size_a - size_b
        ) <= apriori_size_difference_fraction_tolerance * max(size_a, size_b)

        if verbose:
            print(
                f"Size comparison - Shifted Cohen's D: {cohensD:.3f}, threshold: {d}, similar_size: {similar_size}"
            )
            print(
                f"  Original means - A: {mean_a:.1f}, B: {mean_b:.1f}, diff: {mean_diff:.1f}"
            )
            print(f"  Shifted means - A: {shifted_mean_a:.1f}, B: {shifted_mean_b:.1f}")
            print(
                f"  Tolerances - A: {size_tolerance_a:.1f}, B: {size_tolerance_b:.1f}"
            )
            print(f"  Population A sizes: {population_a.tolist()}")
            print(f"  Population B sizes: {population_b.tolist()}")

    # Use inserted sequences from inversions for k-mer similarity comparison
    # For inversions, we need to get the inverted sequences (inserted_sequence from SVpatternInversion)
    inverted_sequences_a = []
    inverted_sequences_b = []

    for svPattern in a.svPatterns:
        if isinstance(svPattern, SVpatterns.SVpatternInversion):
            seq = svPattern.get_sequence()
            inverted_sequences_a.append(seq)

    for svPattern in b.svPatterns:
        if isinstance(svPattern, SVpatterns.SVpatternInversion):
            seq = svPattern.get_sequence()
            inverted_sequences_b.append(seq)

    if size_a < 100 and size_b < 100:
        similarity = 1.0  # don't compare k-mers for small inversions, they are too short to be meaningful
    else:
        similarity: float = kmer_similarity_of_groups(
            group_a=inverted_sequences_a, group_b=inverted_sequences_b
        )

    similar_inversion_kmers: bool = similarity >= min_kmer_overlap

    if verbose:
        print(
            f"K-mer similarity: {similarity:.3f}, threshold: {min_kmer_overlap}, similar_kmers: {similar_inversion_kmers}"
        )

    if not similar_size:
        if verbose:
            print("Cannot merge inversions: sizes are not similar enough")
        return False
    if not similar_inversion_kmers:
        if verbose:
            print("Cannot merge inversions: k-mer similarity is too low")
        return False

    if verbose:
        print("Can merge inversions: all criteria passed")
    return True


def merge_inversions(
    inversions: list[SVcomposite],
    d: float,
    near: int,
    min_kmer_overlap: float,
    apriori_size_difference_fraction_tolerance: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    """Merge inversions that overlap on the reference and have similar sizes."""
    # 1) test if they have svPatterns
    if len(inversions) == 0:
        log.debug("merge_inversions: No inversions to merge.")
        return []
    for svComposite in inversions:
        if not svComposite.svPatterns or len(svComposite.svPatterns) == 0:
            raise ValueError(
                "merge_inversions: SVcomposites must contain at least one SVpattern to merge."
            )
        # 2) check every svPattern if it contains svPrimitives
        for svp in svComposite.svPatterns:
            if not svp.SVprimitives or len(svp.SVprimitives) == 0:
                raise ValueError(
                    "merge_inversions: SVpattern must contain at least one SVprimitive to merge."
                )
        # 3) test every svPattern if it has a genotypeMeasurement with supporting_reads_start
        for svPattern in svComposite.svPatterns:
            for svPrimitive in svPattern.SVprimitives:
                if (
                    not svPrimitive.genotypeMeasurement
                    or not svPrimitive.genotypeMeasurement.supporting_reads_start
                ):
                    raise ValueError(
                        "merge_inversions: SVprimitive must have a genotypeMeasurement with supporting_reads_start to merge."
                    )
    if not all(
        issubclass(sv.sv_type, SVpatterns.SVpatternInversion) for sv in inversions
    ):
        raise ValueError(
            "All SVcomposites must be of type SVpatternInversion to merge them."
        )

    # check if all svPatterns are inversions
    # create a set of sv_types and check if all are subclasses of SVpatternInversion
    present_sv_types = {sv.sv_type for sv in inversions}
    if not all(
        issubclass(sv_type, SVpatterns.SVpatternInversion)
        for sv_type in present_sv_types
    ):
        raise ValueError(
            "All SVcomposites must be a subclass of SVpatternInversion to merge them."
        )
    uf = UnionFind(range(len(inversions)))
    log.debug(f"Starting merge of {len(inversions)} inversions.")
    for i in range(len(inversions)):
        for j in range(i + 1, len(inversions)):
            if can_merge_svComposites_inversions(
                a=inversions[i],
                b=inversions[j],
                d=d,
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
                verbose=verbose,
            ):
                uf.union(i, j)
                log.debug(
                    f"Connected inversions: {inversions[i].svPatterns[0].consensusID} and {inversions[j].svPatterns[0].consensusID}"
                )
    result: list[SVcomposite] = []
    for cc in uf.get_connected_components(allow_singletons=True):
        log.debug(
            f"Merging group of size {len(cc)} with consensusIDs: {[inversions[idx].svPatterns[0].consensusID for idx in cc]}"
        )
        result.append(
            SVcomposite.from_SVpatterns(
                svPatterns=[
                    svPattern for idx in cc for svPattern in inversions[idx].svPatterns
                ]
            )
        )
    return result


def can_merge_svComposites_breakends(
    a: SVcomposite,
    b: SVcomposite,
    apriori_size_difference_fraction_tolerance: float,
    near: int,
    d: float = 2.0,
    min_kmer_overlap: float = 0.7,
    verbose: bool = False,
) -> bool:
    """Check if two breakend SVcomposites can be merged based on proximity and sequence context similarity."""
    # 1) test if they have svPatterns
    if (
        not a.svPatterns
        or not b.svPatterns
        or len(a.svPatterns) == 0
        or len(b.svPatterns) == 0
    ):
        raise ValueError(
            "can_merge_svComposites_breakends: SVcomposites must contain at least one SVpattern to merge."
        )
    # 2) check every svPattern if it contains svPrimitives
    for svp in a.svPatterns + b.svPatterns:
        if not svp.SVprimitives or len(svp.SVprimitives) == 0:
            raise ValueError(
                "can_merge_svComposites_breakends: SVpattern must contain at least one SVprimitive to merge."
            )
    # 3) test every svPattern if it has a genotypeMeasurement with supporting_reads_start
    for svPattern in a.svPatterns + b.svPatterns:
        for svPrimitive in svPattern.SVprimitives:
            if (
                not svPrimitive.genotypeMeasurement
                or not svPrimitive.genotypeMeasurement.supporting_reads_start
            ):
                raise ValueError(
                    "can_merge_svComposites_breakends: SVprimitive must have a genotypeMeasurement with supporting_reads_start to merge."
                )

    are_near: bool = a.overlaps_any(
        b, tolerance_radius=near
    )  # either spatially close or share at least one repeatID
    if not are_near:
        if verbose:
            print(
                f"Cannot merge breakends: SVcomposites are not near (tolerance_radius={near}) ({a.get_regions()} vs {b.get_regions()})"
            )
        return False

    # For breakends, we compare sequence contexts using k-mer similarity
    # Get sequence contexts from both breakends
    contexts_a = (
        a.get_inserted_sequences()
    )  # calls get_sequence on the SVpattern, which routes to get_sequence_context
    contexts_b = b.get_inserted_sequences()

    if len(contexts_a) == 0 or len(contexts_b) == 0:
        if verbose:
            print(
                "Cannot merge breakends: one or both SVcomposites have no sequence context"
            )
        return False

    similarity: float = kmer_similarity_of_groups(
        group_a=contexts_a, group_b=contexts_b
    )

    similar_context_kmers: bool = similarity >= min_kmer_overlap

    if verbose:
        print(
            f"K-mer similarity: {similarity:.3f}, threshold: {min_kmer_overlap}, similar_kmers: {similar_context_kmers}"
        )

    if not similar_context_kmers:
        if verbose:
            print("Cannot merge breakends: k-mer similarity is too low")
        return False

    if verbose:
        print("Can merge breakends: all criteria passed")
    return True


def merge_breakends(
    breakends: list[SVcomposite],
    d: float,
    near: int,
    min_kmer_overlap: float,
    apriori_size_difference_fraction_tolerance: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    """Merge breakends that are spatially close and have similar sequence contexts."""
    # 1) test if they have svPatterns
    if len(breakends) == 0:
        return []
    for svComposite in breakends:
        if not svComposite.svPatterns or len(svComposite.svPatterns) == 0:
            raise ValueError(
                "merge_breakends: SVcomposites must contain at least one SVpattern to merge."
            )
        # 2) check every svPattern if it contains svPrimitives
        for svp in svComposite.svPatterns:
            if not svp.SVprimitives or len(svp.SVprimitives) == 0:
                raise ValueError(
                    "merge_breakends: SVpattern must contain at least one SVprimitive to merge."
                )
        # 3) test every svPattern if it has a genotypeMeasurement with supporting_reads_start
        for svPattern in svComposite.svPatterns:
            for svPrimitive in svPattern.SVprimitives:
                if (
                    not svPrimitive.genotypeMeasurement
                    or not svPrimitive.genotypeMeasurement.supporting_reads_start
                ):
                    raise ValueError(
                        "merge_breakends: SVprimitive must have a genotypeMeasurement with supporting_reads_start to merge."
                    )
    if not all(
        issubclass(sv.sv_type, SVpatterns.SVpatternSingleBreakend) for sv in breakends
    ):
        raise ValueError(
            "All SVcomposites must be of type SVpatternSingleBreakend to merge them."
        )

    if len(breakends) == 0:
        return []
    uf = UnionFind(range(len(breakends)))
    for i in range(len(breakends)):
        for j in range(i + 1, len(breakends)):
            if verbose:
                print(
                    f"Checking breakends {i} and {j}: {breakends[i].svPatterns[0].samplename}:{breakends[i].svPatterns[0].consensusID} vs {breakends[j].svPatterns[0].samplename}:{breakends[j].svPatterns[0].consensusID}"
                )
            if can_merge_svComposites_breakends(
                a=breakends[i],
                b=breakends[j],
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
                d=d,
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                verbose=verbose,
            ):
                uf.union(i, j)
            else:
                if verbose:
                    print(f"Not merging breakends {i} and {j}")
    result: list[SVcomposite] = []
    for cc in uf.get_connected_components(
        allow_singletons=True
    ):  # connected components of svComposites to merge
        result.append(
            SVcomposite.from_SVpatterns(
                svPatterns=[
                    svPattern for idx in cc for svPattern in breakends[idx].svPatterns
                ]
            )
        )
    return result


def _complexes_overlap(complex_a: SVcomposite, complex_b: SVcomposite, near: int) -> bool:
    """Check if two complex SVcomposites overlap in all their breakends within a given tolerance.
    It is assumed that each sv complex is composed of exactly one sv pattern at this point."""
    if len(complex_a.svPatterns) != 1 or len(complex_b.svPatterns) != 1:
        raise ValueError(
            f"_complexes_overlap: SVcomposites must contain exactly one SVpattern each to check for overlap. Got {len(complex_a.svPatterns)} and {len(complex_b.svPatterns)}."
        )
    
    primitives_a = sorted(complex_a.svPatterns[0].SVprimitives, key=lambda p: (p.chr, p.ref_start))
    primitives_b = sorted(complex_b.svPatterns[0].SVprimitives, key=lambda p: (p.chr, p.ref_start))

    if len(primitives_a) != len(primitives_b):
        return False  # Different number of breakends, cannot overlap

    for prim_a, prim_b in zip(primitives_a, primitives_b):
        if prim_a.chr != prim_b.chr:
            return False  # Different chromosomes, cannot overlap
        start_a, end_a = sorted((prim_a.ref_start, prim_a.ref_end))
        start_b, end_b = sorted((prim_b.ref_start, prim_b.ref_end))

        # Check if the intervals overlap within the tolerance
        if end_a + near < start_b or end_b + near < start_a:
            return False  # No overlap within tolerance

    return True  # All breakends overlap within tolerance


def _complexes_kmer_similarity(
    complex_a: SVcomposite,
    complex_b: SVcomposite,
    min_kmer_overlap: float,
) -> bool:
    """Calculate k-mer similarity between two complex SVcomposites based on their breakend sequences."""
    # each complex has exactly one sv pattern at this point
    if len(complex_a.svPatterns) != 1 or len(complex_b.svPatterns) != 1:
        raise ValueError(
            f"_complexes_kmer_similarity: SVcomposites must contain exactly one SVpattern each to calculate k-mer similarity. Got {len(complex_a.svPatterns)} and {len(complex_b.svPatterns)}."
        )
    if not issubclass(complex_a.sv_type, SVpatterns.SVpatternComplex):
        raise ValueError(
            f"_complexes_kmer_similarity: SVcomposite 'a' must be of type SVpatternComplex. Got {complex_a.sv_type}."
        )
    if not issubclass(complex_b.sv_type, SVpatterns.SVpatternComplex):
        raise ValueError(
            f"_complexes_kmer_similarity: SVcomposite 'b' must be of type SVpatternComplex. Got {complex_b.sv_type}."
        )
    # each sv pattern has saved sequence contexts (400 bp of consenus sequence, just like an inserrtion)
    # they are indexed by their sv primitive ID
    # this is not necessarily in the same order as the sv primitves sorted by reference chr, start_pos.
    # so it is necessary to keep track of the indexes.
    # sort sv primitives by chr, start_pos to get a consistent order. keep track of their old ID, so they can be mapped back to the sequences.
    
    # generate tuples of (original sv primtive index, sv primitive, sequence context)
    data_a = sorted(
        [
            (idx, sv_prim, complex_a.svPatterns[0].get_sequence_context(svp_index=idx))
            for idx, sv_prim in enumerate(complex_a.svPatterns[0].SVprimitives)
        ],        key=lambda x: (x[1].chr, x[1].ref_start),
    )
    data_b = sorted(
        [
            (idx, sv_prim, complex_b.svPatterns[0].get_sequence_context(svp_index=idx))
            for idx, sv_prim in enumerate(complex_b.svPatterns[0].SVprimitives)
        ],        key=lambda x: (x[1].chr, x[1].ref_start),
    )
    # now calculate kmer similarity for each group
    group_a = [seq for _, _, seq in data_a]
    group_b = [seq for _, _, seq in data_b]
    similarity: float = kmer_similarity_of_groups(
        group_a=group_a, group_b=group_b
    )
    return similarity >= min_kmer_overlap
    

def merge_complexes(
    complexes: list[SVcomposite],
    near: int,
    min_kmer_overlap: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    """Merge complex SVcomposites that are spatially close in every underlying sv primitive (breakend), and if their kmer similarities are high enough."""
    # any combination of complexes can be merged here, so every combination of complexes must be tested.
    uf: UnionFind = UnionFind(range(len(complexes)))
    # if two complexes match, connect them. Every connected component will be merged into one complex at the end.
    for i in range(len(complexes)):
        for j in range(i + 1, len(complexes)):
            if _complexes_overlap(complex_a=complexes[i], complex_b=complexes[j], near=near):
                if _complexes_kmer_similarity(
                    complex_a=complexes[i],
                    complex_b=complexes[j],
                    min_kmer_overlap=min_kmer_overlap,
                ):
                    uf.union(i, j)
                    if verbose:
                        print(
                            f"Merging complexes {i} and {j} based on overlap and k-mer similarity."
                        )
                else:
                    if verbose:
                        print(
                            f"Not merging complexes {i} and {j}: k-mer similarity too low."
                        )
            else:
                if verbose:
                    print(
                        f"Not merging complexes {i} and {j}: no overlap within tolerance."
                    )
    result: list[SVcomposite] = []
    for cc in uf.get_connected_components(allow_singletons=True):
        result.append(
            SVcomposite.from_SVpatterns(
                svPatterns=[
                    svPattern for idx in cc for svPattern in complexes[idx].svPatterns
                ]
            )
        )

    return result

# %%


# horizontal merge
def generate_svComposites_from_dbs(
    input: list[str | Path],
    sv_types: set[type[SVpatterns.SVpatternType]],
    candidate_regions_filter: dict[str, set[int]] | None = None,
) -> list[SVcomposite]:
    svComposites: list[SVcomposite] = []
    for p in input:
        # Get metadata to determine samplename
        metadata = svirltile.get_metadata(Path(p))
        samplename = metadata.get("samplename")

        # Determine if we should filter by crIDs for this sample
        crIDs_to_query: set[int] | None = None
        if (
            candidate_regions_filter is not None
            and samplename in candidate_regions_filter
        ):
            crIDs_to_query = candidate_regions_filter[samplename]
            log.info(
                f"Querying {len(crIDs_to_query)} candidate regions for sample {samplename}"
            )

        # Read SVpatterns from database, optionally filtered by crIDs
        svPatterns = SVpatterns.read_svPatterns_from_db(
            database=Path(p), crIDs=crIDs_to_query
        )
        # DEBUG START
        # if there are sv patterns with consensusID "7.x", then print their consensusID and sv_type
        # for svp in svPatterns:
        #     if svp.consensusID.startswith("7."):
        #         print(
        #             f"generate_svComposites_from_dbs:: {svp.consensusID} - {svp.get_sv_type()}"
        #         )
        # note: 7.0, 7.1, 7.2 are here
        # DEBUG END

        # count the type of each svPattern
        # then print the counts
        sv_type_counts: dict[str, int] = {}
        for svp in svPatterns:
            sv_type_counts[svp.get_sv_type()] = (
                sv_type_counts.get(svp.get_sv_type(), 0) + 1
            )
        log.info(
            f"Retrieved {len(svPatterns)} SVpatterns for sample {samplename}: {sv_type_counts}"
        )

        svComposites.extend(
            svPatterns_to_horizontally_merged_svComposites(
                svPatterns, sv_types=sv_types
            )
        )

        # after horizontal merging, count again teh types of sv composites
        sv_composite_type_counts: dict[str, int] = {}
        for svc in svComposites:
            sv_composite_type_counts[svc.sv_type.get_sv_type()] = (
                sv_composite_type_counts.get(svc.sv_type.get_sv_type(), 0) + 1
            )
        log.debug(
            f"After horizontal merging, total SVcomposites for sample {samplename}: {sv_composite_type_counts}"
        )

    return svComposites

def merge_svComposites_across_chromosomes(
    svComposites: list[SVcomposite],
    apriori_size_difference_fraction_tolerance: float,
    max_cohens_d: float,
    near: int,
    min_kmer_overlap: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    """
    Merge SVcomposites across chromosomes without parallelization.
    This function processes all chromosomes sequentially.
    """
    log.debug(f"Processing {len(svComposites)} SVcomposites across chromosomes")
    
    # build overlap trees for all chromosomes
    overlap_trees: dict[str, IntervalTree] = defaultdict(IntervalTree)
    for idx, svComposite in enumerate(svComposites):
        for svprimitive in svComposite.svPatterns[0].SVprimitives: # assuming there is only one sv pattern per sv composite at this point
            chr = svprimitive.chr
            start = svprimitive.ref_start
            end = svprimitive.ref_end
            if start > end:
                start, end = end, start
            if end - start < near:
                start -= near
                end += near
            if start < 0:
                start = 0
            overlap_trees[chr].addi(start, end, {idx})
    log.debug(f"Merging overlap trees for all chromosomes")
    uf = UnionFind(range(len(svComposites)))
    
    for _chr_name, tree in overlap_trees.items():
        tree.merge_overlaps(data_reducer=lambda x, y: x.union(y))

    for _chr_name, tree in overlap_trees.items():
        for interval in tree:
            indices = sorted(interval.data)
            if len(indices) > 1:
                for i in range(len(indices)):
                    for j in range(i + 1, len(indices)):
                        uf.union_by_name(indices[i], indices[j])
    
    # vertical merging
    connected_components = uf.get_connected_components(allow_singletons=True)
    merged_svComposites = []
    log.debug(f"Vertically merging {len(connected_components)} groups across chromosomes")
    for cc in tqdm(connected_components, desc=f"All Chromosomes"):
        if verbose:
            log.info(f"Merging group with {len(cc)} svComposites: {cc}")
        svComposite_group = [svComposites[idx] for idx in cc]

        if len(svComposite_group) > 1:
            log.debug(
                f"Vertically merging group of size {len(svComposite_group)} with consensusIDs: {[svc.svPatterns[0].consensusID for svc in svComposite_group]}"
            )

            merged = vertically_merged_svComposites_from_group(
                group=svComposite_group,
                d=max_cohens_d,
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
                verbose=verbose,
            )

            merged_svComposites.extend(merged)
        else:
            merged_svComposites.append(svComposite_group[0])
    log.debug(
        f"Merged {len(svComposites)} SVcomposites into {len(merged_svComposites)} across chromosomes"
    )
    return merged_svComposites
    
        

def merge_svComposites_for_chromosome(
    chr_name: str,
    svComposites: list[SVcomposite],
    apriori_size_difference_fraction_tolerance: float,
    max_cohens_d: float,
    near: int,
    min_kmer_overlap: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    """
    Merge SVcomposites for a single chromosome.
    This function is designed to be run in parallel for different chromosomes.
    """
    log.debug(f"Processing chromosome {chr_name} with {len(svComposites)} SVcomposites")

    # Build overlap trees for this chromosome only
    overlaptree = IntervalTree()
    for idx, svComposite in enumerate(svComposites):
        for svp in svComposite.svPatterns:
            genome_intervals = svp.get_consensus_aln_intervals_from_svPrimitives()
            for interval_chr, start, end in genome_intervals:
                # Skip intervals not on this chromosome
                if interval_chr != chr_name:
                    continue

                if start > end:
                    start, end = end, start
                if end - start < near:
                    start -= near
                    end += near
                if start < 0:
                    start = 0
                overlaptree.addi(start, end, {idx})

    # Merge overlapping intervals
    log.debug(f"Merging overlap tree for chromosome {chr_name}")
    uf = UnionFind(range(len(svComposites)))
    overlaptree.merge_overlaps(data_reducer=lambda x, y: x.union(y))

    for interval in overlaptree:
        indices = list(set(interval.data))
        if len(indices) > 1:
            for i in range(len(indices)):
                for j in range(i + 1, len(indices)):
                    # Don't merge if same sample and consensus are identical - this is vertical merging (across samples and consensuses)
                    if (
                        svComposites[indices[i]].svPatterns[0].consensusID
                        == svComposites[indices[j]].svPatterns[0].consensusID
                        and svComposites[indices[i]].svPatterns[0].samplename
                        == svComposites[indices[j]].svPatterns[0].samplename
                    ):
                        continue
                    else:
                        uf.union_by_name(indices[i], indices[j])

    # Get connected components and merge groups
    connected_components = uf.get_connected_components(allow_singletons=True)
    merged_svComposites = []

    log.debug(
        f"Vertically merging {len(connected_components)} groups for chromosome {chr_name}"
    )
    for cc in tqdm(connected_components, desc=f"Chr {chr_name}"):
        if verbose:
            log.info(f"Merging group with {len(cc)} svComposites: {cc}")
        svComposite_group = [svComposites[idx] for idx in cc]

        if len(svComposite_group) > 1:
            log.debug(
                f"Vertically merging group of size {len(svComposite_group)} with consensusIDs: {[svc.svPatterns[0].consensusID for svc in svComposite_group]}"
            )
            merged = vertically_merged_svComposites_from_group(
                group=svComposite_group,
                d=max_cohens_d,
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
                verbose=verbose,
            )
            merged_svComposites.extend(merged)
        else:
            merged_svComposites.append(svComposite_group[0])

    log.debug(
        f"Chromosome {chr_name}: Merged {len(svComposites)} SVcomposites into {len(merged_svComposites)}"
    )
    return merged_svComposites

def merge_svComposites(
    svComposites: list[SVcomposite],
    apriori_size_difference_fraction_tolerance: float,
    max_cohens_d: float,
    near: int,
    min_kmer_overlap: float,
    verbose: bool = False,
    threads: int = 1,
) -> list[SVcomposite]:
    """
    Merge SVcomposites across samples, optionally using parallel processing per chromosome.
    """
    # Group SVcomposites by chromosome
    log.debug(f"Grouping {len(svComposites)} SVcomposites by chromosome...")
    chr_groups: dict[str, list[SVcomposite]] = defaultdict(list)

    # not part of this merging should be all types of SV composites that are of a translocation type
    # since the work is parallelized across chromosomes.

    cpx_groups: list[SVcomposite] = []
    # how to handle the cpx groups?
    # since they can have parts hopping to homologous regions, the chromosome pattern
    # can easily be disturbed. It could work much better to align all associated consensus
    # sequences (the core sequences) in a all-vs-all alignment and then cluster them based partial overlaps.
    # The alignment needs to be very strict and allow only small indels and mismatches.
    # the condition of a merge attempt would be that the reciprocal overlap is greater than 50%
    # however, it must be prevented that contradictory merges are made, e.g. A-B and B-C but A-C is not valid.

    for svComposite in svComposites:
        if svComposite.sv_type in (
            SVpatterns.SVpatternComplex,
            SVpatterns.SVpatternTranslocation,
            SVpatterns.SVpatternInvertedTranslocation,
            SVpatterns.SVpatternInversionTranslocation,
        ):
            cpx_groups.append(svComposite)
        else:
            chr_name = svComposite.ref_start[0]  # Get chromosome from ref_start tuple
            chr_groups[chr_name].append(svComposite)

    log.info(f"Found {len(chr_groups)} chromosomes with variants")
    for chr_name, group in chr_groups.items():
        log.info(f"  {chr_name}: {len(group)} SVcomposites")

    merged_svComposites = []

    if threads > 1:
        # Parallel processing per chromosome
        log.info(f"Processing chromosomes in parallel with {threads} workers")

        from functools import partial

        worker_fn = partial(
            merge_svComposites_for_chromosome,
            apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
            max_cohens_d=max_cohens_d,
            near=near,
            min_kmer_overlap=min_kmer_overlap,
            verbose=verbose,
        )

        with ProcessPoolExecutor(max_workers=threads) as executor:
            # Submit jobs for each chromosome
            futures = {
                executor.submit(worker_fn, chr_name, chr_svComposites): chr_name
                for chr_name, chr_svComposites in chr_groups.items()
            }

            # Collect results as they complete
            from concurrent.futures import as_completed

            for future in tqdm(
                as_completed(futures), total=len(futures), desc="Processing chromosomes"
            ):
                chr_name = futures[future]
                try:
                    chr_merged = future.result()
                    merged_svComposites.extend(chr_merged)
                except Exception as e:
                    log.error(f"Error processing chromosome {chr_name}: {e}")
                    raise
    else:
        # Sequential processing
        log.info("Processing chromosomes sequentially")
        for chr_name, chr_svComposites in tqdm(
            chr_groups.items(), desc="Processing chromosomes"
        ):
            chr_merged = merge_svComposites_for_chromosome(
                chr_name=chr_name,
                svComposites=chr_svComposites,
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
                max_cohens_d=max_cohens_d,
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                verbose=verbose,
            )
            merged_svComposites.extend(chr_merged)

    # cpx_groups can jump across chromosomes, so they cannot be parallelized by chromosome.
    # they are merged here sequentially.
    merged_svComposites.extend(merge_svComposites_across_chromosomes(
        svComposites=cpx_groups,
        apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
        max_cohens_d=max_cohens_d,
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        verbose=verbose,
    ))

    log.info(
        f"Merged {len(svComposites)} SVcomposites into {len(merged_svComposites)} across all chromosomes."
    )
    return merged_svComposites


# %%


@attrs.define
class Genotype:
    samplename: str
    genotype: str  # e.g. "0/0", "0/1", "1/1"
    gt_likelihood: float
    genotype_quality: int  # GQ phred of gt_likelihood
    total_coverage: int  # TC
    ref_reads: int  # DR
    var_reads: int  # DV

    def to_format_field(self, FORMAT_field="GT:GQ:TC:DR:DV:GP") -> str:
        properties = {
            "GT": self.genotype,
            "GQ": self.genotype_quality,
            "TC": self.total_coverage,
            "DR": self.ref_reads,
            "DV": self.var_reads,
            "GP": self.gt_likelihood,
        }
        return ":".join([str(properties[k]) for k in FORMAT_field.split(":")])


def create_wild_type_genotype(samplename: str, total_coverage: int) -> Genotype:
    return Genotype(
        samplename=samplename,
        genotype="0/0",
        gt_likelihood=1.0,
        genotype_quality=60,
        total_coverage=total_coverage,
        ref_reads=total_coverage,
        var_reads=0,
    )


@attrs.define
class SVcall:
    genotypes: dict[str, Genotype]  # samplename -> Genotype
    passing: bool
    chrname: str
    end: int
    start: int
    svtype: str  # e.g. "INS", "DEL", "BND", "INV", "DUP"
    svlen: int  # length of the SV, e.g. for INS it is the length of the inserted sequence, for DEL it is the length of the deleted sequence
    pass_altreads: int
    pass_gq: int
    precise: bool
    mateid: str  # list of mate breakends' IDs
    consensusIDs: list[
        str
    ]  # IDs of the consensus that this SV originates from. Other consensus sequences can also be involved.
    ref_sequence: bytes | None = None
    alt_sequence: bytes | None = None
    sequence_id: str | None = None  # ID for sequence in FASTA file if symbolic
    description: dict[str, str] | None = (
        None  # optional description field for additional annotations; e.g. outer and inner intervals of inversions or duplications
    )

    def to_vcf_line(
        self,
        vcfIDnumber: int,
        samplenames: list[str],
        covtrees: dict[str, dict[str, IntervalTree]],
        refdict: dict[str, str],
        symbolic_threshold: int,
        ONE_BASED: int = 1,
    ) -> str | None:
        info_fields: dict = {
            "PASS_ALTREADS": self.pass_altreads,
            "pass_GQ": self.pass_gq,
            "SVTYPE": self.svtype,
            "END": str(
                self.end + ONE_BASED
            ),  # end is inclusive in vcf, so we need to subtract 1 -> 'END':str(self.end+ONE_BASED-1)
            "SVLEN": str(self.svlen),
            "CONSENSUSIDs": ",".join(self.consensusIDs),
            "MATEID": self.mateid if len(self.mateid) > 0 else "NA",
        }

        # Add sequence ID if this is a symbolic allele
        if self.sequence_id:
            info_fields["SEQ_ID"] = self.sequence_id

        info_line = ";".join([f"{key}={value}" for key, value in info_fields.items()])
        info_line += ";" + ("PRECISE" if self.precise else "IMPRECISE")
        FORMAT_field = "GT:GQ:TC:DR:DV:GP"
        format_content: list[Genotype] = []
        for samplename in samplenames:
            gt: Genotype | None = self.genotypes.get(samplename, None)
            if gt is None:
                tree: IntervalTree | None = covtrees[samplename].get(self.chrname, None)
                found_intervals: set[Interval] = tree[self.start] if tree else set()
                coverage: int = len(found_intervals)
                format_content.append(
                    create_wild_type_genotype(
                        samplename=samplename, total_coverage=coverage
                    )
                )
            else:
                format_content.append(self.genotypes[samplename])
        # now construct the vcf line
        vcfID = f"{self.svtype}.{str(vcfIDnumber)}"
        format_line = [gt.to_format_field(FORMAT_field) for gt in format_content]
        refbase: str | None = refdict.get(f"{self.chrname}:{self.start}", None)
        if refbase is None:
            log.warning(
                f"Reference base not found for {self.chrname}:{self.start}. Please check the reference dictionary."
            )
            return None

        # Determine if we should use symbolic alleles
        ref_seq = pickle.loads(self.ref_sequence) if self.ref_sequence else ""
        alt_seq = pickle.loads(self.alt_sequence) if self.alt_sequence else ""

        use_symbolic = (
            len(ref_seq) > symbolic_threshold or len(alt_seq) > symbolic_threshold
        )

        if use_symbolic:
            # Use symbolic allele notation
            ref_str = refbase  # Just the anchor base
            alt_str = f"<{self.svtype}>"
        else:
            # Use explicit sequences for small variants
            ref_str = refbase + ref_seq
            alt_str = refbase + alt_seq

        return "\t".join([
            str(self.chrname),
            str(self.start + ONE_BASED),
            str(vcfID),  # ID
            ref_str,  # REF
            alt_str,  # ALT
            str(60),  # QUAL
            "PASS" if self.passing else "LowQual",  # FILTER
            info_line,
            FORMAT_field,
            *format_line,
        ])


def get_ref_reads_from_covtrees(
    samplename: str,
    chrname: str,
    start: int,
    end: int,
    covtrees: dict[str, dict[str, IntervalTree]],
    min_radius:int=1,
) -> set[int]:
    if end < start:
        raise ValueError(
            f"get_ref_reads_from_covtrees: end {end} is less than start {start} for sample {samplename} at {chrname}:{start}-{end}."
        )
    min_radius = abs(min_radius)
    if end-start < min_radius * 2:
        difference = min_radius * 2 - (end - start)
        start -= floor(difference * 0.5)
        end += ceil(difference * 0.5)
        if start < 0:
            start = 0
    all_reads: set[int] = {
        int(it.data)
        for it in covtrees[samplename][chrname][start:end]
    }
    return all_reads


def genotype_of_sample(
    samplename: str,
    chrname: str,
    start: int,
    end: int,
    raw_alt_reads: set[int],
    covtrees: dict[str, dict[str, IntervalTree]],
    cn_tracks: dict[str, dict[str, IntervalTree]],
    min_radius:int = 1,
) -> Genotype:
    genotype:Genotype
    if chrname not in covtrees.get(samplename, {}):
        log.warning(
            f"SVcalls_from_SVcomposite: Chromosome {chrname} not found in coverage tree for sample {samplename}. Assigning 0/0 genotype."
        )
        genotype = create_wild_type_genotype(
            samplename=samplename, total_coverage=0
        )
        return genotype
    all_reads: set[int] = get_ref_reads_from_covtrees(
        samplename=samplename,
        chrname=chrname,
        start=start,
        end=end,
        covtrees=covtrees,
        min_radius=min_radius,
    )
    alt_reads: set[int] = all_reads.intersection(raw_alt_reads)

    ref_reads: set[int] = all_reads.difference(alt_reads)

    if len(alt_reads) == 0:
        log.warning(
            f"SVcalls_from_SVcomposite: No alt reads found for sample {samplename} at {chrname}:{start}-{end}. Assigning 0/0 genotype."
        )
        genotype = create_wild_type_genotype(
            samplename=samplename, total_coverage=len(all_reads)
        )
        return genotype

    if len(all_reads) == 0:
        log.warning(
            f"SVcalls_from_SVcomposite: No ref/alt reads found for sample {samplename} at {chrname}:{start}-{end}. Assigning 0/0 genotype."
        )
        genotype = create_wild_type_genotype(
            samplename=samplename, total_coverage=0
        )
        return genotype

    # Query copy number for this locus from CN tracks
    copy_number : int = 2  # Default diploid - if no copy number track is available
    if samplename in cn_tracks and chrname in cn_tracks[samplename]:
        try:
            # Query overlapping intervals from the CN track IntervalTree
            overlapping_cn = cn_tracks[samplename][chrname][start:end]
            if overlapping_cn:
                # Take the most common CN in the overlapping intervals
                cn_values = [interval.data for interval in overlapping_cn]
                copy_number = (
                    max(set(cn_values), key=cn_values.count) if cn_values else 2
                )
                log.debug(
                    f"Copy number at {chrname}:{start}-{end} for {samplename}: CN={copy_number} (from {len(cn_values)} overlapping bins)"
                )
            else:
                log.debug(
                    f"No CN data overlapping {chrname}:{start}-{end} for {samplename}, using default CN=2"
                )
        except Exception as e:
            log.warning(
                f"Error querying copy number for {samplename} at {chrname}:{start}-{end}: {e}. Using default CN=2"
            )
            # default copy number remains
    else:
        log.debug(
            f"No CN tracks available for sample {samplename} at {chrname}, using default CN=2"
        )

    # first compute the total gt for this locus and this sample (for all SVcomposites, that once overlapped)
    gt_likelihoods: dict[str, float] = genotype_likelihood(
        n_alt_reads=len(alt_reads), n_total_reads=len(all_reads), cn=copy_number
    )

    gt = max(gt_likelihoods.items(), key=lambda x: x[1])[0]
    # if not gt in ('0/0', '0/1', '1/1'):
    #     raise ValueError(f"Invalid genotype {gt} for sample {samplename}, in region: {chrname,start,end}, with {len(alt_reads)} alt reads and {len(all_reads)} total reads.")

    genotype = Genotype(
        samplename=samplename,
        genotype=gt,
        gt_likelihood=gt_likelihoods[gt],
        genotype_quality=(
            int(-10 * np.log10(1.0 - gt_likelihoods[gt]))
            if gt_likelihoods[gt] < 1.0
            else 60
        ),
        total_coverage=len(all_reads),
        ref_reads=len(ref_reads),
        var_reads=len(alt_reads),
    )
    return genotype


def SVcalls_from_SVcomposite(
    svComposite: SVcomposite,
    covtrees: dict[str, dict[str, IntervalTree]],
    cn_tracks: dict[str, dict[str, IntervalTree]],  # saplename -> chrname -> IntervalTree with copy number
    find_leftmost_reference_position: bool,
) -> list[SVcall]:
    # intra-alignment fragment variants (closed locus), e.g. INS, DEL, INV, DUP
    #   have one chr, start, end on the reference
    # inter-alignment fragment variants (interrupted locus), e.g translocations, break points
    #   have multiple chr, start, end on the reference

    # the proceedure is:
    # 1) get the reference interval(s) for the svComposite
    # 2) for each sample:
    #    a) get all reads covering the reference interval(s) from the coverage tree
    #    b) get alt reads from the svComposite
    #    c) compute ref reads = all reads - alt reads
    #    d) compute genotype likelihoods based on alt reads, total reads, and copy number
    #    e) assign genotype with highest likelihood
    # 3) create SVcall object with genotypes for all samples
    
    # the two cases (closed locus vs interrupted locus) are handled a little differently,
    # as two or more SVcalls are generated for the interrupted locus case. Two for each connected break point.
    # They are then linked via the mate ID field.
    # for closed locus variants, only one SVcall is generated.
    
    # necessary functions:
    # - get reads covering locus from covtree (region(chr, start, end), samplename, covtree) -> tuple[set(alt readIDs), set(ref readIDs)]
    # - get gt and gt likelihood from alt reads, total reads, copy number -> Genotype
    # - get copy number of sample in locus (region(chr, start, end), samplename, cn_tracks) -> int
    # - genotype likelihood calculation (n_alt_reads, n_total_reads, copy_number) -> dict[genotype: likelihood]
    
    
    
    # separate cases:
    # 1) insertion & deletion
    if not svComposite.sv_type in SUPPORTED_SV_TYPES:
        log.warning(
            f"SVcalls_from_SVcomposite: svComposite of type {svComposite.sv_type.__name__} is not supported for SVcall generation. Supported types are: {', '.join(SUPPORTED_SV_TYPE_STRINGS)}"
        )
        return []
    all_alt_reads: dict[str, set[int]] = (
        svComposite.get_alt_readnamehashes_per_sample()
    )  # {samplename: {readname, ...}}

    if issubclass(svComposite.sv_type, SVpatterns.SVpatternDeletion) \
            or issubclass(svComposite.sv_type, SVpatterns.SVpatternInsertion) \
            or issubclass(svComposite.sv_type, SVpatterns.SVpatternInversion) \
            or issubclass(svComposite.sv_type, SVpatterns.SVpatternSingleBreakend):
        return [svcall_object_from_svcomposite(
            svComposite=svComposite,
            covtrees=covtrees,
            cn_tracks=cn_tracks,
            find_leftmost_reference_position=find_leftmost_reference_position,
            all_alt_reads=all_alt_reads)]
    elif issubclass(svComposite.sv_type, SVpatterns.SVpatternComplex):
        # other types that are represented as connected break ends will follow!
        return svcall_objects_from_svcomposite(
            svComposite=svComposite,
            covtrees=covtrees,
            cn_tracks=cn_tracks,
            all_alt_reads=all_alt_reads)
    log.warning(
        f"SVcalls_from_SVcomposite: svComposite of type {svComposite.sv_type.__name__} in regions {svComposite.get_regions()} is not yet implemented for SVcall generation."
    )
    return []

def svcall_object_from_svcomposite(
        svComposite:SVcomposite,
        covtrees:dict[str, dict[str, IntervalTree]],
        cn_tracks:dict[str, dict[str, IntervalTree]],
        find_leftmost_reference_position:bool,
        all_alt_reads:dict[str, set[int]]) -> SVcall:
    chrname, start, end = get_svComposite_interval_on_reference(
            svComposite=svComposite,
            find_leftmost_reference_position=find_leftmost_reference_position,
        )
    svlen: int = abs(svComposite.get_size())
        # Get read IDs supporting the SV. arguments: samplename, chrname, start, end, covtree
    consensusIDs: list[str] = list({
            svPattern.samplenamed_consensusID for svPattern in svComposite.svPatterns
        })
        
        # 
    genotypes: dict[str, Genotype] = {
            samplename:genotype_of_sample(
                samplename=samplename,
                chrname=chrname,
                start=start,
                end=end,
                raw_alt_reads=all_alt_reads[samplename],
                covtrees=covtrees,
                cn_tracks=cn_tracks,
            ) for samplename in all_alt_reads.keys()
        }
        # start and end need to be adjusted based on sv_type and need to be determined for the whole sv composite   
        
    if issubclass(svComposite.sv_type, SVpatterns.SVpatternDeletion):
        end = start + svlen
    elif issubclass(svComposite.sv_type, SVpatterns.SVpatternInsertion):
        end = start + 1

    alt_seq: bytes|None = pickle.dumps(svComposite.get_alt_sequence()) if svComposite.get_alt_sequence() else None
        
    ref_seq: bytes|None = pickle.dumps(svComposite.get_ref_sequence()) if svComposite.get_ref_sequence() else None

        # quality filters
    pass_altreads: bool = (
            max(genotypes.items(), key=lambda x: x[1].var_reads)[1].var_reads >= 3
        )
    pass_gq = True
    passing: bool = pass_altreads and pass_gq

        # call is not precise if it is in a repeat. Check for repeatIDs
    precise: bool = not bool(
            sum(
                len(svPrimitive.repeatIDs)
                for svPattern in svComposite.svPatterns
                for svPrimitive in svPattern.SVprimitives
            )
        )

    return SVcall(
                genotypes=genotypes,
                passing=passing,
                chrname=chrname,
                end=end,
                start=start,
                svtype=svComposite.sv_type.get_sv_type(),
                svlen=svlen,
                pass_altreads=pass_altreads,
                pass_gq=pass_gq,
                precise=precise,
                mateid="",
                consensusIDs=consensusIDs,
                ref_sequence=ref_seq,
                alt_sequence=alt_seq,
            )

def svcall_objects_from_svcomposite(
        svComposite:SVcomposite,
        covtrees:dict[str, dict[str, IntervalTree]],
        cn_tracks:dict[str, dict[str, IntervalTree]],
        all_alt_reads:dict[str, set[int]]) -> list[SVcall]:
    """Generate SVcall objects from a complex SVcomposite with multiple breakends.
    The catch is that each breakend will generate its own SVcall object, and they will be linked via the mateid field.
    """
    consensusIDs : list[str] = list({svPattern.samplenamed_consensusID for svPattern in svComposite.svPatterns})
    results : list[SVcall] = []
    # genotypes are a bit more tricky now. Since there are multiple breakends, we need to get genotypes for each breakend.
    for (samplename, consensusID), group in svComposite.get_all_groups().items():
        if len(group) > 1 or len(group[0].SVprimitives) != 2 or not all(svp.sv_type in (3,4) for svp in group[0].SVprimitives):
            raise NotImplementedError(
                f"svcall_objects_from_svcomposite: Currently only single breakends per sample-consensusID are supported. Found {len(group)} breakends for sample {samplename}, consensusID {consensusID} in SVcomposite: {svComposite}"
            )
        # each svPattern represents a pair of BNDs.
        


# SVcall to breakends parsing
# SVcall to ins/del parsing

def get_svComposite_interval_on_reference(
    svComposite: SVcomposite, find_leftmost_reference_position: bool
) -> tuple[str, int, int]:
    if svComposite.sv_type not in SUPPORTED_SV_TYPES:
        raise ValueError(
            f"get_svComposite_indel_interval called with svComposite that is neither {', '.join(t.__name__ for t in SUPPORTED_SV_TYPES)}: {svComposite}"
        )
    weighted_regions: list[tuple[tuple[str, int, int], int]] = []
    for svPattern in svComposite.svPatterns:
        if svPattern.get_sv_type() not in SUPPORTED_SV_TYPE_STRINGS:
            raise ValueError(
                f"All svPatterns of a SVcomposite need to be of supported types {', '.join(SUPPORTED_SV_TYPE_STRINGS)}. The SVcomposite is: {svComposite}"
            )
        region: tuple[str, int, int] = svPattern.get_reference_region()
        if find_leftmost_reference_position:
            weight = min(
                svprimitive.ref_start for svprimitive in svPattern.SVprimitives
            )
        else:
            weight = len(svPattern.get_supporting_reads() * svPattern.get_size())
        weighted_regions.append((region, weight))
    # pick the winning region
    if find_leftmost_reference_position:
        return min(
            weighted_regions, key=lambda x: x[1]
        )[
            0
        ]  # reports the leftmost SVpattern, instead of the one with most supporting reads * size. This might be better aligned with the giab SV benchmark, but should be discussed in the paper.
    else:
        return max(weighted_regions, key=lambda x: x[1])[0]


# %% VCF file stuff


def generate_header(
    reference: Path, samplenames: list[str], fasta_path: Path | None = None
) -> list[str]:
    reference = Path(reference)
    header = [
        "##fileformat=VCFv4.2",
        f"##fileDate={datetime.now().strftime('%Y%m%d')}",
        f"##reference=file://{str(reference.absolute())}",
    ]

    # Add FASTA reference if provided
    if fasta_path:
        header.append(f"##sequences=file://{str(fasta_path.absolute())}")

    # add contigs to header
    # load reference index as pd dataframe
    ref_index = read_csv(
        reference.with_suffix(reference.suffix + ".fai"), sep="\t", header=None
    )
    # compute lengths of contigs
    # add contigs to header
    for i in range(ref_index.shape[0]):
        header.append(
            f"##contig=<ID={str(ref_index.iloc[i, 0])},length={ref_index.iloc[i, 1]}>"
        )
    header.append(
        '##FILTER=<ID=LowQual,Description="Poor quality and insufficient number of informative reads.">'
    )
    header.append(
        '##FILTER=<ID=PASS,Description="high quality and sufficient number of informative reads.">'
    )
    header.append(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">'
    )
    header.append(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, INS=Insertion, DUP=Duplication, INV=Inversion">'
    )
    header.append(
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">'
    )
    header.append(
        '##INFO=<ID=PASS_ALTREADS,Number=1,Type=String,Description="Passed alt reads threshold">'
    )
    header.append(
        '##INFO=<ID=pass_GQ,Number=1,Type=String,Description="Passed Genotype precision threshold">'
    )
    header.append(
        '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">'
    )
    header.append(
        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">'
    )
    header.append(
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">'
    )
    header.append(
        '##INFO=<ID=CONSENSUSIDs,Number=1,Type=String,Description="ID of the consensus that this SV originiates from. Other consensus sequences can also be involved.">'
    )
    header.append(
        '##INFO=<ID=SEQ_ID,Number=1,Type=String,Description="ID of sequence in companion FASTA file for symbolic alleles">'
    )
    header.append('##ALT=<ID=INS,Description="Insertion">')
    header.append('##ALT=<ID=DEL,Description="Deletion">')
    header.append('##ALT=<ID=DUP,Description="Duplication">')
    header.append('##ALT=<ID=INV,Description="Inversion">')
    header.append('##ALT=<ID=BND,Description="Breakend; Translocation">')
    header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.append(
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">'
    )
    header.append('##FORMAT=<ID=TC,Number=1,Type=Integer,Description="Total coverage">')
    header.append(
        '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">'
    )
    header.append(
        '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">'
    )
    header.append(
        '##FORMAT=<ID=GP,Number=1,Type=Float,Description="Genotype probability">'
    )

    header.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samplenames)
    )
    return header


def genotype_likelihood(
    n_alt_reads: int, n_total_reads: int, cn: int
) -> dict[str, float]:
    """
    Compute genotype likelihoods given observed alt/total reads and copy number.

    Args:
        n_alt_reads: Number of reads supporting the alternate allele
        n_total_reads: Total number of reads covering the locus
        cn: Copy number at this locus (1, 2, 3, 4, etc.)

    Returns:
        Dictionary mapping genotype strings to their likelihoods
    """
    if n_total_reads == 0:
        return {"0/0": 1.0, "0/1": 0.0, "1/1": 0.0}

    # Generate all possible genotypes for this copy number
    # For CN=n, we can have 0 to n copies of the alt allele
    genotype_probs = {}

    if cn == 1:
        # Haploid (e.g., male X/Y chromosomes, or deletions reducing CN to 1)
        # Possible genotypes: 0 (ref) or 1 (alt)
        genotype_probs = {
            "0": 0.01,  # ~0% alt reads expected
            "1": 0.50,  # ~50% alt reads expected (bias towards variant calling)
        }
    elif cn == 2:
        # Diploid (normal autosomal)
        # Possible genotypes: 0/0, 0/1, 1/1
        genotype_probs = {
            "0/0": 0.01,  # ~0% alt reads expected
            "0/1": 0.45,  # ~45% alt reads expected
            "1/1": (
                0.90  # ~90% alt reads expected (compensate n_total_reads overestimation)
            ),
        }
    elif cn == 3:
        # Triploid (duplication)
        # Possible genotypes: 0/0/0, 0/0/1, 0/1/1, 1/1/1
        genotype_probs = {
            "0/0/0": 0.01,  # ~0% alt reads
            "0/0/1": 1.0 / 3.0,  # ~33% alt reads
            "0/1/1": 2.0 / 3.0,  # ~67% alt reads
            "1/1/1": 0.99,  # ~100% alt reads
        }
    elif cn == 4:
        # Tetraploid
        # Possible genotypes: 0/0/0/0, 0/0/0/1, 0/0/1/1, 0/1/1/1, 1/1/1/1
        genotype_probs = {
            "0/0/0/0": 0.01,  # ~0% alt reads
            "0/0/0/1": 0.25,  # ~25% alt reads
            "0/0/1/1": 0.50,  # ~50% alt reads
            "0/1/1/1": 0.75,  # ~75% alt reads
            "1/1/1/1": 0.99,  # ~100% alt reads
        }
    else:
        # General case for CN > 4 or CN == 0
        if cn == 0:
            # Homozygous deletion - no copies, should have no reads
            return {"0": 1.0}

        # For higher CN, generate genotypes dynamically
        for n_alt_copies in range(cn + 1):
            genotype = "/".join(["1" if i < n_alt_copies else "0" for i in range(cn)])
            expected_alt_fraction = n_alt_copies / cn if cn > 0 else 0.0
            # Use small epsilon for 0 and 1 to avoid edge effects
            if expected_alt_fraction == 0.0:
                expected_alt_fraction = 0.01
            elif expected_alt_fraction == 1.0:
                expected_alt_fraction = 0.99
            genotype_probs[genotype] = expected_alt_fraction

    # Compute binomial probabilities for each genotype
    likelihoods = {
        genotype: binom.pmf(n_alt_reads, n_total_reads, p)
        for genotype, p in genotype_probs.items()
    }

    # Normalize to get posterior probabilities
    total_likelihood = sum(likelihoods.values())

    if total_likelihood <= 0.0:
        log.warning(
            f"Total likelihood is zero or negative for n_alt_reads={n_alt_reads}, "
            f"n_total_reads={n_total_reads}, CN={cn}. Returning uniform probabilities."
        )
        # Return uniform distribution over all genotypes
        uniform_prob = 1.0 / len(genotype_probs)
        return dict.fromkeys(genotype_probs.keys(), uniform_prob)

    probabilities = {
        genotype: float(likelihood / total_likelihood)
        for genotype, likelihood in likelihoods.items()
    }

    # Replace any NaN values with 0.0
    for genotype, probability in probabilities.items():
        if np.isnan(probability):
            probabilities[genotype] = 0.0

    return probabilities


# %%


def reference_bases_by_merged_svComposites(
    svComposites: list[SVcomposite],
    reference: Path,
    find_leftmost_reference_position: bool,
) -> dict[str, str]:
    # Collect all positions that need reference bases

    positions = []
    for svComposite in svComposites:
        chrname, start, end = get_svComposite_interval_on_reference(
            svComposite=svComposite,
            find_leftmost_reference_position=find_leftmost_reference_position,
        )
        positions.append((chrname, start))

    unique_positions = list(dict.fromkeys(positions))
    if not unique_positions:
        return {}
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False
    ) as tmp_regions:
        for chrname, start in unique_positions:
            print(f"{chrname}:{start}-{start}", file=tmp_regions)
        tmp_regions_path = tmp_regions.name
    try:
        cmd_faidx = f"samtools faidx --region-file {tmp_regions_path} {str(reference)}"
        result = subprocess.run(
            split(cmd_faidx), capture_output=True, text=True, check=True
        )

        dict_reference_bases = {}
        lines = result.stdout.strip().split("\n")
        i = 0
        while i < len(lines):
            if lines[i].startswith(">"):
                # Parse header line like ">chr1:123-124"
                header = lines[i][1:]  # Remove '>'
                if ":" in header and "-" in header:
                    chrom_pos = header.split("-")[0]  # Get "chr1:123" part
                    if i + 1 < len(lines) and not lines[i + 1].startswith(">"):
                        base = lines[i + 1].upper()
                        dict_reference_bases[chrom_pos] = base
                        i += 2
                    else:
                        i += 1
                else:
                    i += 1
            else:
                i += 1

    finally:
        # Clean up temporary file
        Path(tmp_regions_path).unlink()

    return dict_reference_bases


def write_sequences_to_fasta(
    svCalls: list[SVcall], output_path: Path, symbolic_threshold: int
) -> None:
    """Write large SV sequences to a FASTA file."""
    with open(output_path, "w") as f:
        for svcall in svCalls:
            if not svcall.sequence_id:
                continue

            # Write REF sequence if large
            if svcall.ref_sequence:
                ref_seq = pickle.loads(svcall.ref_sequence)
                if len(ref_seq) > symbolic_threshold:
                    print(f">{svcall.sequence_id}_REF", file=f)
                    # Write sequence in 80-character lines
                    for i in range(0, len(ref_seq), 80):
                        print(ref_seq[i : i + 80], file=f)

            # Write ALT sequence if large
            if svcall.alt_sequence:
                alt_seq = pickle.loads(svcall.alt_sequence)
                if len(alt_seq) > symbolic_threshold:
                    print(f">{svcall.sequence_id}_ALT", file=f)
                    for i in range(0, len(alt_seq), 80):
                        print(alt_seq[i : i + 80], file=f)


def write_svCalls_to_vcf(
    svCalls: list[SVcall],
    samplenames: list[str],
    reference: Path | str,
    covtrees: dict[str, dict[str, IntervalTree]],
    refdict: dict[str, str],
    output: Path,
    symbolic_threshold: int,
) -> None:
    # Determine FASTA output path
    if str(output).endswith(".vcf.gz"):
        fasta_path = Path(str(output).replace(".vcf.gz", ".variants.fasta"))
    elif str(output).endswith(".vcf"):
        fasta_path = Path(str(output).replace(".vcf", ".variants.fasta"))
    else:
        fasta_path = output.with_suffix(".variants.fasta")

    # Assign sequence IDs to SVcalls that need them
    for idx, svCall in enumerate(svCalls):
        ref_seq = pickle.loads(svCall.ref_sequence) if svCall.ref_sequence else ""
        alt_seq = pickle.loads(svCall.alt_sequence) if svCall.alt_sequence else ""

        if len(ref_seq) > symbolic_threshold or len(alt_seq) > symbolic_threshold:
            svCall.sequence_id = f"{svCall.svtype}.{idx}"

    # Write sequences to FASTA
    write_sequences_to_fasta(svCalls, fasta_path, symbolic_threshold=symbolic_threshold)

    tmp_result_unsorted = tempfile.NamedTemporaryFile(delete=True, suffix=".vcf")
    with open(tmp_result_unsorted.name, "w") as f:
        # write header
        header = generate_header(
            reference=Path(reference), samplenames=samplenames, fasta_path=fasta_path
        )
        for line in header:
            print(line, file=f)
        # write SVcalls
        for vcfIDnum, svCall in enumerate(svCalls):
            line = svCall.to_vcf_line(
                ONE_BASED=1,
                samplenames=samplenames,
                covtrees=covtrees,
                vcfIDnumber=vcfIDnum,
                refdict=refdict,
                symbolic_threshold=symbolic_threshold,
            )
            if line is not None:
                print(line, file=f)
    if str(output).endswith(".vcf"):
        # sort and copy to result path
        cmd_sort = f"bcftools sort {str(tmp_result_unsorted.name)} -o {str(output)}"
        subprocess.check_call(split(cmd_sort))
    elif str(output).endswith(".vcf.gz"):
        tmp_sorted = tempfile.NamedTemporaryFile(delete=True, suffix=".vcf")
        cmd_sort = (
            f"bcftools sort {str(tmp_result_unsorted.name)} -o {str(tmp_sorted.name)}"
        )
        subprocess.check_call(split(cmd_sort))
        with open(output, "wb") as f:
            cmd_zip = f"bgzip -c {str(tmp_sorted.name)}"
            subprocess.check_call(split(cmd_zip), stdout=f)
        cmd_index = f"tabix -f -0 -p vcf {str(output)}"
        subprocess.check_call(split(cmd_index))


# %%


def check_if_all_svtypes_are_supported(sv_types: list[str]) -> None:
    unsupported_sv_types = []
    for svtype in sv_types:
        if svtype not in SUPPORTED_SV_TYPE_STRINGS:
            unsupported_sv_types.append(svtype)
    if len(unsupported_sv_types) > 0:
        raise ValueError(
            f"Unsupported SV types: {', '.join(unsupported_sv_types)}. Supported SV types are: {', '.join(SUPPORTED_SV_TYPE_STRINGS)}"
        )


def load_copynumber_tracks_from_svirltiles(
    svirltile_paths: list[Path | str], samplenames: list[str]
) -> dict[str, dict[str, IntervalTree]]:
    """
    Load copy number tracks from svirltile databases.

    Args:
        svirltile_paths: List of paths to svirltile databases
        samplenames: List of sample names corresponding to each database

    Returns:
        Dictionary mapping samplename -> chromosome -> IntervalTree with CN data
    """
    from ..signalprocessing.copynumber_tracks import \
        load_copynumber_trees_from_db

    cn_tracks = {}
    for _i, (path, samplename) in enumerate(
        zip(svirltile_paths, samplenames, strict=True)
    ):
        try:
            cn_tracks[samplename] = load_copynumber_trees_from_db(Path(path))
            log.info(f"Loaded copy number tracks for sample {samplename} from {path}")
        except Exception as e:
            log.warning(
                f"Could not load copy number tracks for sample {samplename} from {path}: {e}"
            )
            log.warning(f"Using empty copy number tracks for {samplename}")
            cn_tracks[samplename] = {}

    return cn_tracks


def create_dummy_covtrees_from_reference(
    reference: Path, samplenames: list[str], default_coverage: int = 30
) -> dict[str, dict[str, IntervalTree]]:
    """
    Create dummy coverage trees with uniform coverage for all samples.
    Reads chromosome lengths from reference .fai file.

    Args:
        reference: Path to reference genome file (will look for .fai file)
        samplenames: List of sample names to create covtrees for
        default_coverage: Uniform coverage value to use (default: 30)

    Returns:
        dict mapping samplename to dict of chr_name to IntervalTree with uniform coverage
    """
    # Find the .fai file - handle various reference extensions
    ref_path = Path(reference)

    # Try common patterns for finding the .fai file
    if ref_path.suffix in [".mmi", ".fa", ".fasta"]:
        # Replace extension with .fa.fai or try adding .fai
        fai_candidates = [
            ref_path.with_suffix(".fa.fai"),
            ref_path.with_suffix(".fasta.fai"),
            Path(str(ref_path) + ".fai"),
        ]
    else:
        # Just add .fai
        fai_candidates = [Path(str(ref_path) + ".fai")]

    fai_path = None
    for candidate in fai_candidates:
        if candidate.exists():
            fai_path = candidate
            break

    if fai_path is None:
        raise FileNotFoundError(
            f"Could not find .fai index file for reference {reference}. Tried: {fai_candidates}"
        )

    log.info(f"Reading chromosome lengths from {fai_path}")

    # Parse the .fai file to get chromosome names and lengths
    chr_lengths = {}
    with open(fai_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chr_name = parts[0]
                chr_length = int(parts[1])
                chr_lengths[chr_name] = chr_length

    log.info(
        f"Creating dummy covtrees with uniform coverage of {default_coverage} for {len(samplenames)} samples and {len(chr_lengths)} chromosomes"
    )

    # Create covtrees for all samples
    covtrees = {}
    for samplename in samplenames:
        sample_covtree = {}
        for chr_name, chr_length in chr_lengths.items():
            # Create a single interval covering the entire chromosome with uniform coverage
            tree = IntervalTree()
            tree.addi(0, chr_length, default_coverage)
            sample_covtree[chr_name] = tree
        covtrees[samplename] = sample_covtree

    return covtrees


def multisample_sv_calling(
    input: list[Path | str],
    output: Path,
    reference: Path,
    threads: int,
    max_cohens_d: float,
    near: int,
    min_kmer_overlap: float,
    sv_types: list[str],
    min_sv_size: int,
    symbolic_threshold: int,
    apriori_size_difference_fraction_tolerance: float,
    find_leftmost_reference_position: bool,
    verbose: bool = False,
    tmp_dir_path: Path | str | None = None,
    candidate_regions_file: Path | str | None = None,
    skip_covtrees: bool = False,
) -> None:
    check_if_all_svtypes_are_supported(sv_types=sv_types)
    samplenames = [svirltile.get_metadata(Path(path))["samplename"] for path in input]

    # Create covtrees - either real or dummy uniform coverage
    if skip_covtrees:
        log.info("Skipping covtree computation, using uniform coverage of 30")
        covtrees: dict[str, dict[str, IntervalTree]] = (
            create_dummy_covtrees_from_reference(
                reference=reference, samplenames=samplenames, default_coverage=30
            )
        )
    else:
        log.info("Computing covtrees from sample data")
        covtrees: dict[str, dict[str, IntervalTree]] = {
            samplenames[i]: covtree(path_db=input[i]) for i in range(len(input))
        }

    # Load copy number tracks from svirltile databases
    log.info("Loading copy number tracks from svirltile databases")
    cn_tracks: dict[str, dict[str, IntervalTree]] = (
        load_copynumber_tracks_from_svirltiles(
            svirltile_paths=input, samplenames=samplenames
        )
    )

    # Parse candidate regions file if provided
    candidate_regions_filter: dict[str, set[int]] | None = None
    if candidate_regions_file is not None:
        candidate_regions_filter = parse_candidate_regions_file(
            Path(candidate_regions_file)
        )
        log.info(
            f"Filtering SVpatterns to candidate regions from {candidate_regions_file}."
        )

    # Convert string sv_types to SVpattern types
    sv_types_set: set[type[SVpatterns.SVpatternType]] = {
        SUPPORTED_SV_TYPE_STRINGS_INVERSE[sv]
        for sv in sv_types
        if sv in SUPPORTED_SV_TYPE_STRINGS_INVERSE
    }
    # debug - check what types are in sv_types_set
    if verbose:
        log.info(
            f"SV types to be processed: {', '.join([t.__name__ for t in sv_types_set])}"
        )

    data: list[SVcomposite] = generate_svComposites_from_dbs(
        input=input,
        sv_types=sv_types_set,
        candidate_regions_filter=candidate_regions_filter,
    )
    if verbose:
        svtype_counts: dict[str, int] = {}
        for svComposite in data:
            svtype = svComposite.sv_type.__name__
            if svtype not in svtype_counts:
                svtype_counts[svtype] = 0
            svtype_counts[svtype] += 1
        log.info(f"SVcomposite counts by type before merging: {svtype_counts}")

    # if tmp dir is provided, dump all svComposites to a compressed json file
    if tmp_dir_path is not None:
        save_svComposites_to_json(
            data=data, output_path=Path(tmp_dir_path) / "all_svComposites.json.gz"
        )

    # --- vertical merging of svComposites across samples and consensus sequences --- #
    merged: list[SVcomposite] = merge_svComposites(
        apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
        svComposites=data,
        max_cohens_d=max_cohens_d,
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        threads=threads,
        verbose=verbose,
    )

    if verbose:
        svtype_counts_merged: dict[str, int] = {}
        for svComposite in merged:
            svtype = svComposite.sv_type.__name__
            if svtype not in svtype_counts_merged:
                svtype_counts_merged[svtype] = 0
            svtype_counts_merged[svtype] += 1
        log.info(f"SVcomposite counts by type after merging: {svtype_counts_merged}")

    merged = [
        svComposite
        for svComposite in merged
        if abs(svComposite.get_size()) >= min_sv_size
    ]

    # if tmp dir is provided, dump all merged svComposites to a compressed json file
    if tmp_dir_path is not None:
        save_svComposites_to_json(
            data=merged, output_path=Path(tmp_dir_path) / "merged_svComposites.json.gz"
        )
    svCalls: list[SVcall] = [
        svcall
        for svComposite in merged
        for svcall in SVcalls_from_SVcomposite(
            svComposite,
            covtrees=covtrees,
            cn_tracks=cn_tracks,
            find_leftmost_reference_position=find_leftmost_reference_position,
        )
    ]
    if verbose:
        # check if there are SVcalls from the consensus with ID 15.0 or 15.1
        for svCall in svCalls:
            if "15.0" in svCall.consensusIDs or "15.1" in svCall.consensusIDs:
                log.info(
                    f"   ------> Found SVcall with consensusID {svCall.consensusIDs} and type {svCall.svtype} and location {svCall.chrname}:{svCall.start}-{svCall.end}"
                )

    ref_bases_dict: dict[str, str] = reference_bases_by_merged_svComposites(
        svComposites=merged,
        reference=reference,
        find_leftmost_reference_position=find_leftmost_reference_position,
    )
    write_svCalls_to_vcf(
        svCalls=svCalls,
        output=output,
        reference=reference,
        samplenames=samplenames,
        covtrees=covtrees,
        refdict=ref_bases_dict,
        symbolic_threshold=symbolic_threshold,
    )


def run(args) -> None:
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s - %(levelname)s - %(message)s",
        force=True,
    )

    multisample_sv_calling(
        input=args.input,
        output=args.output,
        reference=args.reference,
        threads=args.threads,
        max_cohens_d=args.max_cohens_d,
        near=args.near,
        min_kmer_overlap=args.min_kmer_overlap,
        sv_types=args.sv_types,
        min_sv_size=args.min_sv_size,
        apriori_size_difference_fraction_tolerance=args.apriori_size_difference_fraction_tolerance,
        symbolic_threshold=args.symbolic_threshold,
        find_leftmost_reference_position=args.find_leftmost_reference_position,
        tmp_dir_path=args.tmp_dir_path,
        verbose=args.verbose,
        candidate_regions_file=args.candidate_regions_file,
        skip_covtrees=args.skip_covtrees,
    )


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--input",
        help="Paths to the per sample svirltiles.",
        nargs="+",
        required=True,
        type=os.path.abspath,
    )
    parser.add_argument(
        "--output", help="Path to output VCF file.", required=True, type=os.path.abspath
    )
    parser.add_argument(
        "--reference",
        help="Path to reference genome file.",
        required=True,
        type=os.path.abspath,
    )
    parser.add_argument(
        "--threads", help="Number of threads to use (default: 4).", type=int, default=1
    )
    parser.add_argument(
        "--sv-types",
        help=f"List of structural variant types to include. Default: all supported types. Allowed are {SUPPORTED_SV_TYPE_STRINGS_LIST}.",
        nargs="+",
        default=SUPPORTED_SV_TYPE_STRINGS_LIST,
    )
    parser.add_argument(
        "--max_cohens_d",
        help="Maximum Cohen's d value for merging SVs (default: 2.0).",
        type=float,
        default=2.0,
    )
    parser.add_argument(
        "--near",
        help="Maximum distance for merging SVs (default: 150).",
        type=int,
        default=150,
    )
    parser.add_argument(
        "--min_kmer_overlap",
        help="Minimum k-mer overlap for merging SVs (default: 0.7).",
        type=float,
        default=0.7,
    )
    parser.add_argument(
        "--min-sv-size",
        help="Minimum SV size to include (default: 50).",
        type=int,
        default=50,
    )
    parser.add_argument(
        "--apriori-size-difference-fraction-tolerance",
        help="Size difference fraction tolerance for merging SVs (default: 0.2). Decrease for stronger separation of haplotypes",
        type=float,
        default=0.2,
    )
    parser.add_argument(
        "--symbolic-threshold",
        help="Sequence length threshold for using symbolic alleles in VCF (default: 100000). Sequences longer than this will be written to a companion FASTA file.",
        type=int,
        default=100000,
    )
    parser.add_argument(
        "--find-leftmost-reference-position",
        help="When determining the reference position of an SVcomposite, use the leftmost position of all underlying SVpatterns instead of the one with most supporting reads * size. This might be better aligned with the giab SV benchmark, but should be discussed in the paper.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--tmp-dir-path",
        help="Path to temporary directory (default: system temp dir).",
        type=os.path.abspath,
        default=None,
    )
    parser.add_argument(
        "--candidate-regions-file",
        help="Optional TSV file with samplename and comma-separated candidate region IDs (crID) to filter SVpatterns. Format: samplename<TAB>crID1,crID2,crID3",
        type=os.path.abspath,
        default=None,
    )
    parser.add_argument(
        "--skip-covtrees",
        help="Skip computing coverage trees from sample data. Instead, use uniform coverage of 30 across all chromosomes. This significantly speeds up execution when genotype coverage information is not critical.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--log-level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument(
        "--verbose", help="Enable verbose output.", action="store_true", default=False
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Multiple sample SV calling from precomputed svirltile files per sample."
    )
    add_arguments(parser)
    return parser


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================


def main():
    """Main entry point for the consensus script."""
    parser = get_parser()
    args = parser.parse_args()

    # Convert string log level to logging constant
    log_level = getattr(logging, args.log_level.upper())

    if args.logfile:
        # Configure logging to file
        logging.basicConfig(
            level=log_level,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            handlers=[logging.FileHandler(str(args.logfile)), logging.StreamHandler()],
        )
    else:
        # Configure logging to console only
        logging.basicConfig(
            level=log_level,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            handlers=[logging.StreamHandler()],
        )
    run(args)
    return


if __name__ == "__main__":
    main()


# %%


def save_svComposites_to_json(data: list[SVcomposite], output_path: Path | str) -> None:
    """
    Save SVcomposites to a JSON file.

    Compression is automatically detected from file extension:
    - .json.gz: compressed with gzip
    - .json: uncompressed

    Args:
        data: List of SVcomposites to serialize
        output_path: Path to output JSON file

    Returns:
        None
    """
    output_path = Path(output_path)

    # Create parent directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Auto-detect compression from extension
    compress = str(output_path).endswith(".gz")

    # Serialize all SVcomposites
    serialized_data = [sv.unstructure() for sv in data]

    # Write to file
    if compress:
        # Ensure filename ends with .json.gz
        if not str(output_path).endswith(".json.gz"):
            output_path = Path(str(output_path).replace(".json", "") + ".json.gz")

        with gzip.open(output_path, "wt", encoding="utf-8") as f:
            json.dump(serialized_data, f, indent=2)
        log.info(f"Saved {len(data)} SVcomposites to compressed JSON: {output_path}")
    else:
        # Ensure filename ends with .json
        if str(output_path).endswith(".gz"):
            output_path = Path(str(output_path).replace(".gz", ""))
        if not str(output_path).endswith(".json"):
            output_path = Path(str(output_path) + ".json")

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(serialized_data, f, indent=2)
        log.info(f"Saved {len(data)} SVcomposites to JSON: {output_path}")


def load_svComposites_from_json(input_path: Path | str) -> list[SVcomposite]:
    """
    Load SVcomposites from a JSON file.

    Compression is automatically detected from file extension:
    - .json.gz: decompressed from gzip
    - .json: read as plain text

    Args:
        input_path: Path to input JSON file

    Returns:
        List of SVcomposites
    """
    input_path = Path(input_path)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # Auto-detect compression from extension
    compressed = str(input_path).endswith(".gz")

    try:
        if compressed:
            with gzip.open(input_path, "rt", encoding="utf-8") as f:
                serialized_data = json.load(f)
            log.info(
                f"Loaded {len(serialized_data)} SVcomposites from compressed JSON: {input_path}"
            )
        else:
            with open(input_path, "r", encoding="utf-8") as f:
                serialized_data = json.load(f)
            log.info(
                f"Loaded {len(serialized_data)} SVcomposites from JSON: {input_path}"
            )

        # Deserialize SVcomposites
        svcomposites = [
            SVpatterns.converter.structure(item, SVcomposite)
            for item in serialized_data
        ]

        return svcomposites

    except json.JSONDecodeError as e:
        log.error(f"Failed to parse JSON from {input_path}: {e}")
        raise
    except Exception as e:
        log.error(f"Error loading SVcomposites from {input_path}: {e}")
        raise


def extract_test_svComposites(
    data: list[SVcomposite],
    consensus_ids: list[str],
    sv_types: list[str] | None = None,
    output_path: str | Path | None = None,
) -> list[SVcomposite]:
    # Filter by consensus IDs
    filtered_svComposites = [
        svComposite
        for svComposite in data
        if {svp.consensusID for svp in svComposite.svPatterns}.intersection(
            set(consensus_ids)
        )
    ]
    # Optional filter by SV types
    if sv_types is not None:
        # Convert string sv_types to pattern type set for filtering
        sv_types_set = {
            SUPPORTED_SV_TYPE_STRINGS_INVERSE[sv]
            for sv in sv_types
            if sv in SUPPORTED_SV_TYPE_STRINGS_INVERSE
        }
        filtered_svComposites = [
            svComposite
            for svComposite in filtered_svComposites
            if svComposite.sv_type in sv_types_set
        ]
    # Save to file if output_path is provided
    if output_path is not None:
        output_path = Path(output_path)
        # Create parent directory if it doesn't exist
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save the data
        with open(output_path, "wb") as f:
            pickle.dump(filtered_svComposites, f)

        print(f"Saved {len(filtered_svComposites)} SVcomposites to {output_path}")

    return filtered_svComposites


# #%%
# # search in the data all variants that come from the consensus with ID 19.0 and 19.2


# test_svComposites = extract_test_svComposites(
#     data=data,
#     consensus_ids=["19.0", "19.2"],
#     output_path= BASE_PATH / "test_vertically_merged_svComposites_from_group__large_dels.json.gz")

# test_insertions = extract_test_svComposites(
#     data=data,
#     consensus_ids=["0.0", "0.1", "0.2", "0.3", "0.4", "0.5"],
#     sv_types=["INS"],
#     output_path=BASE_PATH / "test_vertically_merged_svComposites_from_group__insertions.json.gz")

# test_svComposites = extract_test_svComposites(
#     data=data,
#     consensus_ids=["1.0", "1.1", "1.2", "1.3"],
#     sv_types=["DEL"],
#     output_path=BASE_PATH / "test_vertically_merged_svComposites_from_group__small_dels.json.gz")

# test_svComposites = extract_test_svComposites(
#     data=data,
#     consensus_ids=[f"6.{i}" for i in range(10)],
#     sv_types=["INS"],
#     output_path=BASE_PATH / "test_vertically_merged_svComposites_from_group__many_ins.json.gz")
