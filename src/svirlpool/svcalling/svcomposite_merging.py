"""Functions for merging and deciding whether two SVcomposites should be merged.

This module contains the vertical merging logic:
- can_merge_svComposites_* functions: decide if two SVcomposites can be merged
- merge_* functions: merge lists of SVcomposites by SV type
- merge_svComposites: top-level entry point for merging across samples/chromosomes
"""

import logging
import math
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from intervaltree import IntervalTree
from tqdm import tqdm

from ..localassembly import SVpatterns
from ..util.datastructures import UnionFind
from ..util.util import kmer_similarity_of_groups
from . import merge_audit
from .SVcomposite import SVcomposite
from .svcomposite_utils import (
    _crIDs_from_svcomposite,
    _regions_str_from_svcomposite,
    _svcomposite_log_id,
    _svcomposite_short_id,
)

log = logging.getLogger(__name__)


#: The vertical merge's size test, set once from the CLI before any worker
#: process is forked. Two sizes are similar when
#:
#:     |size_a - size_b| <= max(SIZE_TOLERANCE_FLOOR, fraction * ref(size_a, size_b))
#:
#: with ref the smaller size ('min'), the larger ('max') or their harmonic mean
#: ('hmean'), and `fraction` the caller's tolerance
#: (--apriori-size-difference-fraction-tolerance, or SIBLING_SIZE_TOLERANCE for
#: two haplotypes of one sample). Nothing else enters the size decision.
#:
#: Above the floor the three references are one test in three
#: parametrisations: each accepts exactly the pairs with max/min <= R, where
#: R = 1 + f for 'min', 1 / (1 - f) for 'max' and f + sqrt(f^2 + 1) for
#: 'hmean'. They differ only in which ratio a given f stands for.
#:
#: Two further allowances were removed after the svp_merging size-gate study
#: (2026-09-29): a sequence-complexity allowance (1 - complexity) * |size|,
#: which rejected identical sizes and was unbounded, and a "noise" gate, Cohen's
#: d of the per-read distortion populations, which never improved the trio
#: benchmark and whose scale overstates the consensus' size error 2-8x.
SIZE_TOLERANCE_REFERENCES: tuple[str, ...] = ("min", "max", "hmean")
SIZE_TOLERANCE_REFERENCE: str = "max"
SIZE_TOLERANCE_FLOOR: float = 0.0


def size_tolerance(size_a: float, size_b: float, fraction: float) -> float:
    """The largest |size_a - size_b| at which the two sizes still count as similar."""
    x, y = abs(size_a), abs(size_b)
    if SIZE_TOLERANCE_REFERENCE == "min":
        ref = min(x, y)
    elif SIZE_TOLERANCE_REFERENCE == "max":
        ref = max(x, y)
    elif SIZE_TOLERANCE_REFERENCE == "hmean":
        ref = 2.0 * x * y / (x + y) if x + y > 0 else 0.0
    else:
        raise ValueError(
            f"SIZE_TOLERANCE_REFERENCE must be one of {SIZE_TOLERANCE_REFERENCES}, "
            f"got {SIZE_TOLERANCE_REFERENCE!r}"
        )
    return max(SIZE_TOLERANCE_FLOOR, fraction * ref)


def sizes_similar(size_a: float, size_b: float, fraction: float) -> bool:
    """Whether two SV sizes agree within the floored relative tolerance."""
    return abs(abs(size_a) - abs(size_b)) <= size_tolerance(size_a, size_b, fraction)


def _similar_size_detail(
    a: SVcomposite,
    b: SVcomposite,
    apriori_size_difference_fraction_tolerance: float,
) -> dict:
    """The size test of the vertical merge, with its inputs for --merge-audit.

    Keys: similar, size_tolerance (bp), size_reference, size_floor.
    """
    size_a = a.get_size()
    size_b = b.get_size()
    tolerance = size_tolerance(
        size_a, size_b, apriori_size_difference_fraction_tolerance
    )
    return {
        "similar": abs(abs(size_a) - abs(size_b)) <= tolerance,
        "size_tolerance": tolerance,
        "size_reference": SIZE_TOLERANCE_REFERENCE,
        "size_floor": SIZE_TOLERANCE_FLOOR,
    }


def _similar_size(
    a: SVcomposite,
    b: SVcomposite,
    apriori_size_difference_fraction_tolerance: float,
) -> bool:
    """Whether the sizes of two SVcomposites are similar (see `size_tolerance`)."""
    return _similar_size_detail(a, b, apriori_size_difference_fraction_tolerance)[
        "similar"
    ]


#: Experimental, set once from the CLI before any worker process is forked.
#: When on (the default), vertical merging treats the consensus assemblies of one sample and
#: candidate region as the haplotypes they are: two patterns of the SAME assembly
#: are distinct events and never merge; two patterns of SIBLING assemblies (same
#: sample, same crID, different consensus) are two haplotypes' alleles and merge
#: only if their sizes agree within SIBLING_SIZE_TOLERANCE (see `size_tolerance`).
HAPLOTYPE_AWARE_MERGE: bool = True
SIBLING_SIZE_TOLERANCE: float = 0.1


def _assembly_relation(a: SVcomposite, b: SVcomposite) -> str:
    """'same_assembly', 'sibling_assembly' or 'independent'."""
    cons_a = {(p.samplename, p.consensusID) for p in a.svPatterns}
    cons_b = {(p.samplename, p.consensusID) for p in b.svPatterns}
    if cons_a & cons_b:
        return "same_assembly"
    crs_a = {(s, c.split(".")[0]) for s, c in cons_a}
    crs_b = {(s, c.split(".")[0]) for s, c in cons_b}
    if crs_a & crs_b:
        return "sibling_assembly"
    return "independent"


def _audit_value(v):
    """JSON-safe: numpy scalars to Python, nan to None."""
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return None if math.isnan(v) else float(v)
    return v


def _audit_pair(
    svtype: str, a: SVcomposite, b: SVcomposite, decision: bool, reason: str, **features
) -> None:
    """One vertical-merge decision for --merge-audit; a no-op when it is off."""
    if not merge_audit.enabled():
        return
    merge_audit.record(
        stage="vertical",
        svtype=svtype,
        relation=_assembly_relation(a, b),
        decision=decision,
        reason=reason,
        chr=a.ref_start[0],
        pos_a=int(a.ref_start[1]),
        pos_b=int(b.ref_start[1]),
        size_a=int(a.get_size()),
        size_b=int(b.get_size()),
        patterns_a=merge_audit.pattern_keys(a),
        patterns_b=merge_audit.pattern_keys(b),
        **{k: _audit_value(v) for k, v in features.items()},
    )


def _audit_size_features(size_detail: dict) -> dict:
    """The size test's tolerance and its inputs, for --merge-audit only."""
    if not merge_audit.enabled():
        return {}
    return {k: v for k, v in size_detail.items() if k != "similar"}


def _haplotype_gate(a: SVcomposite, b: SVcomposite, svtype: str) -> bool | None:
    """Decide the pair outright under HAPLOTYPE_AWARE_MERGE, or return None."""
    if not HAPLOTYPE_AWARE_MERGE:
        return None
    relation = _assembly_relation(a, b)
    if relation == "same_assembly":
        log.debug(
            f"VERTICAL_MERGE|{svtype}|REJECT_SAME_ASSEMBLY	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}"
        )
        _audit_pair(svtype, a, b, False, "same_assembly")
        return False
    if relation == "sibling_assembly":
        size_a, size_b = abs(a.get_size()), abs(b.get_size())
        ok = sizes_similar(size_a, size_b, SIBLING_SIZE_TOLERANCE)
        log.debug(
            f"VERTICAL_MERGE|{svtype}|{'ACCEPT' if ok else 'REJECT_SIBLING_SIZE'}	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"size_a={size_a}	size_b={size_b}	sibling_tolerance={SIBLING_SIZE_TOLERANCE}	"
            f"size_reference={SIZE_TOLERANCE_REFERENCE}	size_floor={SIZE_TOLERANCE_FLOOR}"
        )
        if not ok:
            _audit_pair(
                svtype,
                a,
                b,
                False,
                "sibling_size",
                sibling_tolerance=SIBLING_SIZE_TOLERANCE,
                size_tolerance=size_tolerance(size_a, size_b, SIBLING_SIZE_TOLERANCE),
                size_reference=SIZE_TOLERANCE_REFERENCE,
                size_floor=SIZE_TOLERANCE_FLOOR,
            )
            return False
        # sizes agree: the remaining checks (proximity, k-mers) still apply
    return None


def can_merge_svComposites_insertions(
    a: SVcomposite,
    b: SVcomposite,
    apriori_size_difference_fraction_tolerance: float,
    near: int,
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

    _gate = _haplotype_gate(a, b, "INS")
    if _gate is False:
        return False
    _sibling = HAPLOTYPE_AWARE_MERGE and _assembly_relation(a, b) == "sibling_assembly"

    are_near: bool = a.overlaps_any(
        b, tolerance_radius=max(near, max(a.get_size(), b.get_size()))
    )  # either spatially close or share at least one repeatID
    if not are_near:
        log.debug(
            f"VERTICAL_MERGE|INS|REJECT_NOT_NEAR	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"regions_a={_regions_str_from_svcomposite(a)}	regions_b={_regions_str_from_svcomposite(b)}	tolerance={near}"
        )
        if verbose:
            print(
                f"Cannot merge insertions: SVcomposites are not near (tolerance_radius={near}) ({a.get_regions()} vs {b.get_regions()})"
            )
        _audit_pair("INS", a, b, False, "not_near", near=near)
        return False

    # - have similar sizes (floored relative tolerance, see `size_tolerance`)

    size_a = a.get_size()
    size_b = b.get_size()

    _size = _similar_size_detail(a, b, apriori_size_difference_fraction_tolerance)
    similar_size = _size["similar"]

    if verbose:
        print(
            f"Size comparison: size_a={size_a}, size_b={size_b}, "
            f"tolerance={_size['size_tolerance']:.1f} bp, similar_size={similar_size}"
        )

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

    _features = {
        **_audit_size_features(_size),
        "kmer_sim": similarity,
        "kmer_compared": not (size_a < 100 and size_b < 100),
    }
    if _sibling:
        similar_size = True  # decided by _haplotype_gate
    if not similar_size:
        _audit_pair("INS", a, b, False, "size", **_features)
        log.debug(
            f"VERTICAL_MERGE|INS|REJECT_SIZE	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"size_a={size_a}	size_b={size_b}	regions_a={_regions_str_from_svcomposite(a)}	regions_b={_regions_str_from_svcomposite(b)}"
        )
        if verbose:
            print("Cannot merge insertions: sizes are not similar enough")
        return False
    if not similar_insertion_kmers:
        log.debug(
            f"VERTICAL_MERGE|INS|REJECT_KMER	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"kmer_sim={similarity:.3f}	threshold={min_kmer_overlap}	regions_a={_regions_str_from_svcomposite(a)}	regions_b={_regions_str_from_svcomposite(b)}"
        )
        if verbose:
            print("Cannot merge insertions: k-mer similarity is too low")
        _audit_pair("INS", a, b, False, "kmer", **_features)
        return False

    log.debug(
        f"VERTICAL_MERGE|INS|ACCEPT	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
        f"size_a={size_a}	size_b={size_b}	kmer_sim={similarity:.3f}	regions_a={_regions_str_from_svcomposite(a)}	regions_b={_regions_str_from_svcomposite(b)}"
    )
    if verbose:
        print("Can merge insertions: all criteria passed")
    _audit_pair("INS", a, b, True, "accept", **_features)
    return True


def can_merge_svComposites_deletions(
    a: SVcomposite,
    b: SVcomposite,
    apriori_size_difference_fraction_tolerance: float,
    near: int,
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

    _gate = _haplotype_gate(a, b, "DEL")
    if _gate is False:
        return False
    _sibling = HAPLOTYPE_AWARE_MERGE and _assembly_relation(a, b) == "sibling_assembly"

    are_near: bool = a.overlaps_any(
        b, tolerance_radius=near
    )  # eiter spatially close or share at least one repeatID
    if not are_near:
        log.debug(
            f"VERTICAL_MERGE|DEL|REJECT_NOT_NEAR	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"regions_a={_regions_str_from_svcomposite(a)}|regions_b={_regions_str_from_svcomposite(b)}|tolerance={near}"
        )
        if verbose:
            print(
                f"Cannot merge deletions: SVcomposites are not near (tolerance_radius={near}) ({a.get_regions()} vs {b.get_regions()})"
            )
        _audit_pair("DEL", a, b, False, "not_near", near=near)
        return False

    size_a = a.get_size()
    size_b = b.get_size()

    _size = _similar_size_detail(a, b, apriori_size_difference_fraction_tolerance)
    similar_size = _size["similar"]

    if verbose:
        print(
            f"Size comparison: size_a={size_a}, size_b={size_b}, "
            f"tolerance={_size['size_tolerance']:.1f} bp, similar_size={similar_size}"
        )

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

    _features = {
        **_audit_size_features(_size),
        "kmer_sim": similarity,
        "kmer_compared": not (size_a < 100 and size_b < 100),
    }
    if _sibling:
        similar_size = True  # decided by _haplotype_gate
    if not similar_size:
        _audit_pair("DEL", a, b, False, "size", **_features)
        log.debug(
            f"VERTICAL_MERGE|DEL|REJECT_SIZE	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"size_a={size_a}|size_b={size_b}|regions_a={_regions_str_from_svcomposite(a)}|regions_b={_regions_str_from_svcomposite(b)}"
        )
        if verbose:
            print("Cannot merge deletions: sizes are not similar enough")
        return False
    if not similar_deletion_kmers:
        log.debug(
            f"VERTICAL_MERGE|DEL|REJECT_KMER	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"kmer_sim={similarity:.3f}|threshold={min_kmer_overlap}|regions_a={_regions_str_from_svcomposite(a)}|regions_b={_regions_str_from_svcomposite(b)}"
        )
        if verbose:
            print("Cannot merge deletions: k-mer similarity is too low")
        _audit_pair("DEL", a, b, False, "kmer", **_features)
        return False

    log.debug(
        f"VERTICAL_MERGE|DEL|ACCEPT	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
        f"size_a={size_a}|size_b={size_b}|kmer_sim={similarity:.3f}|regions_a={_regions_str_from_svcomposite(a)}|regions_b={_regions_str_from_svcomposite(b)}"
    )
    if verbose:
        print("Can merge deletions: all criteria passed")
    _audit_pair("DEL", a, b, True, "accept", **_features)
    return True


def vertically_merged_svComposites_from_group(
    group: list[SVcomposite],
    near: int,
    min_kmer_overlap: float,
    apriori_size_difference_fraction_tolerance: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    # The group is divided into subgroups with shared SV types
    all_crIDs: set[int] = set()
    for svc in group:
        all_crIDs.update(_crIDs_from_svcomposite(svc))
    log.debug(
        "VERTICAL_MERGE|GROUP_START|group_size=%d|crIDs=%s",
        len(group),
        sorted(all_crIDs),
    )
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
    adjacencies = {
        i: svc
        for i, svc in enumerate(group)
        if issubclass(svc.sv_type, SVpatterns.SVpatternAdjacency)
    }
    log.debug(
        "VERTICAL_MERGE|GROUP_CATEGORIZE|INS=%d|DEL=%d|INV=%d|BND=%d|ADJ=%d|crIDs=%s",
        len(insertions),
        len(deletions),
        len(inversions),
        len(breakends),
        len(adjacencies),
        sorted(all_crIDs),
    )

    used_indices = set().union(
        insertions.keys(),
        deletions.keys(),
        inversions.keys(),
        breakends.keys(),
        adjacencies.keys(),
    )
    others = [
        svc for i, svc in enumerate(group) if i not in used_indices
    ]  # BND, DUP, etc.
    if others:
        other_types: dict[str, int] = defaultdict(int)
        for svc in others:
            other_types[svc.sv_type.get_sv_type()] += 1
        log.debug(
            "VERTICAL_MERGE|GROUP_OTHERS|counts=%s|crIDs=%s",
            dict(other_types),
            sorted(all_crIDs),
        )

    merged_insertions = merge_insertions(
        insertions=list(insertions.values()),
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        verbose=verbose,
        apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
    )
    merged_deletions = merge_deletions(
        deletions=list(deletions.values()),
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        verbose=verbose,
        apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
    )
    merged_inversions = merge_inversions(
        inversions=list(inversions.values()),
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        verbose=verbose,
        apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
    )
    merged_breakends = merge_breakends(
        breakends=list(breakends.values()),
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        verbose=verbose,
        apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
    )
    merged_cpx = merge_adjacencies(
        adjacencies=list(adjacencies.values()),
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        verbose=verbose,
    )

    result = (
        merged_insertions
        + merged_deletions
        + merged_inversions
        + merged_breakends
        + merged_cpx
        + others
    )
    log.debug(
        "VERTICAL_MERGE|GROUP_DONE|in=%d|out=%d|crIDs=%s",
        len(group),
        len(result),
        sorted(all_crIDs),
    )
    return result


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
    log.debug(f"MERGE_INSERTIONS|START	n_composites={len(insertions)}")
    uf = UnionFind(range(len(insertions)))
    for i in range(len(insertions)):
        for j in range(i + 1, len(insertions)):
            if can_merge_svComposites_insertions(
                a=insertions[i],
                b=insertions[j],
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                verbose=verbose,
            ):
                uf.union(i, j)
    result: list[SVcomposite] = []
    for cc in uf.get_connected_components(
        allow_singletons=True
    ):  # connected components of svComposites to merge
        svpatterns_to_merge = [
            svPattern for idx in cc for svPattern in insertions[idx].svPatterns
        ]
        old_svcs = [insertions[idx] for idx in cc]
        new_svc = SVcomposite.from_SVpatterns(svPatterns=svpatterns_to_merge)
        if len(cc) > 1:
            log.debug(
                f"TRANSFORMED::merge_insertions::vertical_merge:(to merged SVcomposite) "
                f"svComposites={[_svcomposite_log_id(svc) for svc in old_svcs]} "
                f"-->   {_svcomposite_log_id(new_svc)}"
            )
        result.append(new_svc)
    log.debug(
        f"MERGE_INSERTIONS|DONE	input={len(insertions)}	output={len(result)}"
    )
    return result


def merge_deletions(
    deletions: list[SVcomposite],
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
    log.debug(f"MERGE_DELETIONS|START	n_composites={len(deletions)}")
    uf = UnionFind(range(len(deletions)))
    for i in range(len(deletions)):
        for j in range(i + 1, len(deletions)):
            if verbose:
                print(
                    f"Checking deletions {i} and {j}: {deletions[i].svPatterns[0].samplename}:{deletions[i].svPatterns[0].consensusID} vs {deletions[j].svPatterns[0].samplename}:{deletions[j].svPatterns[0].consensusID}"
                )
            if can_merge_svComposites_deletions(
                a=deletions[i],
                b=deletions[j],
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
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
        old_svcs = [deletions[idx] for idx in cc]
        new_svc = SVcomposite.from_SVpatterns(svPatterns=svPatterns_to_merge)
        if len(cc) > 1:
            log.debug(
                f"TRANSFORMED::merge_deletions::vertical_merge:(to merged SVcomposite)"
                f"svComposites={[_svcomposite_log_id(svc) for svc in old_svcs]} "
                f"-->   {_svcomposite_log_id(new_svc)}"
            )
        else:
            log.debug(
                f"TRANSFORMED::merge_deletions::vertical_merge: (no change): {_svcomposite_log_id(new_svc)}"
            )
        result.append(new_svc)
    log.debug(f"MERGE_DELETIONS|DONE	input={len(deletions)}	output={len(result)}")
    return result


def can_merge_svComposites_inversions(
    a: SVcomposite,
    b: SVcomposite,
    apriori_size_difference_fraction_tolerance: float,
    near: int,
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

    are_near: bool = a.overlaps_any(
        b, tolerance_radius=near
    )  # eiter spatially close or share at least one repeatID
    if not are_near:
        log.debug(
            f"VERTICAL_MERGE|INV|REJECT_NOT_NEAR	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"regions_a={_regions_str_from_svcomposite(a)}	regions_b={_regions_str_from_svcomposite(b)}	tolerance={near}"
        )
        return False

    if (
        a.get_representative_SVpattern().SVprimitives[0].chr
        != b.get_representative_SVpattern().SVprimitives[0].chr
        or a.get_representative_SVpattern().SVprimitives[-1].chr
        != b.get_representative_SVpattern().SVprimitives[-1].chr
    ):
        log.debug(
            f"VERTICAL_MERGE|INV|REJECT_DIFF_CHR	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"chr_a=({a.get_representative_SVpattern().SVprimitives[0].chr},{a.get_representative_SVpattern().SVprimitives[-1].chr})	"
            f"chr_b=({b.get_representative_SVpattern().SVprimitives[0].chr},{b.get_representative_SVpattern().SVprimitives[-1].chr})"
        )
        return False

    # in the case that inversions of different subclasses are merged, it is necessary to measure their inner sizes to find a good size comparison,
    # since outer break ends can align very well while the inverted interval can be of very different sizes.
    # this allows to merge e.g. inverted dels with inverted dups given a stronger local noise.
    size_a = a.get_size()
    size_b = b.get_size()

    _size = _similar_size_detail(a, b, apriori_size_difference_fraction_tolerance)
    similar_size = _size["similar"]

    if verbose:
        print(
            f"Size comparison: size_a={size_a}, size_b={size_b}, "
            f"tolerance={_size['size_tolerance']:.1f} bp, similar_size={similar_size}"
        )

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
        log.debug(
            f"VERTICAL_MERGE|INV|REJECT_SIZE	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"size_a={size_a}	size_b={size_b}	regions_a={_regions_str_from_svcomposite(a)}	regions_b={_regions_str_from_svcomposite(b)}"
        )
        if verbose:
            print("Cannot merge inversions: sizes are not similar enough")
        return False
    if not similar_inversion_kmers:
        log.debug(
            f"VERTICAL_MERGE|INV|REJECT_KMER	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"kmer_sim={similarity:.3f}	threshold={min_kmer_overlap}	regions_a={_regions_str_from_svcomposite(a)}	regions_b={_regions_str_from_svcomposite(b)}"
        )
        if verbose:
            print("Cannot merge inversions: k-mer similarity is too low")
        return False

    log.debug(
        f"VERTICAL_MERGE|INV|ACCEPT	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
        f"size_a={size_a}	size_b={size_b}	kmer_sim={similarity:.3f}	regions_a={_regions_str_from_svcomposite(a)}	regions_b={_regions_str_from_svcomposite(b)}"
    )
    if verbose:
        print("Can merge inversions: all criteria passed")
    return True


def merge_inversions(
    inversions: list[SVcomposite],
    near: int,
    min_kmer_overlap: float,
    apriori_size_difference_fraction_tolerance: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    """Merge inversions that overlap on the reference and have similar sizes."""
    # 1) test if they have svPatterns
    if len(inversions) == 0:
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
    log.debug(f"MERGE_INVERSIONS|START	n_composites={len(inversions)}")
    for i in range(len(inversions)):
        for j in range(i + 1, len(inversions)):
            if can_merge_svComposites_inversions(
                a=inversions[i],
                b=inversions[j],
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
                verbose=verbose,
            ):
                uf.union(i, j)
    result: list[SVcomposite] = []
    for cc in uf.get_connected_components(allow_singletons=True):
        old_svcs = [inversions[idx] for idx in cc]
        new_svc = SVcomposite.from_SVpatterns(
            svPatterns=[
                svPattern for idx in cc for svPattern in inversions[idx].svPatterns
            ]
        )
        if len(cc) > 1:
            log.debug(
                f"TRANSFORMED::merge_inversions::vertical_merge:(to merged SVcomposite) "
                f"svComposites={[_svcomposite_log_id(svc) for svc in old_svcs]} "
                f"-->   {_svcomposite_log_id(new_svc)}"
            )
        else:
            log.debug(
                f"TRANSFORMED::merge_inversions::vertical merge: (no change): {_svcomposite_log_id(old_svcs[0])}"
            )
        result.append(new_svc)
    log.debug(
        f"MERGE_INVERSIONS|DONE	input={len(inversions)}	output={len(result)}"
    )
    return result


def can_merge_svComposites_breakends(
    a: SVcomposite,
    b: SVcomposite,
    apriori_size_difference_fraction_tolerance: float,
    near: int,
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
        log.debug(
            f"VERTICAL_MERGE|BND|REJECT_NOT_NEAR	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"regions_a={_regions_str_from_svcomposite(a)}	regions_b={_regions_str_from_svcomposite(b)}	tolerance={near}"
        )
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
        log.debug(
            f"VERTICAL_MERGE|BND|REJECT_NO_CONTEXT	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"contexts_a_len={len(contexts_a)}	contexts_b_len={len(contexts_b)}"
        )
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
        log.debug(
            f"VERTICAL_MERGE|BND|REJECT_KMER	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
            f"kmer_sim={similarity:.3f}	threshold={min_kmer_overlap}	regions_a={_regions_str_from_svcomposite(a)}	regions_b={_regions_str_from_svcomposite(b)}"
        )
        if verbose:
            print("Cannot merge breakends: k-mer similarity is too low")
        return False

    log.debug(
        f"VERTICAL_MERGE|BND|ACCEPT	a={_svcomposite_short_id(a)}	b={_svcomposite_short_id(b)}	"
        f"kmer_sim={similarity:.3f}	regions_a={_regions_str_from_svcomposite(a)}	regions_b={_regions_str_from_svcomposite(b)}"
    )
    if verbose:
        print("Can merge breakends: all criteria passed")
    return True


def merge_breakends(
    breakends: list[SVcomposite],
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
    log.debug(f"MERGE_BREAKENDS|START	n_composites={len(breakends)}")
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
        merged_crIDs: set[int] = set()
        for idx in cc:
            merged_crIDs.update(_crIDs_from_svcomposite(breakends[idx]))
        new_svc = SVcomposite.from_SVpatterns(
            svPatterns=[
                svPattern for idx in cc for svPattern in breakends[idx].svPatterns
            ]
        )
        if len(cc) > 1:
            old_svcs = [breakends[idx] for idx in cc]
            log.debug(
                f"TRANSFORMED::merge_breakends::vertical_merge:(to merged SVcomposite)  "
                f"svComposites={[_svcomposite_log_id(svc) for svc in old_svcs]} "
                f"-->   {_svcomposite_log_id(new_svc)}"
            )
        else:
            log.debug(
                f"TRANSFORMED::merge_breakends::vertical_merge: (no change): {_svcomposite_log_id(new_svc)}"
            )
        result.append(new_svc)
    log.debug(f"MERGE_BREAKENDS|DONE	input={len(breakends)}	output={len(result)}")
    return result


def _complexes_overlap(
    complex_a: SVcomposite, complex_b: SVcomposite, near: int
) -> bool:
    """Check if two complex SVcomposites overlap in all their breakends within a given tolerance.
    It is assumed that each sv complex is composed of exactly one sv pattern at this point."""
    if len(complex_a.svPatterns) != 1 or len(complex_b.svPatterns) != 1:
        raise ValueError(
            f"_complexes_overlap: SVcomposites must contain exactly one SVpattern each to check for overlap. Got {len(complex_a.svPatterns)} and {len(complex_b.svPatterns)}."
        )

    primitives_a = sorted(
        complex_a.svPatterns[0].SVprimitives, key=lambda p: (p.chr, p.ref_start)
    )
    primitives_b = sorted(
        complex_b.svPatterns[0].SVprimitives, key=lambda p: (p.chr, p.ref_start)
    )

    if len(primitives_a) != len(primitives_b):
        return False  # Different number of breakends, cannot overlap

    for prim_a, prim_b in zip(primitives_a, primitives_b, strict=True):
        if prim_a.chr != prim_b.chr:
            return False  # Different chromosomes, cannot overlap
        start_a, end_a = sorted((prim_a.ref_start, prim_a.ref_end))
        start_b, end_b = sorted((prim_b.ref_start, prim_b.ref_end))

        # Check if the intervals overlap within the tolerance
        if end_a + near < start_b or end_b + near < start_a:
            return False  # No overlap within tolerance

    return True  # All breakends overlap within tolerance


def _adjacencies_kmer_similarity(
    complex_a: SVcomposite,
    complex_b: SVcomposite,
    min_kmer_overlap: float,
) -> bool:
    """Calculate k-mer similarity between two complex SVcomposites based on their breakend sequences."""
    # each complex has exactly one sv pattern at this point
    if len(complex_a.svPatterns) != 1 or len(complex_b.svPatterns) != 1:
        raise ValueError(
            f"_adjacencies_kmer_similarity: SVcomposites must contain exactly one SVpattern each to calculate k-mer similarity. Got {len(complex_a.svPatterns)} and {len(complex_b.svPatterns)}."
        )
    if not issubclass(complex_a.sv_type, SVpatterns.SVpatternAdjacency):
        raise ValueError(
            f"_adjacencies_kmer_similarity: SVcomposite 'a' must be of type SVpatternAdjacency. Got {complex_a.sv_type}."
        )
    if not issubclass(complex_b.sv_type, SVpatterns.SVpatternAdjacency):
        raise ValueError(
            f"_adjacencies_kmer_similarity: SVcomposite 'b' must be of type SVpatternAdjacency. Got {complex_b.sv_type}."
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
        ],
        key=lambda x: (x[1].chr, x[1].ref_start),
    )
    data_b = sorted(
        [
            (idx, sv_prim, complex_b.svPatterns[0].get_sequence_context(svp_index=idx))
            for idx, sv_prim in enumerate(complex_b.svPatterns[0].SVprimitives)
        ],
        key=lambda x: (x[1].chr, x[1].ref_start),
    )
    # now calculate kmer similarity for each group
    group_a = [seq for _, _, seq in data_a]
    group_b = [seq for _, _, seq in data_b]
    similarity: float = kmer_similarity_of_groups(group_a=group_a, group_b=group_b)
    return similarity >= min_kmer_overlap


def merge_adjacencies(
    adjacencies: list[SVcomposite],
    near: int,
    min_kmer_overlap: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    """Merge complex SVcomposites that are spatially close in every underlying sv primitive (breakend), and if their kmer similarities are high enough."""
    # any combination of adjacencies can be merged here, so every combination of adjacencies must be tested.
    uf: UnionFind = UnionFind(range(len(adjacencies)))
    # if two adjacencies match, connect them. Every connected component will be merged into one complex at the end.
    for i in range(len(adjacencies)):
        for j in range(i + 1, len(adjacencies)):
            id_a = _svcomposite_short_id(adjacencies[i])
            id_b = _svcomposite_short_id(adjacencies[j])
            if _complexes_overlap(
                complex_a=adjacencies[i], complex_b=adjacencies[j], near=near
            ):
                if _adjacencies_kmer_similarity(
                    complex_a=adjacencies[i],
                    complex_b=adjacencies[j],
                    min_kmer_overlap=min_kmer_overlap,
                ):
                    uf.union(i, j)
                    log.debug("VERTICAL_MERGE|ADJ|ACCEPT|a=%s|b=%s", id_a, id_b)
                else:
                    log.debug("VERTICAL_MERGE|ADJ|REJECT_KMER|a=%s|b=%s", id_a, id_b)
            else:
                log.debug("VERTICAL_MERGE|ADJ|REJECT_NO_OVERLAP|a=%s|b=%s", id_a, id_b)
    result: list[SVcomposite] = []
    for cc in uf.get_connected_components(allow_singletons=True):
        merged_crIDs: set[int] = set()
        merged_consIDs: list[str] = []
        for idx in cc:
            merged_crIDs.update(_crIDs_from_svcomposite(adjacencies[idx]))
            merged_consIDs.extend(p.consensusID for p in adjacencies[idx].svPatterns)
        log.debug(
            "VERTICAL_MERGE|ADJ|MERGED_GROUP|component_size=%d|crIDs=%s|consensusIDs=%s",
            len(cc),
            sorted(merged_crIDs),
            merged_consIDs,
        )
        res: SVcomposite = SVcomposite.from_SVpatterns(
            svPatterns=[
                svPattern for idx in cc for svPattern in adjacencies[idx].svPatterns
            ]
        )
        # add logging that sv patterns were merged into one sv composite, with the crIDs and consensusIDs of the merged sv patterns
        log.debug(
            f"TRANSFORMED::merge_adjacencies::vertical_merge:(to merged SVcomposite) "
            f"svComposites={[adjacencies[idx]._log_id() for idx in cc]}. --> {res._log_id()}"
        )
        result.append(res)

    log.debug("VERTICAL_MERGE|ADJ|DONE|n_merged=%d", len(result))
    return result


def merge_svComposites_across_chromosomes(
    svComposites: list[SVcomposite],
    apriori_size_difference_fraction_tolerance: float,
    near: int,
    min_kmer_overlap: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    """
    Merge SVcomposites across chromosomes without parallelization.
    This function processes all chromosomes sequentially.
    """
    log.debug("CHROMOSOME_MERGE|ACROSS_CHR|START|n_composites=%d", len(svComposites))

    # build overlap trees for all chromosomes
    overlap_trees: dict[str, IntervalTree] = defaultdict(IntervalTree)
    for idx, svComposite in enumerate(svComposites):
        for svprimitive in (
            svComposite.svPatterns[0].SVprimitives
        ):  # assuming there is only one sv pattern per sv composite at this point
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
    log.debug(
        "CHROMOSOME_MERGE|ACROSS_CHR|OVERLAP_TREES_BUILT|n_chromosomes=%d",
        len(overlap_trees),
    )
    uf = UnionFind(range(len(svComposites)))

    for _chr_name, tree in overlap_trees.items():
        tree.merge_overlaps(data_reducer=lambda x, y: x.union(y))

    for _chr_name, tree in overlap_trees.items():
        for interval in tree:
            indices = sorted(interval.data)
            if len(indices) > 1:
                merged_crIDs: set[int] = set()
                for idx in indices:
                    merged_crIDs.update(_crIDs_from_svcomposite(svComposites[idx]))
                log.debug(
                    "CHROMOSOME_MERGE|ACROSS_CHR|OVERLAP_UNION|chr=%s|n_indices=%d|crIDs=%s",
                    _chr_name,
                    len(indices),
                    sorted(merged_crIDs),
                )
                for i in range(len(indices)):
                    for j in range(i + 1, len(indices)):
                        uf.union_by_name(indices[i], indices[j])

    # vertical merging
    connected_components = uf.get_connected_components(allow_singletons=True)
    merged_svComposites = []
    log.debug(
        "CHROMOSOME_MERGE|ACROSS_CHR|VERTICAL_MERGE_START|n_groups=%d",
        len(connected_components),
    )
    for cc in tqdm(connected_components, desc="All Chromosomes"):
        svComposite_group = [svComposites[idx] for idx in cc]
        group_crIDs: set[int] = set()
        for svc in svComposite_group:
            group_crIDs.update(_crIDs_from_svcomposite(svc))

        if len(svComposite_group) > 1:
            log.debug(
                "CHROMOSOME_MERGE|ACROSS_CHR|MERGE_GROUP|size=%d|crIDs=%s|consensusIDs=%s",
                len(svComposite_group),
                sorted(group_crIDs),
                [svc.svPatterns[0].consensusID for svc in svComposite_group],
            )

            merged = vertically_merged_svComposites_from_group(
                group=svComposite_group,
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
                verbose=verbose,
            )

            merged_svComposites.extend(merged)
        else:
            log.debug(
                "CHROMOSOME_MERGE|ACROSS_CHR|SINGLETON|crIDs=%s",
                sorted(group_crIDs),
            )
            merged_svComposites.append(svComposite_group[0])
    log.debug(
        "CHROMOSOME_MERGE|ACROSS_CHR|DONE|in=%d|out=%d",
        len(svComposites),
        len(merged_svComposites),
    )
    return merged_svComposites


def merge_svComposites_for_chromosome(
    chr_name: str,
    svComposites: list[SVcomposite],
    apriori_size_difference_fraction_tolerance: float,
    near: int,
    min_kmer_overlap: float,
    verbose: bool = False,
) -> list[SVcomposite]:
    """
    Merge SVcomposites for a single chromosome.
    This function is designed to be run in parallel for different chromosomes.
    """
    log.debug(
        "CHROMOSOME_MERGE|CHR|START|chr=%s|n_composites=%d", chr_name, len(svComposites)
    )
    # This is the entry point of a worker process when threads > 1, so the
    # per-process rate limiter and counters belong here.

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
    log.debug(
        "CHROMOSOME_MERGE|CHR|OVERLAP_TREE_BUILT|chr=%s|n_intervals=%d",
        chr_name,
        len(overlaptree),
    )
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
                        log.debug(
                            "CHROMOSOME_MERGE|CHR|SKIP_SAME_SAMPLE_CONSENSUS|chr=%s|consensusID=%s|sample=%s",
                            chr_name,
                            svComposites[indices[i]].svPatterns[0].consensusID,
                            svComposites[indices[i]].svPatterns[0].samplename,
                        )
                        continue
                    else:
                        uf.union_by_name(indices[i], indices[j])

    # Get connected components and merge groups
    connected_components = uf.get_connected_components(allow_singletons=True)
    merged_svComposites = []

    log.debug(
        "CHROMOSOME_MERGE|CHR|VERTICAL_MERGE_START|chr=%s|n_groups=%d",
        chr_name,
        len(connected_components),
    )
    for cc in tqdm(connected_components, desc=f"Chr {chr_name}"):
        svComposite_group = [svComposites[idx] for idx in cc]
        group_crIDs: set[int] = set()
        for svc in svComposite_group:
            group_crIDs.update(_crIDs_from_svcomposite(svc))

        if len(svComposite_group) > 1:
            log.debug(
                "CHROMOSOME_MERGE|CHR|MERGE_GROUP|chr=%s|size=%d|crIDs=%s|consensusIDs=%s",
                chr_name,
                len(svComposite_group),
                sorted(group_crIDs),
                [svc.svPatterns[0].consensusID for svc in svComposite_group],
            )
            merged = vertically_merged_svComposites_from_group(
                group=svComposite_group,
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
                verbose=verbose,
            )
            merged_svComposites.extend(merged)
        else:
            log.debug(
                "CHROMOSOME_MERGE|CHR|SINGLETON|chr=%s|crIDs=%s",
                chr_name,
                sorted(group_crIDs),
            )
            merged_svComposites.append(svComposite_group[0])

    log.debug(
        "CHROMOSOME_MERGE|CHR|DONE|chr=%s|in=%d|out=%d",
        chr_name,
        len(svComposites),
        len(merged_svComposites),
    )
    return merged_svComposites


def merge_svComposites(
    svComposites: list[SVcomposite],
    apriori_size_difference_fraction_tolerance: float,
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
        if issubclass(svComposite.sv_type, SVpatterns.SVpatternAdjacency):
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
                near=near,
                min_kmer_overlap=min_kmer_overlap,
                verbose=verbose,
            )
            merged_svComposites.extend(chr_merged)

    # cpx_groups can jump across chromosomes, so they cannot be parallelized by chromosome.
    # they are merged here sequentially.
    merged_svComposites.extend(
        merge_svComposites_across_chromosomes(
            svComposites=cpx_groups,
            apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
            near=near,
            min_kmer_overlap=min_kmer_overlap,
            verbose=verbose,
        )
    )

    log.info(
        f"Merged {len(svComposites)} SVcomposites into {len(merged_svComposites)} across all chromosomes."
    )
    return merged_svComposites
