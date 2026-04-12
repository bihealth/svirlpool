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
from .SVcomposite import SVcomposite
from .svcomposite_merging import merge_svComposites
from .svcomposite_utils import _svcomposite_log_id

log = logging.getLogger(__name__)


# ---- Logging helper functions (imported from svcomposite_utils) ---- #


def _crIDs_from_svpattern(svp: SVpatterns.SVpatternType) -> int:
    """Extract crID from an SVpattern. The crID is the integer prefix of the consensusID (format: crID.subID)."""
    try:
        return int(svp.consensusID.split(".")[0])
    except (ValueError, IndexError):
        return -1


# ---- End logging helpers ---- #


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

# Map SV type strings to lists of pattern types (handles multiple types for "BND")
SUPPORTED_SV_TYPE_STRINGS_INVERSE: dict[str, list[type[SVpatterns.SVpatternType]]] = {}
for pattern_type in SUPPORTED_SV_TYPES:
    sv_type_str = pattern_type.get_sv_type()
    if sv_type_str not in SUPPORTED_SV_TYPE_STRINGS_INVERSE:
        SUPPORTED_SV_TYPE_STRINGS_INVERSE[sv_type_str] = []
    SUPPORTED_SV_TYPE_STRINGS_INVERSE[sv_type_str].append(pattern_type)

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
                    f"DROPPED::parse_candidate_regions_file::MALFORMED LINE: {line}"
                )
                continue
            samplename = parts[0]
            crIDs = set(map(int, parts[1].split(",")))
            result[samplename] = crIDs
    log.info(f"Parsed candidate regions for {len(result)} samples")
    return result


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
    collapse_repeats: bool = True,
) -> list[SVcomposite]:

    result: list[SVcomposite] = []
    if len(svPatterns) == 0:
        return result
    # split the insertiosn and deletions. Add the rest to results immediately.
    # filter svPatterns for supported types
    _svPatterns = [svp for svp in svPatterns if type(svp) in sv_types]
    unsupported = [svp for svp in svPatterns if type(svp) not in sv_types]
    if unsupported:
        for svp in unsupported:
            log.debug(
                f"DROPPED::svPatterns_to_horizontally_merged_svComposites::UNSUPPORTED SV TYPE: {svp._log_id()}"
            )
        log.warning(
            f"DROPPED::svPatterns_to_horizontally_merged_svComposites::UNSUPPORTED SV TYPES: {len(unsupported)} SVpatterns of types {[type(svp).__name__ for svp in unsupported]}"
        )
        unsupported.clear()
    indels = [
        svp
        for svp in _svPatterns
        if type(svp) in (SVpatterns.SVpatternInsertion, SVpatterns.SVpatternDeletion)
    ]
    others = [
        svp
        for svp in _svPatterns
        if type(svp)
        not in (SVpatterns.SVpatternInsertion, SVpatterns.SVpatternDeletion)
    ]
    _svPatterns.clear()

    # convert other svPatterns directly to svComposites
    for svp in others:
        svc = SVcomposite.from_SVpattern(svp)
        log.debug(
            f"TRANSFORMED::svPatterns_to_horizontally_merged_svComposites::from_SVpattern: {svp._log_id()} -> {_svcomposite_log_id(svc)}"
        )
        result.append(svc)

    # sort svPatterns by consensusID and read_start
    indels = sorted(indels, key=lambda x: (x.consensusID, x.read_start))

    groups = groupby(indels, key=lambda x: x.consensusID)

    # loop svPatterns of each group and connect them if they share at least one repeatID
    for _consensusID, group in groups:
        group = list(group)
        crID = int(_consensusID.split(".")[0]) if "." in _consensusID else -1

        if len(group) == 1:
            svc = SVcomposite.from_SVpatterns(group)
            log.debug(
                f"TRANSFORMED::svPatterns_to_horizontally_merged_svComposites::singleton_to_composite: {group[0]._log_id()} -> {_svcomposite_log_id(svc)}"
            )
            result.append(svc)
            continue

        # Create a union-find structure for this group
        uf_group = UnionFind(range(len(group)))

        # this is the point where horizontal merge can be prevented by skipping the union step
        if collapse_repeats:
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    if group[i].repeatIDs.intersection(group[j].repeatIDs):
                        # TODO: Edge case, where indels in duplicated overlapping aligned fragments are concatenated horizontally
                        log.debug(
                            f"HORIZONTAL_MERGE|UNION_BY_REPEATID	crID={crID}	consensusID={_consensusID}	"
                            f"pattern_i={group[i]._log_id()}    "
                            f"pattern_j={group[j]._log_id()}    "
                            f"shared_repeatIDs={group[i].repeatIDs.intersection(group[j].repeatIDs)}"
                        )
                        uf_group.union(i, j)

        # Collect connected components
        connected_components = uf_group.get_connected_components(allow_singletons=True)

        for component in connected_components:
            sv_patterns_in_component = [group[idx] for idx in component]
            if len(sv_patterns_in_component) > 0:
                pattern_ids = [svp._log_id() for svp in sv_patterns_in_component]
                svc = SVcomposite.from_SVpatterns(sv_patterns_in_component)
                log.debug(
                    f"TRANSFORMED::svPatterns_to_horizontally_merged_svComposites::merge_by_repeatID: [{'; '.join(pattern_ids)}] -> {_svcomposite_log_id(svc)}"
                )
                result.append(svc)

    return result


# horizontal merge
def generate_svComposites_from_dbs(
    input: list[str | Path],
    sv_types: set[type[SVpatterns.SVpatternType]],
    candidate_regions_filter: dict[str, set[int]] | None = None,
    collapse_repeats: bool = True,
) -> list[SVcomposite]:
    log.debug("HORIZONTAL_MERGE|LOAD_DBS|n_dbs=%d", len(input))
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
                "HORIZONTAL_MERGE|CR_FILTER|sample=%s|n_crIDs=%d",
                samplename,
                len(crIDs_to_query),
            )

        # Read SVpatterns from database, optionally filtered by crIDs
        svPatterns = SVpatterns.read_svPatterns_from_db(
            database=Path(p), crIDs=crIDs_to_query
        )

        # count the type of each svPattern
        sv_type_counts: dict[str, int] = {}
        for svp in svPatterns:
            sv_type_counts[svp.get_sv_type()] = (
                sv_type_counts.get(svp.get_sv_type(), 0) + 1
            )
        log.info(
            "HORIZONTAL_MERGE|LOADED|sample=%s|n_patterns=%d|type_counts=%s",
            samplename,
            len(svPatterns),
            sv_type_counts,
        )

        svComposites.extend(
            svPatterns_to_horizontally_merged_svComposites(
                svPatterns, sv_types=sv_types, collapse_repeats=collapse_repeats
            )
        )

        # after horizontal merging, count again the types of sv composites
        sv_composite_type_counts: dict[str, int] = {}
        for svc in svComposites:
            sv_composite_type_counts[svc.sv_type.get_sv_type()] = (
                sv_composite_type_counts.get(svc.sv_type.get_sv_type(), 0) + 1
            )
        log.debug(
            "HORIZONTAL_MERGE|AFTER_SAMPLE|sample=%s|n_composites=%d|type_counts=%s",
            samplename,
            len(svComposites),
            sv_composite_type_counts,
        )

    log.debug("HORIZONTAL_MERGE|ALL_DONE|total_composites=%d", len(svComposites))
    return svComposites


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

    def to_log_id(self) -> str:
        return f"{self.svtype}|{self.chrname}:{self.start}-{self.end}|consensusIDs={','.join(self.consensusIDs)}|svlen={self.svlen}|mateid={self.mateid}"

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

        # Special handling for BND (breakend) records
        if self.svtype == "BND":
            # For BND records, alt_seq contains the pre-formatted ALT field
            # Replace the placeholder "N" with the actual reference base
            ref_str = refbase
            alt_str = (
                alt_seq.replace("N", refbase) if isinstance(alt_seq, str) else refbase
            )
        else:
            # Normal handling for other SV types
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
    min_radius: int = 1,
) -> set[int]:
    if end < start:
        raise ValueError(
            f"get_ref_reads_from_covtrees: end {end} is less than start {start} for sample {samplename} at {chrname}:{start}-{end}."
        )
    min_radius = abs(min_radius)
    if end - start < min_radius * 2:
        difference = min_radius * 2 - (end - start)
        start -= floor(difference * 0.5)
        end += ceil(difference * 0.5)
        if start < 0:
            start = 0
    all_reads: set[int] = {
        int(it.data) for it in covtrees[samplename][chrname][start:end]
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
    min_radius: int = 1,
    breakpoint_mode: bool = False,
    breakpoint_margin: int = 100,
) -> Genotype:
    genotype: Genotype
    if chrname not in covtrees.get(samplename, {}):
        log.warning(
            f"SVcalls_from_SVcomposite: Chromosome {chrname} not found in coverage tree for sample {samplename}. Assigning 0/0 genotype."
        )
        genotype = create_wild_type_genotype(samplename=samplename, total_coverage=0)
        return genotype
    if breakpoint_mode:
        # For deletions: query at both breakpoints rather than across the full
        # deletion span.  Alt-supporting reads have effective_intervals ending
        # at the start breakpoint and beginning at the end breakpoint, so an
        # interior-only query misses them for large deletions.
        reads_at_start: set[int] = get_ref_reads_from_covtrees(
            samplename=samplename,
            chrname=chrname,
            start=start,
            end=start,
            covtrees=covtrees,
            min_radius=breakpoint_margin,
        )
        reads_at_end: set[int] = get_ref_reads_from_covtrees(
            samplename=samplename,
            chrname=chrname,
            start=end,
            end=end,
            covtrees=covtrees,
            min_radius=breakpoint_margin,
        )
        all_reads: set[int] = reads_at_start | reads_at_end
    else:
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
        genotype = create_wild_type_genotype(samplename=samplename, total_coverage=0)
        return genotype

    # Query copy number for this locus from CN tracks
    copy_number: int = 2  # Default diploid - if no copy number track is available
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
    log.debug(
        f"GENOTYPE|RESULT    sample={samplename}   region={chrname}:{start}-{end}    GT={gt}    GQ={genotype.genotype_quality}    ref={len(ref_reads)}    alt={len(alt_reads)}    total={len(all_reads)}    CN={copy_number}",
    )
    return genotype


def SVcalls_from_SVcomposite(
    svComposite: SVcomposite,
    covtrees: dict[str, dict[str, IntervalTree]],
    cn_tracks: dict[
        str, dict[str, IntervalTree]
    ],  # saplename -> chrname -> IntervalTree with copy number
    find_leftmost_reference_position: bool,
    symbolic_threshold: int,
) -> list[SVcall]:
    # intra-alignment fragment variants (closed locus), e.g. INS, DEL, INV, DUP
    #   have one chr, start, end on the reference
    # inter-alignment fragment variants (interrupted locus), e.g translocations, break points
    #   have multiple chr, start, end on the reference
    composite_id = _svcomposite_log_id(svComposite)

    if svComposite.sv_type not in SUPPORTED_SV_TYPES:
        log.debug(
            f"DROPPED::SVcalls_from_SVcomposite::SVTYPE NOT SUPPORTED: {composite_id}",
        )
        return []
    all_alt_reads: dict[str, set[int]] = (
        svComposite.get_alt_readnamehashes_per_sample()
    )  # {samplename: {readname, ...}}

    if (
        issubclass(svComposite.sv_type, SVpatterns.SVpatternDeletion)
        or issubclass(svComposite.sv_type, SVpatterns.SVpatternInsertion)
        or issubclass(svComposite.sv_type, SVpatterns.SVpatternInversion)
        or issubclass(svComposite.sv_type, SVpatterns.SVpatternSingleBreakend)
    ):
        res = svcall_object_from_svcomposite(
            svComposite=svComposite,
            covtrees=covtrees,
            cn_tracks=cn_tracks,
            find_leftmost_reference_position=find_leftmost_reference_position,
            all_alt_reads=all_alt_reads,
        )
        log.debug(
            f"TRANSFORMED::SVcalls_from_SVcomposite::svcall_object_from_svcomposite:(to DEL, INS, INV, BND) {composite_id}; TRANSFORMED TO: {res.to_log_id()}",
        )
        return [res]
    elif issubclass(svComposite.sv_type, SVpatterns.SVpatternAdjacency):
        log.debug(
            f"TRANSFORMED::SVcalls_from_SVcomposite::svcall_object_from_svcomposite:(to ADJACENCY) {composite_id}",
        )
        res = svcall_objects_from_Adjacencies(
            svComposite=svComposite,
            covtrees=covtrees,
            cn_tracks=cn_tracks,
            all_alt_reads=all_alt_reads,
            symbolic_threshold=symbolic_threshold,
        )
        # log each res
        for _res in res:
            log.debug(
                f"TRANSFORMED::SVcalls_from_SVcomposite::svcall_objects_from_Adjacencies:(to ADJACENCY) {composite_id}; TRANSFORMED TO: {_res.to_log_id()}",
            )
        return res
    else:
        log.warning(
            f"DROPPED::SVcalls_from_SVcomposite: (SV SUBCLASS NOT YET IMPLEMENTED)  {composite_id}."
        )
        return []


def svcall_object_from_svcomposite(
    svComposite: SVcomposite,
    covtrees: dict[str, dict[str, IntervalTree]],
    cn_tracks: dict[str, dict[str, IntervalTree]],
    find_leftmost_reference_position: bool,
    all_alt_reads: dict[str, set[int]],
) -> SVcall:
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
    is_deletion = issubclass(svComposite.sv_type, SVpatterns.SVpatternDeletion)
    genotypes: dict[str, Genotype] = {
        samplename: genotype_of_sample(
            samplename=samplename,
            chrname=chrname,
            start=start,
            end=end,
            raw_alt_reads=all_alt_reads[samplename],
            covtrees=covtrees,
            cn_tracks=cn_tracks,
            breakpoint_mode=is_deletion,
        )
        for samplename in all_alt_reads.keys()
    }
    # start and end need to be adjusted based on sv_type and need to be determined for the whole sv composite

    if issubclass(svComposite.sv_type, SVpatterns.SVpatternDeletion):
        end = start + svlen
    elif issubclass(svComposite.sv_type, SVpatterns.SVpatternInsertion):
        end = start + 1

    alt_seq: bytes | None = (
        pickle.dumps(svComposite.get_alt_sequence())
        if svComposite.get_alt_sequence()
        else None
    )

    ref_seq: bytes | None = (
        pickle.dumps(svComposite.get_ref_sequence())
        if svComposite.get_ref_sequence()
        else None
    )

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

    gt_summary = {s: g.genotype for s, g in genotypes.items()}
    log.debug(
        f"TRANSFORMED::svcall_object_from_svcomposite::(generated SVCALL): {_svcomposite_log_id(svComposite)}; GENOTYPES: genotypes={gt_summary}",
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


def format_bnd_alt_field(
    ref_base: str,
    mate_chr: str,
    mate_pos: int,
    sv_type: int,
    aln_is_reverse: bool,
    inserted_sequence: str = "",
    use_symbolic: bool = False,
    sequence_id: str | None = None,
    insertion_start: int = 1,
    insertion_end: int | None = None,
) -> str:
    """
    Format the ALT field for a BND (breakend) record according to VCF specification.

    Args:
        ref_base: The reference base at this breakend position (REF field, typically 'N' or actual base)
        mate_chr: Chromosome of the mate breakend
        mate_pos: Position of the mate breakend (1-based for VCF output)
        sv_type: SV type of this breakend (3=BNDL/left, 4=BNDR/right)
        aln_is_reverse: Whether the alignment is on the reverse strand
        inserted_sequence: Optional inserted sequence between the breakends (ignored if use_symbolic=True)
        use_symbolic: If True, use symbolic notation for large insertions with sequence_id
        sequence_id: ID of the inserted sequence (e.g., "ctg1") when using symbolic notation
        insertion_start: Start position in the symbolic sequence (1-based, default=1)
        insertion_end: End position in the symbolic sequence (1-based); if None, uses full length

    Returns:
        Formatted ALT field string according to VCF BND specification

    The four cases are:
    - t[p[  : piece extending to the right of p is joined after t
    - t]p]  : reverse comp piece extending left of p is joined after t
    - ]p]t  : piece extending to the left of p is joined before t
    - [p[t  : reverse comp piece extending right of p is joined before t

    Where:
    - t = ref_base + inserted_sequence (replacement string) or ref_base for symbolic
    - p = mate_chr:mate_pos (mate position) or <sequence_id>:pos for symbolic insertions
    - sv_type: 3 = BNDL (left break), 4 = BNDR (right break)

    For large insertions (use_symbolic=True), uses symbolic notation like:
    - C[<ctg1>:1[  (reference to start of contig)
    - ]<ctg1>:329]A  (reference to end of contig)
    """
    # Prepare the replacement string t
    if use_symbolic:
        # For symbolic notation, t is just the ref base (insertion is referenced separately)
        t = ref_base
        # For symbolic insertions, p points to the contig
        if sequence_id is None:
            raise ValueError("sequence_id must be provided when use_symbolic=True")
        # Determine which end of the insertion to reference based on the breakend
        # First breakend points to start of insertion, second points to end
        p = f"<{sequence_id}>:{insertion_start if sv_type == 4 else (insertion_end or len(inserted_sequence))}"
    else:
        # For explicit sequence, t includes ref base + insertion
        t = ref_base + inserted_sequence
        # Regular mate position
        p = f"{mate_chr}:{mate_pos}"

    # Determine the ALT format based on sv_type and orientation
    # sv_type 3 = BNDL (left breakend), sv_type 4 = BNDR (right breakend)

    if sv_type == 4:  # BNDR (right breakend)
        if not aln_is_reverse:
            # Right break, forward strand: t[p[
            # Piece extending to the right of p is joined after t
            return f"{t}[{p}["
        else:
            # Right break, reverse strand: t]p]
            # Reverse comp piece extending left of p is joined after t
            return f"{t}]{p}]"
    elif sv_type == 3:  # BNDL (left breakend)
        if not aln_is_reverse:
            # Left break, forward strand: ]p]t
            # Piece extending to the left of p is joined before t
            return f"]{p}]{t}"
        else:
            # Left break, reverse strand: [p[t
            # Reverse comp piece extending right of p is joined before t
            return f"[{p}[{t}"
    else:
        raise ValueError(
            f"Invalid sv_type for BND: {sv_type}. Expected 3 (BNDL) or 4 (BNDR)"
        )


def svcall_objects_from_Adjacencies(
    svComposite: SVcomposite,
    covtrees: dict[str, dict[str, IntervalTree]],
    cn_tracks: dict[str, dict[str, IntervalTree]],
    all_alt_reads: dict[str, set[int]],
    symbolic_threshold: int,
) -> list[SVcall]:
    """Generate SVcall objects from a SVcomposite that represents novel adjacencies with two connected break ends of each sample.
    The given svComposite generates two SVcall objects, that both represent one end of the novel adjacency.
    """
    if not issubclass(svComposite.sv_type, SVpatterns.SVpatternAdjacency):
        raise ValueError(
            f"svcall_objects_from_Adjacencies: svComposite of type {svComposite.sv_type.__name__} is not of type SVpatternAdjacency."
        )

    # Collect consensusIDs from all SVpatterns in the composite
    consensusIDs: list[str] = list({
        svPattern.samplenamed_consensusID for svPattern in svComposite.svPatterns
    })

    # Get the representative SVpattern (with highest support)
    representative_pattern = svComposite.get_representative_SVpattern()

    # Verify it's an Adjacency pattern (should always be true due to earlier check)
    if not isinstance(representative_pattern, SVpatterns.SVpatternAdjacency):
        raise ValueError(
            f"Representative pattern is not SVpatternAdjacency, got {type(representative_pattern).__name__}"
        )

    # Adjacency patterns must have exactly 2 SVprimitives (the two break ends)
    if len(representative_pattern.SVprimitives) != 2:
        raise ValueError(
            f"SVpatternAdjacency must have exactly 2 SVprimitives, got {len(representative_pattern.SVprimitives)}"
        )

    # Extract the two break ends' reference positions
    svprimitive_0 = representative_pattern.SVprimitives[0]
    svprimitive_1 = representative_pattern.SVprimitives[1]

    # Break end 0 location
    chr0 = svprimitive_0.chr
    start0 = svprimitive_0.ref_start
    end0 = start0 + 1  # Break ends are essentially point locations, so end = start + 1

    # Break end 1 location
    chr1 = svprimitive_1.chr
    start1 = svprimitive_1.ref_start
    end1 = start1 + 1

    # Get the inserted sequence between the two break ends
    inserted_sequence = representative_pattern.get_sequence()
    if inserted_sequence is None:
        inserted_sequence = ""

    # svlen is the length of the inserted sequence (or 1 if none)
    svlen: int = len(inserted_sequence) if inserted_sequence else 1

    # ref_seq is not important for adjacencies
    ref_seq: bytes | None = None

    # Determine if we should use symbolic notation for large insertions
    # Using the symbolic_threshold parameter passed from command-line arguments
    use_symbolic = len(inserted_sequence) > symbolic_threshold

    # Generate a sequence ID for symbolic notation if needed
    sequence_id: str | None = None
    if use_symbolic:
        # Create a unique sequence ID based on the SVcomposite
        sequence_id = f"ctg_{representative_pattern.samplename}_{representative_pattern.consensusID}_{svprimitive_0.svID}_{svprimitive_1.svID}"

    # Format the ALT field for each breakend according to VCF BND specification
    # Note: We use "N" as placeholder for ref_base since actual base will be retrieved in to_vcf_line
    # Breakend 0 points to breakend 1 (or to start of insertion contig)
    alt_str_0 = format_bnd_alt_field(
        ref_base="N",
        mate_chr=chr1,
        mate_pos=start1 + 1,  # VCF is 1-based, convert from 0-based
        sv_type=svprimitive_0.sv_type,
        aln_is_reverse=svprimitive_0.aln_is_reverse,
        inserted_sequence=inserted_sequence,
        use_symbolic=use_symbolic,
        sequence_id=sequence_id,
        insertion_start=1,
        insertion_end=len(inserted_sequence) if use_symbolic else None,
    )

    # Breakend 1 points to breakend 0 (or to end of insertion contig)
    alt_str_1 = format_bnd_alt_field(
        ref_base="N",
        mate_chr=chr0,
        mate_pos=start0 + 1,  # VCF is 1-based, convert from 0-based
        sv_type=svprimitive_1.sv_type,
        aln_is_reverse=svprimitive_1.aln_is_reverse,
        inserted_sequence=inserted_sequence,
        use_symbolic=use_symbolic,
        sequence_id=sequence_id,
        insertion_start=1,
        insertion_end=len(inserted_sequence) if use_symbolic else None,
    )

    # Store the formatted ALT strings as pickled bytes
    alt_seq_0: bytes = pickle.dumps(alt_str_0)
    alt_seq_1: bytes = pickle.dumps(alt_str_1)

    # Generate genotypes for all samples
    # For break end 0
    genotypes_0: dict[str, Genotype] = {
        samplename: genotype_of_sample(
            samplename=samplename,
            chrname=chr0,
            start=start0,
            end=end0,
            raw_alt_reads=all_alt_reads[samplename],
            covtrees=covtrees,
            cn_tracks=cn_tracks,
        )
        for samplename in all_alt_reads.keys()
    }

    # For break end 1
    genotypes_1: dict[str, Genotype] = {
        samplename: genotype_of_sample(
            samplename=samplename,
            chrname=chr1,
            start=start1,
            end=end1,
            raw_alt_reads=all_alt_reads[samplename],
            covtrees=covtrees,
            cn_tracks=cn_tracks,
        )
        for samplename in all_alt_reads.keys()
    }

    # Quality filters (same logic as svcall_object_from_svcomposite)
    # We use the maximum of both break ends for pass_altreads
    max_var_reads_0 = max(genotypes_0.items(), key=lambda x: x[1].var_reads)[
        1
    ].var_reads
    max_var_reads_1 = max(genotypes_1.items(), key=lambda x: x[1].var_reads)[
        1
    ].var_reads
    pass_altreads: bool = max(max_var_reads_0, max_var_reads_1) >= 3

    pass_gq = True
    passing: bool = pass_altreads and pass_gq

    # Check if call is precise (not in a repeat)
    precise: bool = not bool(
        sum(
            len(svPrimitive.repeatIDs)
            for svPattern in svComposite.svPatterns
            for svPrimitive in svPattern.SVprimitives
        )
    )

    # Generate unique mateid base from representative pattern
    # Format: samplename.consensusID.svID_0.svID_1
    mateid_base = f"{representative_pattern.samplename}.{representative_pattern.consensusID}.{svprimitive_0.svID}.{svprimitive_1.svID}"

    # Create mateid for each break end (they reference each other)
    mateid_0 = f"{mateid_base}.0"
    mateid_1 = f"{mateid_base}.1"

    # Create SVcall for break end 0
    svcall_0 = SVcall(
        genotypes=genotypes_0,
        passing=passing,
        chrname=chr0,
        start=start0,
        end=end0,
        svtype="BND",
        svlen=svlen,
        pass_altreads=pass_altreads,
        pass_gq=pass_gq,
        precise=precise,
        mateid=mateid_1,  # Points to the mate (break end 1)
        consensusIDs=consensusIDs,
        ref_sequence=ref_seq,
        alt_sequence=alt_seq_0,  # BND-formatted ALT field for breakend 0
        sequence_id=sequence_id,  # Set if using symbolic notation for large insertions
    )

    # Create SVcall for break end 1
    svcall_1 = SVcall(
        genotypes=genotypes_1,
        passing=passing,
        chrname=chr1,
        start=start1,
        end=end1,
        svtype="BND",
        svlen=svlen,
        pass_altreads=pass_altreads,
        pass_gq=pass_gq,
        precise=precise,
        mateid=mateid_0,  # Points to the mate (break end 0)
        consensusIDs=consensusIDs,
        ref_sequence=ref_seq,
        alt_sequence=alt_seq_1,  # BND-formatted ALT field for breakend 1
        sequence_id=sequence_id,  # Set if using symbolic notation for large insertions
    )

    composite_id = _svcomposite_log_id(svComposite)
    log.debug(
        "TRANSFORMED::svcall_objects_from_Adjacencies::(BND_PAIR_CREATED)|pass=%s|pass_altreads=%s|precise=%s|"
        "bnd0=%s:%d|bnd1=%s:%d|svlen=%d|%s",
        passing,
        pass_altreads,
        precise,
        chr0,
        start0,
        chr1,
        start1,
        svlen,
        composite_id,
    )

    return [svcall_0, svcall_1]


def get_svComposite_interval_on_reference(
    svComposite: SVcomposite, find_leftmost_reference_position: bool
) -> tuple[str, int, int]:
    if svComposite.sv_type not in SUPPORTED_SV_TYPES:
        raise ValueError(
            f"get_svComposite_indel_interval called with svComposite that is neither {', '.join(t.__name__ for t in SUPPORTED_SV_TYPES)}: {svComposite}"
        )
    # measure time. If the execution of this function takes longer than 2 seconds, print the input svComposite so it can be debugged.
    time_start = datetime.now()
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
            weight = len(svPattern.get_supporting_reads()) * svPattern.get_size()
        weighted_regions.append((region, weight))

    time_end = datetime.now()
    time_diff = (time_end - time_start).total_seconds()
    if time_diff > 2:
        log.warning(
            f"get_svComposite_interval_on_reference took {time_diff} seconds for svComposite: {svComposite}"
        )

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
        '##INFO=<ID=CONSENSUSIDs,Number=.,Type=String,Description="ID of the consensus that this SV originates from. Other consensus sequences can also be involved.">'
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


def _parse_consensusID_parts(
    consensus_id: str,
) -> tuple[str, int, int] | None:
    """Parse a samplenamed consensusID like 'HG002:731.0' into (samplename, crID, subID).
    Returns None if the format is unexpected."""
    try:
        samplename, cr_part = consensus_id.split(":", 1)
        cr_str, sub_str = cr_part.split(".", 1)
        return samplename, int(cr_str), int(sub_str)
    except (ValueError, IndexError):
        return None


def correct_genotypes_for_multi_assembly_loci(
    svCalls: list[SVcall],
) -> list[SVcall]:
    """Correct genotypes where multiple consensus assemblies from the same candidate
    region prove that a locus is heterozygous.

    When a candidate region produces two (or more) distinct consensus assemblies for
    a sample, each assembly represents a different allele.  Variants originating from
    different assemblies of the same crID must therefore be heterozygous (0/1) in a
    diploid context, not homozygous (1/1).

    This function detects such cases by parsing the consensusIDs on each SVcall and
    overrides any 1/1 genotype to 0/1 for the affected sample.
    """
    # Step 1: For each SVcall, collect (samplename, crID) → set of subIDs
    #         Also build an index from (samplename, crID) → list of SVcall indices

    # (samplename, crID) → set of subIDs seen across all SVcalls
    cr_subids: dict[tuple[str, int], set[int]] = defaultdict(set)
    # (samplename, crID) → list of SVcall indices that contain this (samplename, crID)
    cr_svcall_indices: dict[tuple[str, int], list[int]] = defaultdict(list)

    for idx, svcall in enumerate(svCalls):
        for cid in svcall.consensusIDs:
            parts = _parse_consensusID_parts(cid)
            if parts is None:
                continue
            samplename, crID, subID = parts
            key = (samplename, crID)
            cr_subids[key].add(subID)
            cr_svcall_indices[key].append(idx)

    # Step 2: Identify (samplename, crID) pairs with multiple assemblies
    multi_assembly_keys = {key for key, subids in cr_subids.items() if len(subids) >= 2}

    if not multi_assembly_keys:
        return svCalls

    # Step 3: For each affected SVcall + sample, override 1/1 → 0/1
    n_corrections = 0
    for key in multi_assembly_keys:
        samplename, crID = key
        for svcall_idx in cr_svcall_indices[key]:
            svcall = svCalls[svcall_idx]
            gt = svcall.genotypes.get(samplename)
            if gt is None:
                continue
            if gt.genotype == "1/1":
                log.info(
                    "GENOTYPE_CORRECTION|MULTI_ASSEMBLY|%s|crID=%d|subIDs=%s|%s: "
                    "GT 1/1 -> 0/1 (DV=%d, DR=%d, TC=%d)",
                    samplename,
                    crID,
                    sorted(cr_subids[key]),
                    svcall.to_log_id(),
                    gt.var_reads,
                    gt.ref_reads,
                    gt.total_coverage,
                )
                gt.genotype = "0/1"
                # Recompute GQ: we are confident this is het, but reflect that the
                # override is heuristic by assigning a moderate quality.
                gt.genotype_quality = min(gt.genotype_quality, 30)
                n_corrections += 1

    if n_corrections > 0:
        log.info(
            "GENOTYPE_CORRECTION|MULTI_ASSEMBLY|SUMMARY: corrected %d genotype(s) "
            "across %d multi-assembly loci",
            n_corrections,
            len(multi_assembly_keys),
        )

    return svCalls


# %%


def reference_bases_by_merged_svComposites(
    svComposites: list[SVcomposite],
    reference: Path,
    find_leftmost_reference_position: bool,
    tmp_dir_path: Path | str | None = None,
) -> dict[str, str]:
    """
    Retrieve reference bases for SVcomposite positions using streaming to avoid memory issues.
    """
    log.info(f"Collecting positions from {len(svComposites)} SVcomposites...")
    # print how much memory scComposites take
    from sys import getsizeof

    log.info(f"Size of svComposites in memory: {getsizeof(svComposites)} bytes")

    dict_reference_bases: dict[str, str] = {}

    with tempfile.TemporaryDirectory(
        dir=tmp_dir_path, delete=True if tmp_dir_path is None else False
    ) as temp_dir_str:
        temp_dir = Path(temp_dir_str)
        tmp_regions_path = temp_dir / "regions.txt"
        ref_bases_path = temp_dir / "ref_bases.txt"

        # 1. Collect unique positions and write to regions file formatted for samtools
        #    Use a set to avoid duplicates without invoking 'uniq' command
        positions = set()
        for svComposite in tqdm(svComposites, desc="Collecting positions"):
            chrname, start, end = get_svComposite_interval_on_reference(
                svComposite=svComposite,
                find_leftmost_reference_position=find_leftmost_reference_position,
            )
            # Store 0-based start
            positions.add((str(chrname), int(start)))

        # Sort positions by chromosome and then position
        sorted_positions = sorted(positions)

        with open(tmp_regions_path, "w") as f:
            for chrname, start_0 in sorted_positions:
                # Format: chr:start-end (1-based inclusive)
                # region start is the same as region end for a single base: chr:pos-pos
                pos_1 = start_0 + 1
                f.write(f"{chrname}:{pos_1}-{pos_1}\n")

        # 2. Retrieve reference bases with samtools
        cmd_faidx = f"samtools faidx --region-file {tmp_regions_path} {str(reference)}"
        log.info("Retrieving reference bases with samtools faidx...")

        with open(ref_bases_path, "w") as f_out:
            process_faidx = subprocess.Popen(
                split(cmd_faidx), stdout=f_out, stderr=subprocess.PIPE, text=True
            )
            _, stderr_faidx = process_faidx.communicate()
            if process_faidx.returncode != 0:
                log.error(
                    f"samtools faidx failed with return code {process_faidx.returncode}"
                )
                log.error(f"cmd_faidx: {cmd_faidx}")
                log.error(f"stderr: {stderr_faidx}")
                raise subprocess.CalledProcessError(
                    process_faidx.returncode, cmd_faidx, stderr=stderr_faidx
                )

        # 3. Parse the results
        log.info("Parsing retrieved reference bases...")
        with open(ref_bases_path, "r") as f:
            current_header = None
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    # Parse header line like ">chr1:123-123"
                    current_header = line[1:]  # Remove '>'
                elif current_header:
                    # This is the sequence line
                    if ":" in current_header and "-" in current_header:
                        chrom_pos_str = current_header.split("-")[
                            0
                        ]  # Get "chr1:123" part
                        try:
                            # Convert 1-based header coordinate back to 0-based key
                            r_chr, r_start_str = chrom_pos_str.rsplit(":", 1)
                            r_start_1 = int(r_start_str)
                            # Key format: chr:start (0-based)
                            key = f"{r_chr}:{r_start_1 - 1}"

                            base = line.upper()
                            if len(base) == 1:
                                dict_reference_bases[key] = base
                            else:
                                log.warning(
                                    f"Expected single base at {chrom_pos_str}, got '{base}'"
                                )
                        except ValueError:
                            log.warning(
                                f"Failed to parse header or key: {current_header}"
                            )

                    current_header = None  # Reset for next entry
    log.info(f"Successfully retrieved {len(dict_reference_bases)} reference bases")
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
    from ..signalprocessing.copynumber_tracks import load_copynumber_trees_from_db

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
    scale_by_complexity_factor: float = 1.0,
    verbose: bool = False,
    tmp_dir_path: Path | str | None = None,
    candidate_regions_file: Path | str | None = None,
    skip_covtrees: bool = False,
    collapse_repeats: bool = True,
) -> None:
    check_if_all_svtypes_are_supported(sv_types=sv_types)
    samplenames = [svirltile.get_metadata(Path(path))["samplename"] for path in input]

    # check if reference exists
    if not Path(reference).exists():
        raise FileNotFoundError(f"Reference file {reference} does not exist.")
    # check if the input files exist
    if not all(Path(path).exists() for path in input):
        missing_files = [str(path) for path in input if not Path(path).exists()]
        raise FileNotFoundError(
            f"The following input files do not exist: {', '.join(missing_files)}"
        )

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
    sv_types_set: set[type[SVpatterns.SVpatternType]] = set()
    for sv in sv_types:
        if sv in SUPPORTED_SV_TYPE_STRINGS_INVERSE:
            # Add all pattern types that map to this SV type string (e.g., "BND" -> [SVpatternSingleBreakend, SVpatternAdjacency])
            sv_types_set.update(SUPPORTED_SV_TYPE_STRINGS_INVERSE[sv])
    # debug - check what types are in sv_types_set
    if verbose:
        log.info(
            f"SV types to be processed: {', '.join([t.__name__ for t in sv_types_set])}"
        )

    data: list[SVcomposite] = generate_svComposites_from_dbs(
        input=input,
        sv_types=sv_types_set,
        candidate_regions_filter=candidate_regions_filter,
        collapse_repeats=collapse_repeats,
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
        # just save as pickle
        import pickle

        with open(Path(tmp_dir_path) / "all_svComposites.pkl", "wb") as f:
            pickle.dump(data, f)

    # --- vertical merging of svComposites across samples and consensus sequences --- #
    merged: list[SVcomposite] = merge_svComposites(
        apriori_size_difference_fraction_tolerance=apriori_size_difference_fraction_tolerance,
        svComposites=data,
        max_cohens_d=max_cohens_d,
        near=near,
        min_kmer_overlap=min_kmer_overlap,
        scale_by_complexity_factor=scale_by_complexity_factor,
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

    # Filter by minimum SV size and report what was dropped
    dropped_composites = []
    filtered_merged = []
    for svComposite in merged:
        if abs(svComposite.get_size()) >= min_sv_size:
            filtered_merged.append(svComposite)
        else:
            dropped_composites.append(svComposite)
            log.debug(
                f"DROPPED::multisample_sv_calling::MIN SV SIZE NOT REACHED: {abs(svComposite.get_size())} < min_sv_size {min_sv_size}, svComposite={_svcomposite_log_id(svComposite)}"
            )

    if dropped_composites:
        log.info(
            "SIZE_FILTER|SUMMARY|min_sv_size=%d|dropped=%d|kept=%d",
            min_sv_size,
            len(dropped_composites),
            len(filtered_merged),
        )

    merged = filtered_merged

    # if tmp dir is provided, dump all merged svComposites to a compressed json file
    if tmp_dir_path is not None:
        save_svComposites_to_json(
            data=merged, output_path=Path(tmp_dir_path) / "merged_svComposites.json.gz"
        )
    data.clear()  # free memory
    log.info(f"Generating SVcalls from {len(merged)} merged SVcomposites...")
    svCalls: list[SVcall] = [
        svcall
        for svComposite in merged
        for svcall in SVcalls_from_SVcomposite(
            svComposite,
            covtrees=covtrees,
            cn_tracks=cn_tracks,
            find_leftmost_reference_position=find_leftmost_reference_position,
            symbolic_threshold=symbolic_threshold,
        )
    ]

    # Correct genotypes at multi-assembly loci: when a candidate region produced
    # multiple consensus assemblies for a sample, the variants must be heterozygous.
    svCalls = correct_genotypes_for_multi_assembly_loci(svCalls)

    if tmp_dir_path is not None:
        save_svCalls_to_json(
            data=svCalls, output_path=Path(tmp_dir_path) / "svCalls.json.gz"
        )

    cn_tracks.clear()  # free memory
    log.info(
        f"Generated {len(svCalls)} SVcalls from {len(merged)} merged SVcomposites. Now adding ref bases.."
    )

    ref_bases_dict: dict[str, str] = reference_bases_by_merged_svComposites(
        svComposites=merged,
        reference=reference,
        find_leftmost_reference_position=find_leftmost_reference_position,
        tmp_dir_path=tmp_dir_path,
    )

    # save covtrees and ref_bases_dict to tmp dir if provided
    if tmp_dir_path is not None:
        import pickle

        with open(Path(tmp_dir_path) / "covtrees.pkl", "wb") as f:
            pickle.dump(covtrees, f)
        with open(Path(tmp_dir_path) / "ref_bases_dict.pkl", "wb") as f:
            pickle.dump(ref_bases_dict, f)

    log.info("Writing SVcalls to VCF file...")
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
    log_level = getattr(logging, args.log_level)
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    logfile = getattr(args, "logfile", None)
    if logfile:
        file_handler = logging.FileHandler(str(logfile), mode="w")
        file_handler.setLevel(logging.DEBUG)  # always capture full detail in file
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        handlers.append(file_handler)
    logging.basicConfig(
        level=min(log_level, logging.DEBUG) if logfile else log_level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=handlers,
        force=True,
    )
    # If a logfile is set but console level should stay at the requested level,
    # set only the console handler to the requested level.
    if logfile:
        handlers[0].setLevel(log_level)

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
        scale_by_complexity_factor=args.scale_by_complexity_factor,
        tmp_dir_path=args.tmp_dir_path,
        verbose=args.verbose,
        candidate_regions_file=args.candidate_regions_file,
        skip_covtrees=args.skip_covtrees,
        collapse_repeats=not args.dont_collapse_repeats,
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
        help="Size difference fraction tolerance for merging SVs (default: 0.1). Decrease for stronger separation of haplotypes",
        type=float,
        default=0.05,
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
        "--scale-by-complexity-factor",
        help="Weight (0.0 to 1.0) for scaling SV sizes by sequence complexity when merging. 0.0 disables complexity scaling, 1.0 applies full complexity scaling (default: 1.0).",
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "--dont-collapse-repeats",
        help="Disable merging of indels with the same repeatIDs during horizontal merge.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--verbose", help="Enable verbose output.", action="store_true", default=False
    )
    parser.add_argument(
        "--logfile",
        help="Path to log file. If provided, detailed decision logging (including region and crID annotations) is written to this file.",
        type=os.path.abspath,
        default=None,
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
        svcomposites = [SVcomposite.from_unstructured(item) for item in serialized_data]

        return svcomposites

    except json.JSONDecodeError as e:
        log.error(f"Failed to parse JSON from {input_path}: {e}")
        raise
    except Exception as e:
        log.error(f"Error loading SVcomposites from {input_path}: {e}")
        raise


def save_svCalls_to_json(data: list[SVcall], output_path: Path | str) -> None:
    """
    Save SVcalls to a compressed pickle file.

    Args:
        data: List of SVcalls to serialize
        output_path: Path to output pickle file

    Returns:
        None
    """
    output_path = Path(output_path)

    # Create parent directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save the data
    with open(output_path, "wb") as f:
        pickle.dump(data, f)

    log.info(f"Saved {len(data)} SVcalls to {output_path}")


def load_svCalls_from_json(input_path: Path | str) -> list[SVcall]:
    """
    Load SVcalls from a pickle file.

    Args:
        input_path: Path to input pickle file

    Returns:
        List of SVcalls
    """
    input_path = Path(input_path)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    try:
        with open(input_path, "rb") as f:
            svCalls = pickle.load(f)

        log.info(f"Loaded {len(svCalls)} SVcalls from {input_path}")
        return svCalls

    except pickle.UnpicklingError as e:
        log.error(f"Failed to unpickle SVcalls from {input_path}: {e}")
        raise
    except Exception as e:
        log.error(f"Error loading SVcalls from {input_path}: {e}")
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
        sv_types_set = set()
        for sv in sv_types:
            if sv in SUPPORTED_SV_TYPE_STRINGS_INVERSE:
                # Add all pattern types that map to this SV type string
                sv_types_set.update(SUPPORTED_SV_TYPE_STRINGS_INVERSE[sv])
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
