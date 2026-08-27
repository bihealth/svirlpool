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
from collections.abc import Iterable
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
})
#    SVpatterns.SVpatternAdjacency,
#    SVpatterns.SVpatternSingleBreakend,


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


# Ceiling for the phred-scaled genotype quality of *every* genotype this module
# emits, reference and alternate alike.  99 is the near-universal convention
# (GATK, bcftools, freebayes all cap GQ there), and it is where this caller's
# header now declares the ceiling.
#
# GQ is deliberately capped while `PL` is not: GQ is a single headline number in
# a conventional range, and the *uncapped* evidence strength that a consumer
# needs for ranking lives in PL, which is emitted alongside it.  Before v0.3 the
# cap was 60, which a clean heterozygote reached at about 18x -- so GQ stopped
# discriminating inside the depth range this caller actually operates in.
GQ_CEILING: int = 99

# The pre-v0.3 ceiling.  Kept only so that `--legacy-force-wildtype-genotypes`
# reproduces a pre-fix run's genotype fields exactly; it is not a modelling
# choice, it is the literal that used to be hard-coded in two places.
LEGACY_GQ_CEILING: int = 60

# Per-read error rate of the genotype model: the alternate-allele fraction
# expected at a homozygous-reference locus, and (as 1 - e) the fraction expected
# at a homozygous-alternate one.  Exposed as `--genotype-error-rate`.
#
# It is NOT only an error rate.  Because the genotype call is the argmax of
# binomial likelihoods evaluated at fixed expected fractions, this one number
# also fixes the allele-fraction boundaries between genotypes: with e = 0.05 a
# locus is called 0/0 below AF 0.218, 0/1 between 0.218 and 0.782, and 1/1 above
# 0.782, at every depth.  Raising e moves both boundaries inward, trading
# homozygote confidence against heterozygote confidence.  The default is
# unchanged from the value that was hard-coded before v0.3; it has deliberately
# *not* been fitted to any dataset.
DEFAULT_GENOTYPE_ERROR_RATE: float = 0.05

# Genotype string used when there is no evidence at all.  Diploid form, matching
# the default copy number assumed throughout this module.
NO_CALL_GENOTYPE: str = "./."

# Tolerance, in reference bases, for *breakpoint-placement disagreement* between
# the read-to-reference alignments that build the coverage tracks and the
# consensus-to-reference alignment that produces a call position.  The two place
# the same breakpoint independently, and in repetitive sequence they disagree.
#
# The value is the default of ``--near``, the distance within which
# ``merge_svComposites`` already declares two independently placed breakpoints to
# be the *same event*.  Refusing to count a read whose alignment ends 150 bases
# from the reported position, while merging two SVpatterns 150 bases apart into
# one variant, is internally inconsistent; one tolerance governs both.  The
# joint caller couples them by default (``--genotype-breakpoint-margin``
# defaults to ``--near``); this literal is the fallback for direct callers.
#
# Empirically (muc1 plus eight HG002 tiles, 21 loci, 42 supporting reads missed
# by a point query): 27 of the 34 recoverable reads are within 150 bases, at a
# cost of 1.7 bases of extra depth per read recovered; past 150 the cost triples
# for two more reads.  The previous value, 100, was undocumented and recovered
# 13.
DEFAULT_BREAKPOINT_MARGIN: int = 150

# The pre-N15 geometry: 100 bases at a deletion's two breakpoints, and no margin
# at all for an insertion.  Restored by ``--legacy-force-wildtype-genotypes`` so
# that the control arm reproduces a pre-fix run's genotype fields exactly.
LEGACY_BREAKPOINT_MARGIN: int = 100


# The FORMAT emitted by default, and the pre-v0.3 one restored by the legacy
# control.  They differ in more than one field, so the control has to select the
# whole string: pre-v0.3 GP was a single scalar (the posterior of the called
# genotype) rather than the per-genotype vector the VCF specification asks for.
FORMAT_FIELDS: str = "GT:GQ:TC:DR:DV:PL:GP"
LEGACY_FORMAT_FIELDS: str = "GT:GQ:TC:DR:DV:GP"


def genotype_alt_allele_count(genotype: str) -> int:
    """Number of alternate alleles in a genotype string ("0/1" -> 1)."""
    return sum(1 for allele in genotype.split("/") if allele not in ("0", "."))


def ordered_genotypes(genotypes: Iterable[str]) -> list[str]:
    """Genotypes in VCF ``Number=G`` order: by increasing alternate-allele count.

    For the biallelic case this module models, that is exactly the ordering the
    VCF specification defines for ``PL``/``GP`` (0/0, 0/1, 1/1 for a diploid).
    """
    return sorted(genotypes, key=genotype_alt_allele_count)


@attrs.define
class Genotype:
    samplename: str
    genotype: str  # e.g. "0/0", "0/1", "1/1", "./." for a no-call
    gt_likelihood: float | None  # posterior of the called genotype; None -> "."
    genotype_quality: int  # GQ phred of gt_likelihood
    total_coverage: int  # TC
    ref_reads: int  # DR
    var_reads: int  # DV
    # Per-genotype vectors, in VCF Number=G order.  None where the model was not
    # used (a no-call, the legacy force-call, or --single-evidence-gt), in which
    # case the fields are emitted as the VCF missing value.
    phred_likelihoods: dict[str, int] | None = None  # PL
    posteriors: dict[str, float] | None = None  # GP

    def _vector(self, values: dict[str, float] | None, fmt: str) -> str:
        if not values:
            return "."
        return ",".join(format(values[g], fmt) for g in ordered_genotypes(values))

    def to_format_field(self, FORMAT_field: str = FORMAT_FIELDS) -> str:
        properties = {
            "GT": self.genotype,
            "GQ": self.genotype_quality,
            "TC": self.total_coverage,
            "DR": self.ref_reads,
            "DV": self.var_reads,
            # Phred-scaled likelihoods, normalised so the best genotype is 0.
            # Uncapped on purpose: this is the field that still separates two
            # well-supported calls after GQ has hit its ceiling.
            "PL": self._vector(self.phred_likelihoods, "d"),
            # Posterior probabilities, one per genotype, summing to 1.  A
            # no-call has no genotype and therefore no posterior; "." is the VCF
            # missing value, and emitting 1.0 would be the very claim the record
            # is refusing to make.
            "GP": (
                self._vector(self.posteriors, ".6g")
                if self.posteriors
                else ("." if self.gt_likelihood is None else str(self.gt_likelihood))
            ),
        }
        return ":".join([str(properties[k]) for k in FORMAT_field.split(":")])


def phred_from_probability(probability: float, cap: int = GQ_CEILING) -> int:
    """Phred-scale a posterior probability, capped at *cap*.

    Retained for callers that hold a posterior rather than a likelihood vector.
    Note that it cannot see past the float64 resolution of ``1 - p``: the
    subtraction loses all precision once the posterior is within ~1e-16 of 1,
    which for this model happens at about 50x.  `genotype_quality_from_phred`
    computes the same quantity in log space and does not have that limit; it is
    what the genotyper uses.
    """
    if probability >= 1.0:
        return cap
    error = 1.0 - probability
    if error <= 0.0:
        return cap
    return min(int(-10 * np.log10(error)), cap)


def legacy_genotype_quality(probability: float) -> int:
    """The pre-v0.3 phred expression: uncapped, with a literal 60 at saturation.

    Kept verbatim, not as a special case of `phred_from_probability`, because
    ``--legacy-force-wildtype-genotypes`` exists to reproduce a pre-fix run's
    genotype fields *exactly*.  It is the expression this module used on the
    alternate-supported path, defects included: unbounded above (GQ=154 was
    observed at chr11:11,246,978 with DV=23/TC=50), and non-monotone, because a
    posterior that saturates to exactly 1.0 in float64 drops to the literal 60
    while a weaker one nearby scores far higher.
    """
    if probability < 1.0:
        return int(-10 * np.log10(1.0 - probability))
    return LEGACY_GQ_CEILING


def reference_genotype_string(copy_number: int) -> str:
    """The all-reference genotype string for a locus of the given copy number."""
    if copy_number <= 1:
        return "0"
    return "/".join(["0"] * copy_number)


def create_no_call_genotype(samplename: str) -> Genotype:
    """A genotype that refuses to make a call because there is no evidence.

    Emitted where the coverage tracks show that the locus is not covered by any
    read of this sample: a technical dropout must not be reported as an observed
    reference allele.  TC/DR/DV are 0 because nothing was observed, GQ is 0
    because no genotype was called, and GP is missing for the same reason.
    """
    return Genotype(
        samplename=samplename,
        genotype=NO_CALL_GENOTYPE,
        gt_likelihood=None,
        genotype_quality=0,
        total_coverage=0,
        ref_reads=0,
        var_reads=0,
    )


def create_wild_type_genotype(
    samplename: str,
    total_coverage: int,
    copy_number: int = 2,
    legacy_force_call: bool = False,
    error_rate: float = DEFAULT_GENOTYPE_ERROR_RATE,
) -> Genotype:
    """A homozygous-reference call justified by *total_coverage* observed reads.

    Scored by the same model that scores every other genotype -- the binomial
    likelihoods with zero alternate reads -- so the quality grows with depth and
    the record carries the same PL/GP vectors as an alternate-supported call.

    With *legacy_force_call* the pre-v0.3 behaviour is reproduced exactly:
    ``0/0`` with a scalar ``GP = 1.0`` and ``GQ = 60`` regardless of depth,
    including zero depth.  See ``--legacy-force-wildtype-genotypes``.
    """
    if legacy_force_call:
        return Genotype(
            samplename=samplename,
            genotype="0/0",
            gt_likelihood=1.0,
            genotype_quality=LEGACY_GQ_CEILING,
            total_coverage=total_coverage,
            ref_reads=total_coverage,
            var_reads=0,
        )
    if total_coverage <= 0:
        # No observations: absence of evidence is not evidence of the reference.
        return create_no_call_genotype(samplename=samplename)
    log_likelihoods = genotype_log_likelihoods(
        n_alt_reads=0,
        n_total_reads=total_coverage,
        cn=copy_number,
        error_rate=error_rate,
    )
    phred = phred_scaled_likelihoods(log_likelihoods)
    posteriors = posteriors_from_log_likelihoods(log_likelihoods)
    ref_gt = reference_genotype_string(copy_number)
    return Genotype(
        samplename=samplename,
        genotype=ref_gt,
        gt_likelihood=posteriors.get(ref_gt, 0.0),
        genotype_quality=genotype_quality_from_phred(phred),
        total_coverage=total_coverage,
        ref_reads=total_coverage,
        var_reads=0,
        phred_likelihoods=phred,
        posteriors=posteriors,
    )


def create_no_evidence_genotype(
    samplename: str, legacy_force_wildtype: bool = False
) -> Genotype:
    """No-call, unless the legacy force-call control is switched on."""
    if legacy_force_wildtype:
        return create_wild_type_genotype(
            samplename=samplename, total_coverage=0, legacy_force_call=True
        )
    return create_no_call_genotype(samplename=samplename)


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
        return f"{self.svtype}|{self.chrname}:{self.start}-{self.end}|consensusIDs={','.join(sorted(self.consensusIDs))}|svlen={self.svlen}|mateid={self.mateid}"

    def to_vcf_line(
        self,
        vcfIDnumber: int,
        samplenames: list[str],
        covtrees: dict[str, dict[str, IntervalTree]],
        refdict: dict[str, str],
        symbolic_threshold: int,
        ONE_BASED: int = 1,
        legacy_force_wildtype: bool = False,
        error_rate: float = DEFAULT_GENOTYPE_ERROR_RATE,
    ) -> str | None:
        info_fields: dict = {
            "PASS_ALTREADS": self.pass_altreads,
            "pass_GQ": self.pass_gq,
            "SVTYPE": self.svtype,
            "END": str(
                self.end + ONE_BASED
            ),  # end is inclusive in vcf, so we need to subtract 1 -> 'END':str(self.end+ONE_BASED-1)
            "SVLEN": str(self.svlen),
            # sorted: CONSENSUSIDs is an unordered Number=. list, and a
            # canonical order is what makes two runs byte-comparable.
            "CONSENSUSIDs": ",".join(sorted(self.consensusIDs)),
            "MATEID": self.mateid if len(self.mateid) > 0 else "NA",
        }

        # Add sequence ID if this is a symbolic allele
        if self.sequence_id:
            info_fields["SEQ_ID"] = self.sequence_id

        info_line = ";".join([f"{key}={value}" for key, value in info_fields.items()])
        info_line += ";" + ("PRECISE" if self.precise else "IMPRECISE")
        FORMAT_field = LEGACY_FORMAT_FIELDS if legacy_force_wildtype else FORMAT_FIELDS
        format_content: list[Genotype] = []
        for samplename in samplenames:
            gt: Genotype | None = self.genotypes.get(samplename, None)
            if gt is None:
                # This sample contributed no SVpattern to the composite, so no
                # Genotype was computed for it.  Whether that means "reference"
                # or "no data" is decided by the coverage tracks, not assumed.
                tree: IntervalTree | None = covtrees.get(samplename, {}).get(
                    self.chrname, None
                )
                found_intervals: set[Interval] = tree[self.start] if tree else set()
                # Count distinct reads, not alignment fragments: a single read
                # contributes several RAF intervals at a locus it spans, so
                # len(found_intervals) over-counts depth (78 reads were reported
                # as TC=210 at the 18 kb chr6 insertion).  This matches what
                # get_ref_reads_from_covtrees returns for every other sample.
                coverage: int = len({it.data for it in found_intervals})
                if legacy_force_wildtype:
                    # The control must reproduce the pre-fix record byte for
                    # byte, including the over-counted depth.
                    format_content.append(
                        create_wild_type_genotype(
                            samplename=samplename,
                            total_coverage=len(found_intervals),
                            legacy_force_call=True,
                        )
                    )
                elif coverage == 0:
                    format_content.append(
                        create_no_call_genotype(samplename=samplename)
                    )
                else:
                    format_content.append(
                        create_wild_type_genotype(
                            samplename=samplename,
                            total_coverage=coverage,
                            error_rate=error_rate,
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

            if not use_symbolic and ref_seq == alt_seq and (ref_seq or alt_seq):
                # REF == ALT must never reach the VCF: htsjdk/IGV (and bcftools)
                # reject this as "Duplicate allele added to VariantContext".  This
                # typically indicates a bug in how the alt sequence is derived for
                # the given SV type (e.g. missing reverse-complement for an
                # inversion whose consensus slice matches the reference strand).
                raise ValueError(
                    "to_vcf_line: REF sequence is identical to ALT sequence — "
                    "this would produce an invalid VCF record.\n"
                    f"  locus              : {self.chrname}:{self.start}-{self.end} (1-based start: {self.start + ONE_BASED})\n"
                    f"  svtype             : {self.svtype}\n"
                    f"  svlen              : {self.svlen}\n"
                    f"  vcfID              : {vcfID}\n"
                    f"  consensusIDs       : {','.join(self.consensusIDs)}\n"
                    f"  mateid             : {self.mateid or 'NA'}\n"
                    f"  symbolic_threshold : {symbolic_threshold}\n"
                    f"  ref_seq length     : {len(ref_seq)}\n"
                    f"  alt_seq length     : {len(alt_seq)}\n"
                    f"  ref_seq[:80]       : {ref_seq[:80]!r}\n"
                    f"  alt_seq[:80]       : {alt_seq[:80]!r}\n"
                    "  Hint: for inversions, the alt allele must be returned in "
                    "reference-forward orientation as the reverse complement of "
                    "the assembled consensus slice."
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
    margin: int = 0,
) -> set[int]:
    """Read-name hashes of every alignment fragment overlapping the queried window.

    ``min_radius`` is a *minimum width*: a window narrower than ``2 * min_radius``
    is grown symmetrically until it reaches that width.  ``margin`` is a genuine
    margin: it is added to both ends whatever the width, so a read whose nearest
    alignment fragment ends at most ``margin`` bases before ``start`` (or begins
    at most ``margin`` bases after ``end``) is still reported.  See
    ``genotype_of_sample`` for why a breakpoint query needs one.
    """
    if end < start:
        raise ValueError(
            f"get_ref_reads_from_covtrees: end {end} is less than start {start} for sample {samplename} at {chrname}:{start}-{end}."
        )
    margin = abs(margin)
    if margin:
        start = max(0, start - margin)
        end = end + margin
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


def _single_evidence_genotype(
    n_alt_reads: int,
    n_ref_reads: int,
    n_total_reads: int,
    copy_number: int,
) -> tuple[str, float]:
    """Assign a genotype based solely on read presence, ignoring probabilistic model.

    Any alt read(s) make the call; complete absence of noise is assumed.
    Returns the genotype string and a fixed likelihood of 1.0.
    """
    if copy_number == 0:
        return "0", 1.0  # Homozygous reference if no copies are present
    if n_alt_reads == 0:
        gt = "0/0" if copy_number >= 2 else "0"
    elif copy_number == 1:
        gt = "1"
    elif copy_number == 2:
        gt = "1/1" if n_ref_reads == 0 else "0/1"
    else:
        # Higher CN: count alt copies by rounding alt fraction to nearest integer
        n_alt_copies = (
            round(n_alt_reads / n_total_reads * copy_number) if n_total_reads > 0 else 0
        )
        n_alt_copies = max(0, min(n_alt_copies, copy_number))
        gt = "/".join(["1" if i < n_alt_copies else "0" for i in range(copy_number)])
    return gt, 1.0


def copy_number_at_locus(
    samplename: str,
    chrname: str,
    start: int,
    end: int,
    cn_tracks: dict[str, dict[str, IntervalTree]],
) -> int:
    """Query the sample's copy-number track at the locus; 2 (diploid) if unknown."""
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
    return copy_number


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
    breakpoint_margin: int = DEFAULT_BREAKPOINT_MARGIN,
    apply_breakpoint_margin: bool = False,
    single_evidence_gt: bool = False,
    legacy_force_wildtype: bool = False,
    error_rate: float = DEFAULT_GENOTYPE_ERROR_RATE,
) -> Genotype:
    """Genotype one sample at one locus from its coverage tracks.

    The coverage query has to tolerate *breakpoint-placement disagreement*.  A
    read that carries the event has, by construction, sequence at the breakpoint
    that does not align to the reference there; the aligner emits one fragment
    ending before the breakpoint and another resuming after it, and no fragment
    covering it.  The breakpoint is then placed twice and independently — once
    by the read-to-reference alignment that built the coverage track, once by
    the consensus-to-reference alignment that produced the call position — and
    in repetitive sequence the two placements differ.  A point query at the call
    position therefore selects *against* the very reads that support the call.

    ``breakpoint_mode`` (deletions) queries the two breakpoints separately
    rather than the deleted span.  ``apply_breakpoint_margin`` (insertions,
    which have a single breakpoint and no span to query) widens the locus query
    by ``breakpoint_margin`` on both sides.  Both use the same tolerance.
    """
    genotype: Genotype
    if chrname not in covtrees.get(samplename, {}):
        log.warning(
            f"SVcalls_from_SVcomposite: Chromosome {chrname} not found in coverage tree "
            f"for sample {samplename}. Emitting a no-call ({NO_CALL_GENOTYPE})."
        )
        return create_no_evidence_genotype(
            samplename=samplename, legacy_force_wildtype=legacy_force_wildtype
        )
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
        # For insertions the reference footprint is a single point, so there is
        # no interior to query and `breakpoint_mode`'s two-window form would be
        # identical to one widened window.  Widening keeps the query a strict
        # superset of the un-margined one: no read a point query found can be
        # lost.
        all_reads: set[int] = get_ref_reads_from_covtrees(
            samplename=samplename,
            chrname=chrname,
            start=start,
            end=end,
            covtrees=covtrees,
            min_radius=min_radius,
            margin=breakpoint_margin if apply_breakpoint_margin else 0,
        )
    alt_reads: set[int] = all_reads.intersection(raw_alt_reads)

    ref_reads: set[int] = all_reads.difference(alt_reads)

    copy_number: int = copy_number_at_locus(
        samplename=samplename,
        chrname=chrname,
        start=start,
        end=end,
        cn_tracks=cn_tracks,
    )

    if len(alt_reads) == 0:
        # `alt_reads` is `all_reads & raw_alt_reads`, so an empty `all_reads`
        # necessarily lands here too.  The former separate `len(all_reads) == 0`
        # branch below this one was unreachable; it is folded in as the explicit
        # sub-case that follows.
        log.warning(
            f"SVcalls_from_SVcomposite: No alt reads found for sample {samplename} at {chrname}:{start}-{end}."
        )
        if len(all_reads) == 0:
            # Genuinely no coverage: a technical dropout, not an observation of
            # the reference allele.  Do not claim a genotype.
            log.warning(
                f"SVcalls_from_SVcomposite: No reads at all for sample {samplename} at "
                f"{chrname}:{start}-{end}. Emitting a no-call ({NO_CALL_GENOTYPE})."
            )
            return create_no_evidence_genotype(
                samplename=samplename, legacy_force_wildtype=legacy_force_wildtype
            )
        # Coverage present and no alternate support: the reference call is
        # justified, with full reference support and a depth-derived quality.
        return create_wild_type_genotype(
            samplename=samplename,
            total_coverage=len(all_reads),
            copy_number=copy_number,
            legacy_force_call=legacy_force_wildtype,
            error_rate=error_rate,
        )

    # Compute genotype either via the probabilistic binomial model or the
    # single-evidence rule (any alt read counts; no noise assumed).
    phred: dict[str, int] | None = None
    posteriors: dict[str, float] | None = None
    if single_evidence_gt:
        # Read presence alone decides the call; no model was evaluated, so there
        # are no per-genotype likelihoods to report.
        gt, gt_likelihood_val = _single_evidence_genotype(
            n_alt_reads=len(alt_reads),
            n_ref_reads=len(ref_reads),
            n_total_reads=len(all_reads),
            copy_number=copy_number,
        )
        gt_likelihoods: dict[str, float] = {gt: gt_likelihood_val}
        genotype_quality = GQ_CEILING
    elif legacy_force_wildtype:
        # The control reproduces a pre-fix run field for field: the scalar GP,
        # the hard-coded error rate, and the uncapped quality expression.
        gt_likelihoods = legacy_genotype_posteriors(
            n_alt_reads=len(alt_reads),
            n_total_reads=len(all_reads),
            cn=copy_number,
        )
        gt = max(gt_likelihoods.items(), key=lambda x: x[1])[0]
        genotype_quality = legacy_genotype_quality(gt_likelihoods[gt])
    else:
        log_likelihoods = genotype_log_likelihoods(
            n_alt_reads=len(alt_reads),
            n_total_reads=len(all_reads),
            cn=copy_number,
            error_rate=error_rate,
        )
        phred = phred_scaled_likelihoods(log_likelihoods)
        posteriors = posteriors_from_log_likelihoods(log_likelihoods)
        gt_likelihoods = posteriors
        # The call is the most likely genotype -- equivalently the one with
        # PL 0 -- and the quality is how far ahead of the runner-up it is, in
        # log space.  Both exits of this function use that one scale, so a
        # reference call and an alternate call are directly comparable.
        gt = min(phred.items(), key=lambda item: item[1])[0]
        genotype_quality = genotype_quality_from_phred(phred)

    genotype = Genotype(
        samplename=samplename,
        genotype=gt,
        gt_likelihood=gt_likelihoods[gt],
        genotype_quality=genotype_quality,
        total_coverage=len(all_reads),
        ref_reads=len(ref_reads),
        var_reads=len(alt_reads),
        phred_likelihoods=phred,
        posteriors=posteriors,
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
    single_evidence_gt: bool = False,
    min_alt_reads: int = 3,
    breakpoint_margin: int = DEFAULT_BREAKPOINT_MARGIN,
    legacy_force_wildtype: bool = False,
    error_rate: float = DEFAULT_GENOTYPE_ERROR_RATE,
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
            single_evidence_gt=single_evidence_gt,
            min_alt_reads=min_alt_reads,
            breakpoint_margin=breakpoint_margin,
            legacy_force_wildtype=legacy_force_wildtype,
            error_rate=error_rate,
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
            single_evidence_gt=single_evidence_gt,
            min_alt_reads=min_alt_reads,
            legacy_force_wildtype=legacy_force_wildtype,
            error_rate=error_rate,
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
    single_evidence_gt: bool = False,
    min_alt_reads: int = 3,
    breakpoint_margin: int = DEFAULT_BREAKPOINT_MARGIN,
    legacy_force_wildtype: bool = False,
    error_rate: float = DEFAULT_GENOTYPE_ERROR_RATE,
) -> SVcall:
    chrname, start, end = get_svComposite_interval_on_reference(
        svComposite=svComposite,
        find_leftmost_reference_position=find_leftmost_reference_position,
    )
    svlen: int = abs(svComposite.get_size())
    # Get read IDs supporting the SV. arguments: samplename, chrname, start, end, covtree
    consensusIDs: list[str] = sorted({
        svPattern.samplenamed_consensusID for svPattern in svComposite.svPatterns
    })

    # Both indels need the breakpoint-placement tolerance, in the geometry their
    # SV type dictates: a deletion has two distant breakpoints and an interior
    # that carries no alt-read alignments, so each breakpoint is queried
    # separately; an insertion has a single breakpoint and no interior, so the
    # point query is simply widened.  Only the deletion half of this was wired
    # up (N15), which made an insertion's own supporting reads invisible to the
    # genotyper.
    is_deletion = issubclass(svComposite.sv_type, SVpatterns.SVpatternDeletion)
    is_insertion = issubclass(svComposite.sv_type, SVpatterns.SVpatternInsertion)
    if legacy_force_wildtype:
        # The control has to reproduce a pre-fix run's genotype fields exactly,
        # which means restoring the *asymmetric* pre-fix geometry: no margin for
        # an insertion, 100 bases for a deletion.  No single value of
        # --genotype-breakpoint-margin can express that.
        is_insertion = False
        breakpoint_margin = LEGACY_BREAKPOINT_MARGIN
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
            apply_breakpoint_margin=is_insertion,
            breakpoint_margin=breakpoint_margin,
            single_evidence_gt=single_evidence_gt,
            legacy_force_wildtype=legacy_force_wildtype,
            error_rate=error_rate,
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
        max(genotypes.items(), key=lambda x: x[1].var_reads)[1].var_reads
        >= min_alt_reads
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
    single_evidence_gt: bool = False,
    min_alt_reads: int = 3,
    legacy_force_wildtype: bool = False,
    error_rate: float = DEFAULT_GENOTYPE_ERROR_RATE,
) -> list[SVcall]:
    """Generate SVcall objects from a SVcomposite that represents novel adjacencies with two connected break ends of each sample.
    The given svComposite generates two SVcall objects, that both represent one end of the novel adjacency.
    """
    if not issubclass(svComposite.sv_type, SVpatterns.SVpatternAdjacency):
        raise ValueError(
            f"svcall_objects_from_Adjacencies: svComposite of type {svComposite.sv_type.__name__} is not of type SVpatternAdjacency."
        )

    # Collect consensusIDs from all SVpatterns in the composite
    consensusIDs: list[str] = sorted({
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
            single_evidence_gt=single_evidence_gt,
            legacy_force_wildtype=legacy_force_wildtype,
            error_rate=error_rate,
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
            single_evidence_gt=single_evidence_gt,
            legacy_force_wildtype=legacy_force_wildtype,
            error_rate=error_rate,
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
    pass_altreads: bool = max(max_var_reads_0, max_var_reads_1) >= min_alt_reads

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
    reference: Path,
    samplenames: list[str],
    fasta_path: Path | None = None,
    legacy_force_wildtype: bool = False,
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
    header.append(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype; ./. means no '
        'read coverage at this locus for this sample, i.e. no call, not reference">'
    )
    header.append(
        f'##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality of '
        f"the called genotype: the phred-scaled likelihood ratio against the next "
        f"most likely genotype, capped at {GQ_CEILING}. See PL for the uncapped "
        f'values; 0 for a no-call">'
    )
    header.append('##FORMAT=<ID=TC,Number=1,Type=Integer,Description="Total coverage">')
    header.append(
        '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">'
    )
    header.append(
        '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">'
    )
    if legacy_force_wildtype:
        header.append(
            '##FORMAT=<ID=GP,Number=1,Type=Float,Description="Posterior probability of '
            'the called genotype; missing (.) for a no-call">'
        )
    else:
        header.append(
            '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype '
            "likelihoods, normalised so the most likely genotype is 0, in the order "
            "0/0,0/1,1/1 for a diploid locus. Not capped: this is the field that "
            "still separates two well-supported calls after GQ reaches its ceiling. "
            'Missing (.) for a no-call">'
        )
        header.append(
            '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype posterior '
            "probabilities under a uniform prior, summing to 1, in the same order as "
            'PL. Missing (.) for a no-call">'
        )

    header.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samplenames)
    )
    return header


def expected_alt_fractions(
    cn: int, error_rate: float = DEFAULT_GENOTYPE_ERROR_RATE
) -> dict[str, float]:
    """Genotype string -> the alternate-read fraction expected under it.

    One formula for every copy number.  ``n`` alternate copies out of ``cn``
    predict an alternate fraction of ``n / cn``, clamped into
    ``[error_rate, 1 - error_rate]`` so that the two extremes admit read noise.
    Before v0.3 the CN=1..4 cases were written out by hand; the clamped general
    expression reproduces all four exactly, and the hand-written copies had
    drifted apart in style without drifting in value.
    """
    if cn <= 0:
        # Copy number 0 (homozygous deletion) admits a single genotype and no
        # reads.  `genotype_likelihood` short-circuits it; this keeps the shape
        # well defined for any caller that gets here anyway.
        return {"0": error_rate}
    fractions: dict[str, float] = {}
    for n_alt_copies in range(cn + 1):
        # Alleles ascending ("0/1", not "1/0"): the VCF convention, and what the
        # hand-written CN=2..4 tables used.  The pre-v0.3 CN>4 branch built them
        # the other way round and so emitted non-canonical genotype strings.
        genotype = "/".join(["0"] * (cn - n_alt_copies) + ["1"] * n_alt_copies)
        fractions[genotype] = min(max(n_alt_copies / cn, error_rate), 1.0 - error_rate)
    return fractions


def genotype_log_likelihoods(
    n_alt_reads: int,
    n_total_reads: int,
    cn: int,
    error_rate: float = DEFAULT_GENOTYPE_ERROR_RATE,
) -> dict[str, float]:
    """Natural-log binomial likelihood of the observed counts under each genotype.

    Everything downstream of the model works from these logs.  Likelihoods
    themselves underflow (0.05**n reaches the float64 floor by n=220) and, worse,
    a posterior formed from them saturates to exactly 1.0 by about 50x, taking
    the genotype quality with it.  In log space neither happens.
    """
    return {
        genotype: float(binom.logpmf(n_alt_reads, n_total_reads, p))
        for genotype, p in expected_alt_fractions(cn=cn, error_rate=error_rate).items()
    }


# Value substituted for a non-finite phred-scaled likelihood.  Reachable only
# for impossible counts (n_alt > n_total); a finite sentinel keeps the VCF
# parseable rather than emitting "inf".
PL_SENTINEL: int = 99999


def phred_scaled_likelihoods(log_likelihoods: dict[str, float]) -> dict[str, int]:
    """PL: phred-scaled likelihoods, shifted so the most likely genotype is 0.

    Deliberately **uncapped**.  GQ is a headline number in a conventional range
    and saturates; PL is the field that still separates two well-supported calls
    afterwards, and it grows roughly linearly with depth.
    """
    best = max(log_likelihoods.values())
    scale = 10.0 / np.log(10.0)
    phred: dict[str, int] = {}
    for genotype, log_likelihood in log_likelihoods.items():
        value = (best - log_likelihood) * scale
        phred[genotype] = int(round(value)) if np.isfinite(value) else PL_SENTINEL
    return phred


def posteriors_from_log_likelihoods(
    log_likelihoods: dict[str, float],
) -> dict[str, float]:
    """Flat-prior posterior probabilities, by log-sum-exp.

    The prior over genotypes is uniform: this caller makes no allele-frequency
    or Hardy-Weinberg assumption.  Subtracting the maximum before exponentiating
    is what makes this safe at any depth; the pre-v0.3 code divided raw
    ``binom.pmf`` values, which could underflow to an all-zero vector and fall
    back to a silent uniform.
    """
    best = max(log_likelihoods.values())
    weights = {
        genotype: float(np.exp(log_likelihood - best))
        for genotype, log_likelihood in log_likelihoods.items()
    }
    total = sum(weights.values())
    return {genotype: weight / total for genotype, weight in weights.items()}


def genotype_quality_from_phred(
    phred_likelihoods: dict[str, int], cap: int = GQ_CEILING
) -> int:
    """GQ: the phred-scaled likelihood ratio between the best genotype and the next.

    This is the conventional definition (the difference between the two smallest
    PL values) and it is the same quantity the pre-v0.3 code was reaching for
    with ``-10 log10(1 - P(best))`` -- the two agree to within the 3 dB that
    counting one runner-up rather than all of them costs.  The difference is that
    this one is computed in log space, so it neither saturates nor inverts.
    """
    values = sorted(phred_likelihoods.values())
    if len(values) < 2:
        return cap
    return int(min(values[1] - values[0], cap))


def legacy_genotype_posteriors(
    n_alt_reads: int, n_total_reads: int, cn: int
) -> dict[str, float]:
    """The pre-v0.3 posterior arithmetic, reproduced verbatim.

    Raw ``binom.pmf`` values divided by their sum, rather than the log-sum-exp
    the fixed path uses.  The two agree to within one unit in the last place,
    but ``legacy_genotype_quality`` amplifies exactly that last place: at
    DV=23/TC=50 it is the difference between GQ 154 and GQ 153.  The control
    exists to reproduce a pre-fix run *exactly*, so it has to reproduce the
    arithmetic and not merely the formula.
    """
    fractions = expected_alt_fractions(cn=cn, error_rate=DEFAULT_GENOTYPE_ERROR_RATE)
    likelihoods = {
        genotype: binom.pmf(n_alt_reads, n_total_reads, p)
        for genotype, p in fractions.items()
    }
    total_likelihood = sum(likelihoods.values())
    if total_likelihood <= 0.0:
        return dict.fromkeys(likelihoods, 1.0 / len(likelihoods))
    probabilities = {
        genotype: float(likelihood / total_likelihood)
        for genotype, likelihood in likelihoods.items()
    }
    return {
        genotype: (0.0 if np.isnan(probability) else probability)
        for genotype, probability in probabilities.items()
    }


def genotype_likelihood(
    n_alt_reads: int,
    n_total_reads: int,
    cn: int,
    legacy_zero_coverage: bool = False,
    error_rate: float = DEFAULT_GENOTYPE_ERROR_RATE,
) -> dict[str, float]:
    """
    Posterior probability of each genotype given observed alt/total reads and CN.

    Args:
        n_alt_reads: Number of reads supporting the alternate allele
        n_total_reads: Total number of reads covering the locus
        cn: Copy number at this locus (1, 2, 3, 4, etc.)
        legacy_zero_coverage: Reproduce the pre-v0.3 behaviour of returning
            certainty of ``0/0`` when nothing at all was observed.  See
            ``--legacy-force-wildtype-genotypes``.
        error_rate: Alternate fraction expected at a homozygous-reference locus.
            See ``--genotype-error-rate``; this also sets the allele-fraction
            boundaries between genotypes.

    Returns:
        Dictionary mapping genotype strings to their posterior probabilities.

    With no observations at all (``n_total_reads == 0``) the posterior is
    *uniform* over the genotypes admissible at this copy number: zero
    observations support no genotype over any other.  Callers must treat that as
    "no information" and emit a no-call rather than taking the argmax.

    Note that a posterior cannot express confidence beyond the float64
    resolution of 1.0, which this model reaches at about 50x.  Anything that
    needs to rank calls above that should use `phred_scaled_likelihoods`.
    """
    if n_total_reads == 0 and legacy_zero_coverage:
        return {"0/0": 1.0, "0/1": 0.0, "1/1": 0.0}

    if cn == 0:
        # Homozygous deletion - no copies, so no reads and only one genotype.
        return {"0": 1.0}

    fractions = expected_alt_fractions(cn=cn, error_rate=error_rate)

    if n_total_reads == 0:
        # No observations: report an uninformative (uniform) posterior rather
        # than certainty of the reference genotype.
        log.debug(
            f"genotype_likelihood: no reads observed (CN={cn}); returning a uniform, "
            "uninformative posterior."
        )
        return dict.fromkeys(fractions, 1.0 / len(fractions))

    return posteriors_from_log_likelihoods(
        genotype_log_likelihoods(
            n_alt_reads=n_alt_reads,
            n_total_reads=n_total_reads,
            cn=cn,
            error_rate=error_rate,
        )
    )


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
    #         sorted(): the emitted GENOTYPE_CORRECTION log lines are used for
    #         run-to-run diffing, so their order must not depend on set
    #         iteration over (samplename, crID) tuples.
    n_corrections = 0
    for key in sorted(multi_assembly_keys):
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
    legacy_force_wildtype: bool = False,
    error_rate: float = DEFAULT_GENOTYPE_ERROR_RATE,
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
            legacy_force_wildtype=legacy_force_wildtype,
            reference=Path(reference),
            samplenames=samplenames,
            fasta_path=fasta_path,
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
                legacy_force_wildtype=legacy_force_wildtype,
                error_rate=error_rate,
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
    single_evidence_gt: bool = False,
    min_alt_reads: int = 3,
    genotype_breakpoint_margin: int | None = None,
    legacy_force_wildtype: bool = False,
    error_rate: float = DEFAULT_GENOTYPE_ERROR_RATE,
) -> None:
    check_if_all_svtypes_are_supported(sv_types=sv_types)
    if not 0.0 < error_rate < 0.5:
        raise ValueError(
            f"--genotype-error-rate must lie strictly between 0 and 0.5, got "
            f"{error_rate}. It is the alternate fraction expected at a "
            f"homozygous-reference locus; at 0.5 the reference and heterozygous "
            f"hypotheses become identical and the model cannot separate them."
        )
    if error_rate != DEFAULT_GENOTYPE_ERROR_RATE:
        log.warning(
            "--genotype-error-rate=%s moves the allele-fraction boundaries between "
            "genotypes away from their defaults; genotypes are not comparable with "
            "a run at the default of %s.",
            error_rate,
            DEFAULT_GENOTYPE_ERROR_RATE,
        )
    # The genotyping window's tolerance for breakpoint-placement disagreement.
    # It is the same notion as `--near`, so it follows `--near` unless the user
    # separates them deliberately.
    breakpoint_margin: int = (
        near if genotype_breakpoint_margin is None else (genotype_breakpoint_margin)
    )
    log.info(
        "Genotyping breakpoint-placement margin: %d bp (%s)",
        breakpoint_margin,
        "from --near" if genotype_breakpoint_margin is None else "explicit",
    )
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
        log.warning(
            "--skip-covtrees: coverage tracks are fabricated and carry no information "
            "about the samples. Every depth-derived field (TC, DR, DV, GQ) in the "
            "output is meaningless. Investigation use only; never for a released VCF."
        )
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
            single_evidence_gt=single_evidence_gt,
            min_alt_reads=min_alt_reads,
            breakpoint_margin=breakpoint_margin,
            legacy_force_wildtype=legacy_force_wildtype,
            error_rate=error_rate,
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
        legacy_force_wildtype=legacy_force_wildtype,
        error_rate=error_rate,
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
        single_evidence_gt=args.single_evidence_gt,
        min_alt_reads=args.min_alt_reads,
        genotype_breakpoint_margin=args.genotype_breakpoint_margin,
        legacy_force_wildtype=args.legacy_force_wildtype_genotypes,
        error_rate=args.genotype_error_rate,
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
        "--min-alt-reads",
        help="Minimum number of variant-supporting reads for a call to PASS the "
        "quality filter (default: 3). Calls below this are emitted with FILTER=LowQual. "
        "Lower (e.g. 2) to increase recall on low-coverage or noisy loci at some cost "
        "to precision.",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--genotype-breakpoint-margin",
        help="Tolerance, in reference bases, for disagreement between where the "
        "read alignments place a breakpoint and where the consensus alignment "
        "places it, when querying the coverage tracks for a genotype. Reads that "
        "carry an indel have no alignment across its breakpoint, so a point query "
        "at the call position systematically misses them. Defaults to --near, "
        "which is the same tolerance applied when merging SVpatterns into one "
        "variant. Set to 0 to query the exact call position only.",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--apriori-size-difference-fraction-tolerance",
        help="Fraction of the larger of the two sizes that two SVs may "
        "differ by and still be merged (default: 0.06). 0.0 = no tolerance (sizes must "
        "be identical); 1.0 = maximum tolerance, which is vacuous by construction and "
        "reproduces the pre-fix behaviour of an inert size gate. Decrease for stronger "
        "separation of haplotypes.",
        type=float,
        default=0.06,
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
        help=(
            "Weight (0.0 to 1.0) on the size tolerance granted by low sequence "
            "complexity when merging. Low-complexity sequence gives the aligner more "
            "freedom in where it places an indel and how large it calls it, so the "
            "two size populations' means are shifted toward each other by "
            "weight * (1 - mean_complexity) * |size| before Cohen's d is computed. "
            "0.0 grants no complexity allowance; 1.0 grants it in full (default: 1.0)."
        ),
        type=float,
        default=1.0,
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
        "--single-evidence-gt",
        help="Experimental: assign genotypes based solely on read presence rather than a probabilistic model. "
        "Any alt read(s) call the variant; any ref read(s) alongside alt make it heterozygous. "
        "Avoids false HOM calls when a small number of ref reads are present.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--genotype-error-rate",
        help="Alternate-allele fraction expected at a homozygous-reference locus, "
        "and (as 1 - e) the fraction expected at a homozygous-alternate one "
        f"(default: {DEFAULT_GENOTYPE_ERROR_RATE}). This is not only a read error "
        "rate: because the genotype is the most likely of a set of binomial "
        "hypotheses evaluated at fixed expected fractions, this one number also "
        "fixes the allele-fraction boundaries between genotypes. At the default "
        "0.05 a locus is called 0/0 below AF 0.218, 0/1 between 0.218 and 0.782, "
        "and 1/1 above 0.782, at every depth; raising it moves both boundaries "
        "inward, trading homozygote confidence against heterozygote confidence. "
        "Genotypes from runs with different values are not comparable.",
        type=float,
        default=DEFAULT_GENOTYPE_ERROR_RATE,
    )
    parser.add_argument(
        "--legacy-force-wildtype-genotypes",
        help="DEPRECATED CONTROL, kept for one release to reproduce pre-v0.3 "
        "genotype fields exactly, for a controlled re-benchmark. Despite the name "
        "it now restores every pre-fix genotype behaviour, not only the "
        "force-called wild types: (1) force-calls 0/0 with GQ=60 and GP=1.0 "
        "whenever no alternate-supporting read is found, including at loci with no "
        "read coverage at all (without the flag such loci are no-calls, ./. with "
        "TC=0, GQ=0, GP=., and covered reference loci get a depth-derived GQ); "
        "(2) restores the pre-fix coverage-query geometry (no breakpoint margin "
        "for insertions, 100 bases for deletions), overriding "
        "--genotype-breakpoint-margin; (3) restores the pre-fix uncapped, "
        "non-monotone GQ on alternate-supported calls (which could exceed the 60 "
        "the VCF header declares, and dropped back to 60 when the posterior "
        "saturated). Use only to reproduce a pre-fix run for comparison; the "
        "output is not a correct call set.",
        action="store_true",
        default=False,
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
