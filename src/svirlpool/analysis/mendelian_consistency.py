# %%
# script to test Mendelian consistency in a vcf file
#
# This module implements a generalized Mendelian consistency checker that works
# across different copy numbers (CN), including:
#
# - CN=0: Homozygous deletion (marked as uncertain)
# - CN=1: Hemizygous (e.g., male X chromosome): genotypes like "0" or "1"
# - CN=2: Diploid (standard case): genotypes like "0/0", "0/1", "1/1"
# - CN=3: Triploid (duplication): genotypes like "0/0/1", "0/1/1", "1/1/1"
# - CN=4: Tetraploid: genotypes like "0/0/0/1", "0/0/1/1", "0/1/1/1", "1/1/1/1"
#
# Core Algorithm:
# ---------------
# The checker uses an allele-counting approach based on the fundamental rule
# that each allele in the child must be inherited from one of the parents.
#
# 1. Parse genotypes to count ref (0) and alt (non-0) alleles
# 2. Check if parents can collectively provide the child's alleles
#
# For diploid (CN=2):
#   - Child 0/0: Always possible from any parents
#   - Child 0/1: At least one parent must have an alt allele
#   - Child 1/1: Both parents must have at least one alt allele each
#
# For hemizygous (CN=1):
#   - Child inherits from mother only (e.g., X chromosome in males)
#   - Child 0: Mother can be any genotype
#   - Child 1: Mother must have at least one alt (0/1, 1/1, or 1)
#
# For higher CN (3-4):
#   - Conservative check: parents must collectively have enough alt alleles
#   - CN mismatches are flagged as "uncertain" for manual review
#
# Status codes:
# -------------
# - consistent: Follows Mendelian inheritance rules
# - inconsistent: Violates Mendelian rules (e.g., child has alts neither parent has)
# - incomplete: One parent missing, but could be consistent
# - missing: Child genotype missing or both parents missing
# - uncertain: Complex CN patterns requiring manual review
# - non_informative: every trio member was observed and none carries an alt
# - no_call: at least one trio member was not observed at all (./.)
#
# no_call vs. non_informative
# ---------------------------
# These are different claims and are kept apart deliberately.  "Nobody in the
# trio carries this variant" is an observation; "somebody in the trio was not
# observed" is the absence of one.  Folding the second into the first -- which
# is what treating `./.` as `0/0` does -- scores a technical dropout as a
# correctly inherited genotype whenever it lands in the child, and as a
# Mendelian violation whenever it lands in a parent of a carrier child.
#
# Consequently:
#
# - a trio containing any no-call is reported as `no_call` and is excluded from
#   the consistent/inconsistent denominator, so the headline consistency figure
#   means "among trios in which all three genotypes were observed, how many obey
#   Mendelian inheritance";
# - the no-call rate is reported alongside it (`no_call` / all evaluated trios)
#   so the headline cannot be improved by calling nothing;
# - the pre-v0.3 reading, in which `./.` was counted as an observed homozygous
#   reference, stays reachable via `--legacy-missing-as-reference` so that
#   previously published numbers can be reproduced exactly.

import argparse
import json
import logging
import subprocess
import tempfile
from collections.abc import Iterable
from enum import Enum
from pathlib import Path, PosixPath

import vcfpy
from intervaltree import IntervalTree

# logger
log = logging.getLogger(__name__)

# %%


# parse chr,start and GTs
def parse_variant(variant):
    return (variant.CHROM, variant.POS), [call.data["GT"] for call in variant.calls]


# %%
# Allele separators recognised in a GT field.  "|" marks a phased genotype and
# is as much a separator as "/"; reading "0|1" as a single allele made every
# phased heterozygote look like a hemizygous alt.
_ALLELE_SEPARATORS = ("/", "|")


def _split_alleles(gt_string: str) -> list[str]:
    """Split a GT field into its alleles on either separator."""
    normalised = gt_string.replace("|", "/")
    return normalised.split("/")


def parse_genotype(
    gt_string: str | None, missing_as_reference: bool = False
) -> tuple[int, int, int]:
    """
    Parse a genotype string into observed allele counts.

    Returns: (n_ref_alleles, n_alt_alleles, ploidy)

    For an observed genotype ``n_ref + n_alt == ploidy``.  A **no-call** is
    reported as ``(0, 0, ploidy)``: no reference alleles were observed and no
    alternate alleles were observed either.  That state is unreachable for an
    observed genotype, so :func:`is_no_call` can recognise it without a fourth
    tuple element, and a no-call cannot be mistaken for a homozygous reference
    call the way ``(ploidy, 0, ploidy)`` was.

    A genotype in which *any* allele is missing (``0/.``) is a no-call: half an
    observation cannot support a Mendelian determination.  ``None`` -- which is
    what vcfpy decodes a bare ``.`` GT field into -- and the empty string are
    no-calls too.

    Args:
        gt_string: the raw GT field, or None when the field is absent.
        missing_as_reference: reproduce the pre-v0.3 reading, in which a missing
            genotype was counted as an observed homozygous reference.  Only for
            reproducing previously published numbers; see
            ``--legacy-missing-as-reference``.

    Examples:
        "0" → (1, 0, 1)         # hemizygous ref
        "1" → (0, 1, 1)         # hemizygous alt
        "0/0" → (2, 0, 2)       # homozygous ref
        "0/1" → (1, 1, 2)       # heterozygous
        "0|1" → (1, 1, 2)       # phased heterozygous
        "1/1" → (0, 2, 2)       # homozygous alt
        "0/0/0" → (3, 0, 3)     # triploid ref (copy number 3)
        "0/0/1" → (2, 1, 3)     # triploid with 1 alt
        "0/1/1/1" → (1, 3, 4)   # tetraploid with 3 alts
        "1/2" → (0, 2, 2)       # any non-zero allele index counts as alt
        "./." → (0, 0, 2)       # no-call, diploid
        "." / None → (0, 0, 1)  # no-call, hemizygous
        "0/." → (0, 0, 2)       # partial no-call
    """
    if gt_string is None or gt_string == "":
        # vcfpy decodes an absent GT field as None.  Nothing was observed.
        return (0, 0, 1) if missing_as_reference is False else (1, 0, 1)

    if "." in gt_string:
        n_alleles = 1
        for separator in _ALLELE_SEPARATORS:
            if separator in gt_string:
                n_alleles = len(_split_alleles(gt_string))
                break
        if missing_as_reference:
            # Pre-v0.3 behaviour: a missing genotype was counted as an observed
            # homozygous reference call.
            return (n_alleles, 0, n_alleles)
        # No information: neither reference nor alternate alleles were observed.
        return (0, 0, n_alleles)

    alleles = _split_alleles(gt_string)

    # Count reference (0) and alt (non-0) alleles
    n_ref = alleles.count("0")
    # For now, treat any non-0 allele as alt (handles 1, 2, etc.)
    n_alt = sum(1 for a in alleles if a != "0")
    total = len(alleles)

    return (n_ref, n_alt, total)


def is_no_call(parsed_genotype: tuple[int, int, int]) -> bool:
    """True when :func:`parse_genotype` reported "nothing was observed".

    An observed genotype always accounts for every copy -- ``n_ref + n_alt ==
    ploidy`` -- so zero of both is only ever produced by a missing genotype.
    """
    n_ref, n_alt, _ = parsed_genotype
    return n_ref == 0 and n_alt == 0


# %%
def check_diploid_inheritance(
    child_n_alt: int, father: tuple[int, int, int], mother: tuple[int, int, int]
) -> str:
    """
    Check Mendelian inheritance for diploid child (CN=2).
    Child receives 1 allele from father, 1 from mother.

    Rules:
        - child_n_alt = 0 (0/0): Both parents can be any genotype
        - child_n_alt = 1 (0/1): At least one parent must have alt (0/1, 1/1, or hemizygous 1)
        - child_n_alt = 2 (1/1): Both parents must have at least one alt allele each

    Returns: "consistent", "inconsistent", or "incomplete"
    """
    # Both parents present
    _, father_n_alt, _ = father
    _, mother_n_alt, _ = mother

    if child_n_alt == 0:
        # Both alleles are ref - always possible
        return "consistent"

    elif child_n_alt == 1:
        # One alt allele - at least one parent needs an alt
        if father_n_alt >= 1 or mother_n_alt >= 1:
            return "consistent"
        else:
            return "inconsistent"  # 0/0 × 0/0 cannot produce 0/1

    elif child_n_alt == 2:
        # Both alleles are alt - both parents need at least one alt
        if father_n_alt >= 1 and mother_n_alt >= 1:
            return "consistent"
        else:
            return "inconsistent"  # Need alt from both parents

    return "inconsistent"


# %%
def check_hemizygous_inheritance(child_n_alt: int, mother: tuple[int, int, int]) -> str:
    """
    Check Mendelian inheritance for hemizygous child (CN=1).
    For hemizygous regions (e.g., male X chromosome), child inherits from mother only.

    Rules:
        - child = 0: Mother can be 0/0 or 0/1 (or 0)
        - child = 1: Mother must have 0/1 or 1/1 (or 1) - at least one alt

    Returns: "consistent", "inconsistent", or "incomplete"
    """
    _, mother_n_alt, _ = mother

    if child_n_alt == 0:
        # Ref allele - mother can have any genotype with at least one ref
        return "consistent"

    elif child_n_alt == 1:
        # Alt allele - mother must have at least one alt
        if mother_n_alt >= 1:
            return "consistent"
        else:
            return "inconsistent"  # Mother 0/0 cannot give alt to child

    return "inconsistent"


# %%
def check_higher_cn_inheritance(
    child_n_alt: int,
    child_cn: int,
    father: tuple[int, int, int],
    mother: tuple[int, int, int],
) -> str:
    """
    Check Mendelian inheritance for higher copy number children (CN >= 3).

    This is complex because we don't know which parent(s) have the duplication.
    We use a conservative approach: check if it's POSSIBLE given parental alleles.

    Rule: Count total available alt alleles from parents.
    If parents collectively can provide >= child_n_alt, then it's possible.

    Note: This may allow some false positives when CN structure is complex,
    but avoids false negatives. Cases flagged as "uncertain" need manual review.

    Returns: "consistent", "inconsistent", "incomplete", or "uncertain"
    """
    # Both parents present
    _, father_n_alt, father_cn = father
    _, mother_n_alt, mother_cn = mother

    # Check for CN mismatch - suggests de novo or complex inheritance
    if child_cn not in [father_cn, mother_cn, 2]:
        # Child has different CN than both parents (and not standard diploid)
        # This suggests de novo duplication or complex scenario
        log.debug(
            f"CN mismatch: child={child_cn}, father={father_cn}, mother={mother_cn}"
        )

    # Total alt alleles available from parents
    total_parental_alts = father_n_alt + mother_n_alt

    if child_n_alt == 0:
        # All ref - always possible
        return "consistent"

    # Minimum check: Can parents provide enough alt alleles?
    if total_parental_alts >= child_n_alt:
        # Potentially consistent
        return "consistent"
    else:
        # Definitely inconsistent - not enough alt alleles available
        return "inconsistent"


# define an enum for gt inheritance status: consistent, inconsistent, missing
class GTInheritanceStatus(Enum):
    missing = 0
    consistent = 1
    incomplete = 2
    inconsistent = 3
    uncertain = 4  # For complex CN cases that need manual review
    non_informative = 5  # Every trio member observed, none carries an alt
    no_call = 6  # At least one trio member was not observed at all (./.)


# %%
def is_trio_informative(
    child: tuple[int, int, int],
    father: tuple[int, int, int],
    mother: tuple[int, int, int],
) -> bool:
    """
    Check if a trio has at least one member with an *observed* alternate allele.

    Returns False if no member carries an alt allele (n_alt = 0 for all three).

    A no-call carries no observed alt allele either, so it cannot make a trio
    informative.  Callers must therefore test for a no-call *before* calling
    this, or a dropout will be filed as "nobody carries the variant"; see
    :func:`is_variant_inconsistent`.

    Args:
        child: Parsed child genotype (n_ref, n_alt, cn)
        father: Parsed father genotype (n_ref, n_alt, cn)
        mother: Parsed mother genotype (n_ref, n_alt, cn)

    Returns:
        True if at least one member has n_alt > 0, False otherwise
    """
    # Check if anyone has an alternate allele
    for person in [child, father, mother]:
        _, n_alt, _ = person
        if n_alt > 0:
            return True  # At least one person has a variant

    # All members are homozygous reference (including treated missing genotypes)
    return False


# %%
def is_variant_inconsistent(
    variant: vcfpy.Record,
    names_trio: list[str] | None = None,
    missing_as_reference: bool = False,
) -> GTInheritanceStatus:
    """
    Check Mendelian consistency for a variant across all copy numbers (CN 0-4+).

    Uses generalized allele counting approach that works for:
    - Hemizygous (CN=1): e.g., "0" or "1"
    - Diploid (CN=2): e.g., "0/0", "0/1", "1/1"
    - Triploid (CN=3): e.g., "0/0/1", "0/1/1", "1/1/1"
    - Tetraploid (CN=4): e.g., "0/0/0/1", "0/0/1/1", etc.

    A trio in which any member was not observed is reported as
    ``GTInheritanceStatus.no_call``: the inheritance of an allele that was never
    seen is not determinable, and calling it either way is a claim the data does
    not support.  The check precedes the informativeness test, because a
    dropout is not the same observation as "nobody in this trio carries it".

    Args:
        variant: VCF record
        names_trio: List of [child, father, mother] sample names
        missing_as_reference: reproduce the pre-v0.3 reading in which ``./.``
            counted as an observed homozygous reference, so that previously
            published numbers can be reproduced.  ``no_call`` is never returned
            under this flag.

    Returns:
        GTInheritanceStatus enum value
    """
    # names_trio has 3 elements, the names of son, father, mother
    if names_trio is not None:
        assert len(names_trio) == 3, "names_trio must have 3 elements"
        gt_son = None
        gt_father = None
        gt_mother = None
        all_sample_names = {call.sample for call in variant.calls}
        assert all_sample_names.intersection(set(names_trio)) == set(names_trio), (
            "names_trio must be in the vcf file"
        )
        for call in variant.calls:
            if call.sample == names_trio[0]:
                gt_son = call.data["GT"]
            elif call.sample == names_trio[1]:
                gt_father = call.data["GT"]
            elif call.sample == names_trio[2]:
                gt_mother = call.data["GT"]
    else:
        gt_son = variant.calls[0].data["GT"]
        gt_father = variant.calls[1].data["GT"]
        gt_mother = variant.calls[2].data["GT"]

    # Parse genotypes to allele counts
    child = parse_genotype(gt_son, missing_as_reference=missing_as_reference)
    father = parse_genotype(gt_father, missing_as_reference=missing_as_reference)
    mother = parse_genotype(gt_mother, missing_as_reference=missing_as_reference)

    # A member that was not observed makes the trio undeterminable.  This is
    # checked before informativeness: a no-call carries no alt allele, so it
    # would otherwise be filed as "nobody in the trio carries the variant".
    if any(is_no_call(person) for person in (child, father, mother)):
        return GTInheritanceStatus.no_call

    # Check if trio is informative (at least one member has variant)
    if not is_trio_informative(child, father, mother):
        return GTInheritanceStatus.non_informative

    _, child_n_alt, child_cn = child

    # Dispatch based on child's copy number
    if child_cn == 1:
        # Hemizygous inheritance (e.g., male X chromosome)
        status = check_hemizygous_inheritance(child_n_alt, mother)
    elif child_cn == 2:
        # Standard diploid inheritance
        status = check_diploid_inheritance(child_n_alt, father, mother)
    elif child_cn in [3, 4]:
        # Higher copy number (duplications/triplications)
        status = check_higher_cn_inheritance(child_n_alt, child_cn, father, mother)
    elif child_cn == 0:
        # Homozygous deletion - no copies
        # This is consistent if parents can both contribute deletions
        # For now, mark as uncertain - needs special handling
        log.debug(f"Child has CN=0 at {variant.CHROM}:{variant.POS}")
        return GTInheritanceStatus.uncertain
    else:
        # CN > 4 or unexpected
        log.warning(f"Unexpected CN={child_cn} at {variant.CHROM}:{variant.POS}")
        return GTInheritanceStatus.uncertain

    # Map string status to enum
    if status == "consistent":
        return GTInheritanceStatus.consistent
    elif status == "inconsistent":
        return GTInheritanceStatus.inconsistent
    elif status == "incomplete":
        return GTInheritanceStatus.incomplete
    elif status == "uncertain":
        return GTInheritanceStatus.uncertain
    else:
        return GTInheritanceStatus.inconsistent


# %%


def load_regions_from_bed(bed_path: Path) -> dict[str, IntervalTree]:
    """Load a BED file into a dict of IntervalTrees keyed by chromosome."""
    trees: dict[str, IntervalTree] = {}
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if start == end:
                end += 1  # IntervalTree requires non-zero-length intervals
            if chrom not in trees:
                trees[chrom] = IntervalTree()
            trees[chrom].addi(start, end)
    return trees


def variant_overlaps_regions(
    variant: vcfpy.Record, regions: dict[str, IntervalTree]
) -> bool:
    """Return True if the variant overlaps at least one region in *regions*."""
    chrom = variant.CHROM
    if chrom not in regions:
        return False
    start = variant.POS
    svtype = variant.INFO.get("SVTYPE", "")
    if svtype == "DEL" and "END" in variant.INFO:
        end = int(variant.INFO["END"])
    else:
        end = start + 1
    if start == end:
        end += 1
    return bool(regions[chrom].overlap(start, end))


def check_ped(path_ped: PosixPath):
    with open(path_ped, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            assert len(line) == 6, (
                "ped file must have 6 columns: familyID, individualID, fatherID, motherID, sex, phenotype"
            )


def parse_ped(path_ped: PosixPath) -> dict[str, tuple[str, str, str]]:
    check_ped(path_ped)
    "returns a dictionary with the samplename as key and a tuple with father, mother, family as value"
    dict_ped = dict[str, tuple[str, str, str]]()
    with open(path_ped, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            dict_ped[line[1]] = (line[2], line[3], line[0])
    return dict_ped


def parse_variants(path_vcf: PosixPath, passonly: bool = False) -> list[vcfpy.Record]:
    reader = vcfpy.Reader.from_path(path_vcf)
    if passonly:
        return [record for record in reader if record.FILTER == ["PASS"]]
    return list(reader)


def extract_variant_features(variant: vcfpy.Record) -> str:
    """
    Extract features from a variant record for BED file output.

    Returns a comma-separated string of features: SVLEN,SVTYPE,PRECISION
    where PRECISION is either "IMPRECISE" or "PRECISE" based on the INFO field.

    Args:
        variant: VCF record

    Returns:
        Comma-separated string of features
    """
    # Extract SVLEN
    svlen = variant.INFO.get("SVLEN", "NA")
    if isinstance(svlen, list):
        svlen = svlen[0]

    # Extract SVTYPE
    svtype = variant.INFO.get("SVTYPE", "NA")

    # Check IMPRECISE flag
    if "IMPRECISE" in variant.INFO:
        precision = "IMPRECISE"
    else:
        precision = "PRECISE"

    return f"{svlen},{svtype},{precision}"


# Statuses that are not Mendelian verdicts and are therefore kept out of the
# consistent/inconsistent denominator.  `no_call` joins the two pre-existing
# members: the headline consistency figure means "among trios in which all
# three genotypes were observed and at least one carries the variant, how many
# obey Mendelian inheritance".
_NOT_A_VERDICT: frozenset[GTInheritanceStatus] = frozenset({
    GTInheritanceStatus.uncertain,
    GTInheritanceStatus.non_informative,
    GTInheritanceStatus.no_call,
})


def _empty_status_dict() -> dict[GTInheritanceStatus, int]:
    # Explicit order: the historical row order of every emitted table, with
    # no_call appended so existing outputs keep their shape.
    return {
        GTInheritanceStatus.consistent: 0,
        GTInheritanceStatus.inconsistent: 0,
        GTInheritanceStatus.missing: 0,
        GTInheritanceStatus.incomplete: 0,
        GTInheritanceStatus.uncertain: 0,
        GTInheritanceStatus.non_informative: 0,
        GTInheritanceStatus.no_call: 0,
    }


def n_informative_trios(dict_stats: dict[GTInheritanceStatus, int]) -> int:
    """Number of trios that received an actual Mendelian verdict."""
    return sum(
        count for status, count in dict_stats.items() if status not in _NOT_A_VERDICT
    )


def no_call_rate(dict_stats: dict[GTInheritanceStatus, int]) -> float:
    """Fraction of *all* evaluated trios in which somebody was not observed.

    Reported next to the consistency figure on purpose.  Excluding no-calls
    from the consistency denominator is only honest if the size of that
    exclusion is visible: otherwise the metric can be improved by calling less.
    """
    n_total = sum(dict_stats.values())
    if n_total == 0:
        return 0.0
    return dict_stats.get(GTInheritanceStatus.no_call, 0) / n_total


# Size stratification bins: (low_inclusive, high_exclusive, label)
SIZE_BINS: list[tuple[int, int | None, str]] = [
    (0, 100, "<100"),
    (100, 350, "100-350"),
    (350, 1000, "350-1000"),
    (1000, 10000, "1k-10k"),
    (10000, None, ">10k"),
]
SIZE_BIN_LABELS: list[str] = [label for _, _, label in SIZE_BINS]


def get_variant_size(variant: vcfpy.Record) -> int | None:
    """Return the size (in bp) of a structural variant, or None for BNDs.

    Priority:
    1. Absolute value of SVLEN INFO field.
    2. Sequence-resolved fallback:
       - DEL: len(REF) - 1 (anchor base subtracted)
       - INS: len(ALT[0]) - 1
       - INV: max(len(REF), len(ALT[0])) - 1
       - other: max(len(REF), len(ALT[0])) - 1
    """
    svtype = variant.INFO.get("SVTYPE", "")
    if svtype == "BND":
        return None
    svlen = variant.INFO.get("SVLEN", None)
    if svlen is not None:
        if isinstance(svlen, list):
            svlen = svlen[0]
        return abs(int(svlen))
    ref = variant.REF
    alt = variant.ALT[0].value if variant.ALT else ""
    if svtype == "DEL":
        size = len(ref) - 1
    elif svtype == "INS":
        size = len(alt) - 1
    elif svtype == "INV":
        size = max(len(ref), len(alt)) - 1
    else:
        size = max(len(ref), len(alt)) - 1
    return size if size > 0 else None


def get_size_bin(size: int) -> str:
    """Return the label for the bin that *size* (bp) falls into."""
    for low, high, label in SIZE_BINS:
        if size >= low and (high is None or size < high):
            return label
    return SIZE_BIN_LABELS[-1]


# creates the output and results line of one trio
def variants_consistency_stats(
    variants: list[vcfpy.Record],
    output: PosixPath,
    names_trio: list[str] | None = None,
    bed: PosixPath | None = None,
    min_size: int = 0,
    regions: dict[str, IntervalTree] | None = None,
    svtypes: set[str] | None = None,
    size_stratification: bool = True,
    missing_as_reference: bool = False,
) -> tuple[
    dict[str, dict[GTInheritanceStatus, int]],
    dict[str, dict[str, dict[GTInheritanceStatus, int]]] | None,
]:
    """
    Calculate Mendelian consistency statistics for variants in a trio,
    stratified by SVTYPE.

    Args:
        variants: List of VCF records
        output: Output BED file path
        names_trio: Optional list of [child, father, mother] sample names
        bed: Optional path for BED file with features (chr, start, end, features, status)
        min_size: Minimum absolute SVLEN to include
        regions: Optional IntervalTree dict for region filtering
        svtypes: Optional set of SVTYPEs to include. If None, all SVTYPEs are included.
        size_stratification: If True, also return counts stratified by size bin (BNDs excluded).
        missing_as_reference: reproduce the pre-v0.3 reading of ``./.``; see
            ``--legacy-missing-as-reference``.

    Returns:
        Tuple of:
          - Nested dict: svtype ("all", <type1>, ...) → GTInheritanceStatus → count
          - Size-stratified dict: svtype → size_bin → GTInheritanceStatus → count,
            or None when size_stratification is False.
    """
    if names_trio is not None:
        assert len(names_trio) == 3, "names_trio must have 3 elements"

    bed_file = open(bed, "w") if bed else None

    dict_stats: dict[str, dict[GTInheritanceStatus, int]] = {
        "all": _empty_status_dict(),
    }
    size_stats: dict[str, dict[str, dict[GTInheritanceStatus, int]]] | None = (
        {"all": {label: _empty_status_dict() for label in SIZE_BIN_LABELS}}
        if size_stratification
        else None
    )

    with open(output, "w"):
        for variant in variants:
            svtype = variant.INFO.get("SVTYPE", "")
            if svtypes is not None and svtype not in svtypes:
                continue
            if regions is not None and not variant_overlaps_regions(variant, regions):
                continue
            if min_size > 0:
                svlen = variant.INFO.get("SVLEN", [0])
                if isinstance(svlen, list):
                    svlen = svlen[0]
                if abs(int(svlen)) < min_size:
                    continue
            status = is_variant_inconsistent(
                variant=variant,
                names_trio=names_trio,
                missing_as_reference=missing_as_reference,
            )
            dict_stats["all"][status] += 1
            if svtype not in dict_stats:
                dict_stats[svtype] = _empty_status_dict()
            dict_stats[svtype][status] += 1

            if size_stats is not None and svtype != "BND":
                size = get_variant_size(variant)
                if size is not None:
                    bin_label = get_size_bin(size)
                    size_stats["all"][bin_label][status] += 1
                    if svtype not in size_stats:
                        size_stats[svtype] = {
                            label: _empty_status_dict() for label in SIZE_BIN_LABELS
                        }
                    size_stats[svtype][bin_label][status] += 1

            end = int(variant.INFO["END"]) if svtype == "DEL" else variant.POS + 1

            if bed_file:
                features = extract_variant_features(variant)
                print(
                    f"{variant.CHROM}\t{variant.POS}\t{end}\t{status.name},{features}",
                    file=bed_file,
                )

    if bed_file:
        bed_file.close()

    return dict_stats, size_stats


def _status_dict_to_lines(
    dict_stats: dict[GTInheritanceStatus, int], indent: str = "  "
) -> str:
    """Format a single GTInheritanceStatus → count dict as indented lines."""
    n_total = sum(dict_stats.values())
    n_informative = n_informative_trios(dict_stats)
    lines = f"{indent}- total: {n_total:,}\n"
    for status in dict_stats:
        count = dict_stats[status]
        if status is GTInheritanceStatus.no_call:
            # Percentage of *all* evaluated trios, not of the verdict
            # denominator, which no-calls are deliberately outside of.
            pct = count / n_total if n_total > 0 else 0
            lines += f"{indent}- {status.name}: {count:,}; {pct:.2%} of total\n"
        elif status in _NOT_A_VERDICT:
            lines += f"{indent}- {status.name}: {count:,}\n"
        else:
            pct = count / n_informative if n_informative > 0 else 0
            lines += f"{indent}- {status.name}: {count:,}; {pct:.2%}\n"
    lines += f"{indent}- informative (consistency denominator): {n_informative:,}\n"
    lines += f"{indent}- no_call_rate: {no_call_rate(dict_stats):.2%}\n"
    return lines


def dict_stats_to_line(
    samplename: str,
    dict_stats: dict[str, dict[GTInheritanceStatus, int]],
    size_stats: dict[str, dict[str, dict[GTInheritanceStatus, int]]] | None = None,
) -> str:
    """Format per-sample stratified stats as a human-readable block."""
    lines = f"{samplename}:\n"
    svtype_keys = ["all"] + sorted(k for k in dict_stats if k != "all")
    for svtype in svtype_keys:
        lines += f"  [{svtype}]\n"
        lines += _status_dict_to_lines(dict_stats[svtype], indent="  ")
        if size_stats is not None and svtype in size_stats:
            for bin_label in SIZE_BIN_LABELS:
                lines += f"  [{svtype} / {bin_label}]\n"
                lines += _status_dict_to_lines(
                    size_stats[svtype][bin_label], indent="    "
                )
    return lines.rstrip()


def get_trios_from_variants(
    variants: list[vcfpy.Record], ped: PosixPath
) -> list[list[str]]:
    "returns a list of trios that are present in the vcf file and in the pedigree. Each sublist has 3 elements: child, father, mother"
    all_samplenames = {call.sample for variant in variants for call in variant.calls}
    # then parse the pedigree and find all trios that are alos present in all_samplenames
    dict_ped: dict[str, tuple[str, str, str]] = parse_ped(ped)
    trios = []
    for samplename in all_samplenames:
        if samplename in dict_ped:
            father_name = dict_ped[samplename][0]
            mother_name = dict_ped[samplename][1]
            if father_name in all_samplenames and mother_name in all_samplenames:
                trios.append([samplename, father_name, mother_name])
            else:
                log.warning(f"trio of {samplename} is not present in the vcf file")
                if father_name not in all_samplenames:
                    log.warning(f"father {father_name} is not present in the vcf file")
                if mother_name not in all_samplenames:
                    log.warning(f"mother {mother_name} is not present in the vcf file")
        else:
            log.warning(f"{samplename} is not present in the pedigree")
    return trios


def all_dict_stats_to_tsv(
    all_dict_stats: dict[str, dict[str, dict[GTInheritanceStatus, int]]],
    path_tsv: PosixPath,
    all_size_dict_stats: dict[str, dict[str, dict[str, dict[GTInheritanceStatus, int]]]]
    | None = None,
) -> None:
    """Write stratified stats to a TSV.

    Columns: sample, svtype, [size_bin,] status, count, percentage, denominator.

    ``denominator`` names what ``percentage`` is a fraction of, because the two
    are not the same for every row:

    - ``informative`` -- consistent / inconsistent / incomplete / missing, as a
      fraction of the trios that received a Mendelian verdict;
    - ``total`` -- ``no_call``, as a fraction of *all* evaluated trios, i.e. the
      no-call rate;
    - ``NA`` -- uncertain / non_informative, which are not rates.
    """

    def _percentage_and_denominator(
        status: GTInheritanceStatus, count: int, n_informative: int, n_total: int
    ) -> tuple[str, str]:
        if status is GTInheritanceStatus.no_call:
            return (str(count / n_total if n_total > 0 else 0), "total")
        if status in _NOT_A_VERDICT:
            return ("NA", "NA")
        return (str(count / n_informative if n_informative > 0 else 0), "informative")

    with open(path_tsv, "w") as f:
        if all_size_dict_stats is not None:
            print(
                "sample\tsvtype\tsize_bin\tstatus\tcount\tpercentage\tdenominator",
                file=f,
            )
        else:
            print("sample\tsvtype\tstatus\tcount\tpercentage\tdenominator", file=f)
        for sample, svtype_stats in all_dict_stats.items():
            svtype_keys = ["all"] + sorted(k for k in svtype_stats if k != "all")
            for svtype in svtype_keys:
                dict_stats = svtype_stats[svtype]
                n_informative = n_informative_trios(dict_stats)
                n_total = sum(dict_stats.values())
                for status in dict_stats:
                    count = dict_stats[status]
                    percentage_str, denominator = _percentage_and_denominator(
                        status, count, n_informative, n_total
                    )
                    if all_size_dict_stats is not None:
                        print(
                            f"{sample}\t{svtype}\tall\t{status.name}\t{count}"
                            f"\t{percentage_str}\t{denominator}",
                            file=f,
                        )
                    else:
                        print(
                            f"{sample}\t{svtype}\t{status.name}\t{count}"
                            f"\t{percentage_str}\t{denominator}",
                            file=f,
                        )
                # Size-stratified rows
                if all_size_dict_stats is not None and sample in all_size_dict_stats:
                    size_stats_sample = all_size_dict_stats[sample]
                    if svtype in size_stats_sample:
                        for bin_label in SIZE_BIN_LABELS:
                            bin_stats = size_stats_sample[svtype][bin_label]
                            n_inf = n_informative_trios(bin_stats)
                            n_bin_total = sum(bin_stats.values())
                            for status in bin_stats:
                                count = bin_stats[status]
                                pct_str, denominator = _percentage_and_denominator(
                                    status, count, n_inf, n_bin_total
                                )
                                print(
                                    f"{sample}\t{svtype}\t{bin_label}\t{status.name}"
                                    f"\t{count}\t{pct_str}\t{denominator}",
                                    file=f,
                                )


def all_dict_stats_to_latex(
    all_dict_stats: dict[str, dict[str, dict[GTInheritanceStatus, int]]],
    path_latex: PosixPath,
    all_size_dict_stats: dict[str, dict[str, dict[str, dict[GTInheritanceStatus, int]]]]
    | None = None,
) -> None:
    """Write stratified stats to a LaTeX table.

    Columns: status, svtype, [size_bin,] <sample1>, <sample2>, ...
    One row per (status, svtype[, size_bin]) combination showing counts and fractional percentages.
    """
    samples = list(all_dict_stats.keys())
    # Collect all svtype keys across all samples
    all_svtype_keys_set: set[str] = set()
    for svtype_stats in all_dict_stats.values():
        all_svtype_keys_set.update(svtype_stats.keys())
    svtype_keys = ["all"] + sorted(k for k in all_svtype_keys_set if k != "all")
    with open(path_latex, "w") as f:
        if all_size_dict_stats is not None:
            header = r"status & svtype & size\_bin & " + " & ".join(samples) + r"\\"
        else:
            header = r"status & svtype & " + " & ".join(samples) + r"\\"
        print(header, file=f)
        for status in GTInheritanceStatus:
            for svtype in svtype_keys:
                counts = [
                    all_dict_stats[s][svtype][status]
                    if svtype in all_dict_stats[s]
                    else 0
                    for s in samples
                ]
                totals = [
                    sum(all_dict_stats[s][svtype].values())
                    if svtype in all_dict_stats[s]
                    else 0
                    for s in samples
                ]
                abs_cells = " & ".join(f"{c:,}" for c in counts)
                pct_cells = " & ".join(
                    f"{c / t:.2f}" if t > 0 else "NA"
                    for c, t in zip(counts, totals, strict=False)
                )
                if all_size_dict_stats is not None:
                    print(f"{status.name} & {svtype} & all & {abs_cells} \\\\", file=f)
                    print(f"  & & & {pct_cells} \\\\", file=f)
                else:
                    print(f"{status.name} ({svtype}) & {abs_cells} \\\\", file=f)
                    print(f"  & & {pct_cells} \\\\", file=f)
                # Size-stratified rows
                if all_size_dict_stats is not None:
                    for bin_label in SIZE_BIN_LABELS:
                        bin_counts = [
                            (
                                all_size_dict_stats[s][svtype][bin_label][status]
                                if s in all_size_dict_stats
                                and svtype in all_size_dict_stats[s]
                                else 0
                            )
                            for s in samples
                        ]
                        bin_totals = [
                            (
                                sum(all_size_dict_stats[s][svtype][bin_label].values())
                                if s in all_size_dict_stats
                                and svtype in all_size_dict_stats[s]
                                else 0
                            )
                            for s in samples
                        ]
                        bc_cells = " & ".join(f"{c:,}" for c in bin_counts)
                        bp_cells = " & ".join(
                            f"{c / t:.2f}" if t > 0 else "NA"
                            for c, t in zip(bin_counts, bin_totals, strict=False)
                        )
                        print(
                            f"{status.name} & {svtype} & {bin_label} & {bc_cells} \\\\",
                            file=f,
                        )
                        print(f"  & & & {bp_cells} \\\\", file=f)


def _status_dict_to_json_entry(d: dict[GTInheritanceStatus, int]) -> dict:
    """One JSON object per stratum: per-status counts plus the no-call rate.

    Each status carries the ``denominator`` its ``percentage`` is taken over,
    so a reader never has to guess whether a number is a share of the verdicts
    or of every evaluated trio.
    """
    n_informative = n_informative_trios(d)
    n_total = sum(d.values())
    entry: dict = {}
    for status, count in d.items():
        if status is GTInheritanceStatus.no_call:
            entry[status.name] = {
                "count": count,
                "percentage": count / n_total if n_total > 0 else 0,
                "denominator": "total",
            }
        elif status in _NOT_A_VERDICT:
            entry[status.name] = {
                "count": count,
                "percentage": None,
                "denominator": None,
            }
        else:
            entry[status.name] = {
                "count": count,
                "percentage": count / n_informative if n_informative > 0 else 0,
                "denominator": "informative",
            }
    # Grouped under "summary" so the rest of the object stays a clean
    # status -> {count, percentage, denominator} map.
    entry["summary"] = {
        "total": n_total,
        "informative": n_informative,
        "no_call_rate": no_call_rate(d),
    }
    return entry


def all_dict_stats_to_json(
    all_dict_stats: dict[str, dict[str, dict[GTInheritanceStatus, int]]],
    path_json: PosixPath,
    all_size_dict_stats: dict[str, dict[str, dict[str, dict[GTInheritanceStatus, int]]]]
    | None = None,
) -> None:
    """Write stratified stats to a JSON file."""
    out: dict = {}
    for sample, svtype_stats in all_dict_stats.items():
        sample_obj: dict = {}
        svtype_keys = _sort_svtype_keys(svtype_stats.keys())
        for svtype in svtype_keys:
            d = svtype_stats[svtype]
            entry = _status_dict_to_json_entry(d)
            sample_obj[svtype] = entry
            # Size-stratified sub-entries
            if (
                all_size_dict_stats is not None
                and sample in all_size_dict_stats
                and svtype in all_size_dict_stats[sample]
            ):
                size_entry: dict = {}
                for bin_label in SIZE_BIN_LABELS:
                    bin_d = all_size_dict_stats[sample][svtype][bin_label]
                    size_entry[bin_label] = _status_dict_to_json_entry(bin_d)
                sample_obj[f"{svtype}_by_size"] = size_entry
        out[sample] = sample_obj
    with open(path_json, "w") as f:
        json.dump(out, f, indent=2)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

_STACK_ORDER = [
    GTInheritanceStatus.consistent,
    GTInheritanceStatus.inconsistent,
    GTInheritanceStatus.incomplete,
    GTInheritanceStatus.missing,
    GTInheritanceStatus.uncertain,
    GTInheritanceStatus.non_informative,
    GTInheritanceStatus.no_call,
]

_STATUS_HATCH = {
    GTInheritanceStatus.consistent: "",
    GTInheritanceStatus.inconsistent: "///",
    GTInheritanceStatus.incomplete: "---",
    GTInheritanceStatus.missing: "|||",
    GTInheritanceStatus.uncertain: "xxx",
    GTInheritanceStatus.non_informative: "...",
    GTInheritanceStatus.no_call: "\\\\\\",
}

_STATUS_LABEL = {
    GTInheritanceStatus.consistent: "consistent",
    GTInheritanceStatus.inconsistent: "inconsistent",
    GTInheritanceStatus.incomplete: "incomplete",
    GTInheritanceStatus.missing: "missing",
    GTInheritanceStatus.uncertain: "uncertain",
    GTInheritanceStatus.non_informative: "non-informative",
    GTInheritanceStatus.no_call: "no-call",
}

_SVTYPE_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
_SVTYPE_LINESTYLES = ["-", "--", "-.", ":"]


def _sort_svtype_keys(keys: Iterable[str], *, include_all: bool = True) -> list[str]:
    """Return SV-type keys in canonical order: [all,] DEL, INS, <others alphabetically>."""
    key_set = set(keys)
    others = sorted(k for k in key_set if k not in ("all", "DEL", "INS"))
    ordered: list[str] = []
    if include_all and "all" in key_set:
        ordered.append("all")
    for k in ["DEL", "INS"]:
        if k in key_set:
            ordered.append(k)
    ordered.extend(others)
    return ordered


def _lighten_color(rgba, amount=0.45):
    return (
        rgba[0] + (1 - rgba[0]) * amount,
        rgba[1] + (1 - rgba[1]) * amount,
        rgba[2] + (1 - rgba[2]) * amount,
        rgba[3],
    )


def _darken_color(rgba, amount=0.3):
    return (
        rgba[0] * (1 - amount),
        rgba[1] * (1 - amount),
        rgba[2] * (1 - amount),
        rgba[3],
    )


def _consistency_rate(status_dict: dict[GTInheritanceStatus, int]) -> float:
    """Consistency over fully observed, informative trios, as a percentage.

    No-calls are outside the denominator; :func:`no_call_rate` reports how large
    that exclusion is.
    """
    n_cons = status_dict.get(GTInheritanceStatus.consistent, 0)
    n_inf = n_informative_trios(status_dict)
    return n_cons / n_inf * 100 if n_inf > 0 else float("nan")


def _save_fig(fig, prefix: str, suffix: str) -> None:
    fig.savefig(f"{prefix}_{suffix}.svg", bbox_inches="tight")
    fig.savefig(f"{prefix}_{suffix}.png", dpi=150, bbox_inches="tight")
    log.info(f"Saved {prefix}_{suffix}.svg and .png")


def plot_mendelian_results(
    all_dict_stats: dict[str, dict[str, dict[GTInheritanceStatus, int]]],
    all_size_dict_stats: dict[str, dict[str, dict[str, dict[GTInheritanceStatus, int]]]]
    | None,
    fig_prefix: str,
    plot_color: str = "#1f77b4",
) -> None:
    """Generate stacked bar charts and size-stratified scatter plots.

    Produces SVG and PNG for each figure:
      {fig_prefix}_stacked_combined  – ALL / DEL / INS only, with consistency % labels
      {fig_prefix}_stacked_per_svtype – one subplot per sample, bars = SV types
      {fig_prefix}_size_scatter – rows=samples, cols=SV types, consistency rate vs. size bin
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    import numpy as np

    base_rgba = mcolors.to_rgba(plot_color)
    light_rgba = _lighten_color(base_rgba)
    edge_rgba = _darken_color(base_rgba)

    samples = sorted(all_dict_stats.keys())
    all_svt: set[str] = set()
    for stats in all_dict_stats.values():
        all_svt.update(stats.keys())
    n_samples = len(samples)

    # ---- Plot 1: combined stacked bar (ALL, DEL, INS only) ----
    combined_keys = [k for k in ["all", "DEL", "INS"] if k in all_svt or k == "all"]
    n_comb = len(combined_keys)
    bar_width = 0.7 / max(n_samples, 1)
    x = np.arange(n_comb, dtype=float)
    fig1, ax1 = plt.subplots(figsize=(max(6, n_comb * 2.5 + 2), 6))

    legend_handles: list[tuple] = []
    for si, sample in enumerate(samples):
        offset = (si - (n_samples - 1) / 2) * bar_width
        bottoms = np.zeros(n_comb)
        for status in _STACK_ORDER:
            values = np.array(
                [
                    all_dict_stats[sample].get(svt, _empty_status_dict()).get(status, 0)
                    for svt in combined_keys
                ],
                dtype=float,
            )
            fc = base_rgba if status == GTInheritanceStatus.consistent else light_rgba
            bars = ax1.bar(
                x + offset,
                values,
                bar_width,
                bottom=bottoms,
                color=fc,
                edgecolor=edge_rgba,
                hatch=_STATUS_HATCH[status],
                linewidth=0.5,
            )
            # Label consistency bars with percentage
            if status == GTInheritanceStatus.consistent:
                for bi, (bar, val) in enumerate(zip(bars, values, strict=False)):
                    if val > 0:
                        rate = _consistency_rate(
                            all_dict_stats[sample].get(
                                combined_keys[bi], _empty_status_dict()
                            )
                        )
                        if not (rate != rate):  # not NaN
                            ax1.text(
                                bar.get_x() + bar.get_width() / 2,
                                bottoms[bi] + val / 2,
                                f"{rate:.1f}%",
                                ha="center",
                                va="center",
                                fontsize="x-small",
                                fontweight="bold",
                                color="white",
                            )
            if si == 0:
                legend_handles.append((bars[0], _STATUS_LABEL[status]))
            bottoms += values

    ax1.set_xticks(x)
    ax1.set_xticklabels(combined_keys)
    ax1.set_ylabel("Count")
    ax1.set_xlabel("SV Type")
    title1 = "Mendelian Consistency (ALL / DEL / INS)"
    if n_samples > 1:
        title1 += f"  (samples L\u2192R: {', '.join(samples)})"
    ax1.set_title(title1)
    ax1.legend(
        [h for h, _ in legend_handles],
        [label for _, label in legend_handles],
        loc="upper right",
        fontsize="small",
    )
    plt.tight_layout()
    _save_fig(fig1, fig_prefix, "stacked_combined")
    plt.close(fig1)

    # ---- Plot 2: per-svtype stacked bars (one subplot per sample) ----
    svtype_noall = _sort_svtype_keys(all_svt, include_all=False)
    n_svt_noall = len(svtype_noall)
    n_rows2 = max(n_samples, 1)
    fig2, axes2 = plt.subplots(
        n_rows2,
        1,
        figsize=(max(6, n_svt_noall * 1.2 + 2), 4.5 * n_rows2),
        squeeze=False,
    )

    for si, sample in enumerate(samples):
        ax = axes2[si][0]
        x_s = np.arange(n_svt_noall, dtype=float)
        bottoms = np.zeros(n_svt_noall)
        for status in _STACK_ORDER:
            values = np.array(
                [
                    all_dict_stats[sample].get(svt, _empty_status_dict()).get(status, 0)
                    for svt in svtype_noall
                ],
                dtype=float,
            )
            fc = base_rgba if status == GTInheritanceStatus.consistent else light_rgba
            bars = ax.bar(
                x_s,
                values,
                0.6,
                bottom=bottoms,
                color=fc,
                edgecolor=edge_rgba,
                hatch=_STATUS_HATCH[status],
                linewidth=0.5,
                label=_STATUS_LABEL[status] if si == 0 else None,
            )
            # Label consistency bars with percentage
            if status == GTInheritanceStatus.consistent:
                for bi, (bar, val) in enumerate(zip(bars, values, strict=False)):
                    if val > 0:
                        rate = _consistency_rate(
                            all_dict_stats[sample].get(
                                svtype_noall[bi], _empty_status_dict()
                            )
                        )
                        if not (rate != rate):  # not NaN
                            ax.text(
                                bar.get_x() + bar.get_width() / 2,
                                bottoms[bi] + val / 2,
                                f"{rate:.1f}%",
                                ha="center",
                                va="center",
                                fontsize="x-small",
                                fontweight="bold",
                                color="white",
                            )
            bottoms += values
        ax.set_xticks(x_s)
        ax.set_xticklabels(svtype_noall)
        ax.set_ylabel("Count")
        ax.set_title(sample)

    for idx in range(n_samples, n_rows2):
        axes2[idx][0].set_visible(False)

    fig2.legend(
        *axes2[0][0].get_legend_handles_labels(),
        loc="lower right",
        fontsize="small",
    )
    plt.tight_layout()
    _save_fig(fig2, fig_prefix, "stacked_per_svtype")
    plt.close(fig2)

    # ---- Plot 3: size-stratified scatter matrix (rows=samples, cols=SV types) ----
    if all_size_dict_stats:
        size_svt: set[str] = set()
        for ss in all_size_dict_stats.values():
            size_svt.update(ss.keys())
        size_svtype_keys = _sort_svtype_keys(size_svt)
        n_sv_cols = len(size_svtype_keys)
        n_samp_rows = max(n_samples, 1)

        x_bins = np.arange(len(SIZE_BIN_LABELS), dtype=float)

        fig3, axes3 = plt.subplots(
            n_samp_rows,
            n_sv_cols,
            figsize=(4.5 * n_sv_cols, 3.5 * n_samp_rows),
            squeeze=False,
        )

        for si, sample in enumerate(samples):
            size_data = all_size_dict_stats.get(sample, {})
            for ti, svtype in enumerate(size_svtype_keys):
                ax = axes3[si][ti]
                if svtype not in size_data:
                    ax.set_visible(False)
                    continue
                rates = np.array([
                    _consistency_rate(size_data[svtype][bl]) for bl in SIZE_BIN_LABELS
                ])
                valid = ~np.isnan(rates)
                ax.plot(
                    x_bins[valid],
                    rates[valid],
                    marker="o",
                    linestyle="-",
                    color=plot_color,
                    linewidth=1.5,
                    markersize=7,
                )
                # Add percentage labels to each point
                for xi, yi in zip(x_bins[valid], rates[valid], strict=False):
                    ax.text(
                        xi,
                        yi + 2,
                        f"{yi:.1f}%",
                        ha="center",
                        va="bottom",
                        fontsize="x-small",
                        fontweight="bold",
                    )
                ax.set_xticks(x_bins)
                ax.set_xticklabels(
                    SIZE_BIN_LABELS, fontsize="x-small", rotation=30, ha="right"
                )
                ax.set_ylim(-5, 105)
                ax.grid(axis="y", alpha=0.3)
                if si == 0:
                    ax.set_title(svtype)
                if ti == 0:
                    ax.set_ylabel(f"{sample}\nConsistency (%)", fontsize="small")
                else:
                    ax.set_ylabel("")
                if si == n_samp_rows - 1:
                    ax.set_xlabel("SV Size", fontsize="small")

        # Hide any leftover axes
        for si in range(n_samples, n_samp_rows):
            for ti in range(n_sv_cols):
                axes3[si][ti].set_visible(False)

        plt.tight_layout()
        _save_fig(fig3, fig_prefix, "size_scatter")
        plt.close(fig3)


def add_mndl_format_header(header_lines: list[str]) -> list[str]:
    """Add MNDL FORMAT header if not present.

    Args:
        header_lines: List of VCF header lines

    Returns:
        Updated list of header lines with MNDL FORMAT added
    """
    mndl_header = (
        '##FORMAT=<ID=MNDL,Number=1,Type=String,Description="Mendelian consistency">'
    )

    # Check if MNDL header already exists
    if any("ID=MNDL" in line for line in header_lines):
        return header_lines

    # Find the position to insert (after other FORMAT lines, before #CHROM line)
    insert_pos = len(header_lines) - 1  # Default: before last line (#CHROM)

    for i, line in enumerate(header_lines):
        if line.startswith("#CHROM"):
            insert_pos = i
            break

    # Insert the MNDL header
    result = header_lines[:insert_pos] + [mndl_header] + header_lines[insert_pos:]
    return result


def annotate_vcf_with_mendelian_consistency(
    input_vcf: PosixPath,
    output_vcf: PosixPath,
    ped: PosixPath,
    passonly: bool = False,
    min_size: int = 0,
    regions: dict[str, IntervalTree] | None = None,
    svtypes: set[str] | None = None,
    missing_as_reference: bool = False,
) -> None:
    """Annotate VCF file with Mendelian consistency information.

    Args:
        input_vcf: Path to input VCF file
        output_vcf: Path to output VCF file (can be .vcf or .vcf.gz)
        ped: Path to pedigree file
        passonly: If True, only annotate PASS variants
        missing_as_reference: reproduce the pre-v0.3 reading of ``./.``; see
            ``--legacy-missing-as-reference``.
    """
    log.info(f"Annotating VCF file {input_vcf} with Mendelian consistency")

    # Parse pedigree to get trios
    dict_ped = parse_ped(ped)

    # Read VCF file
    reader = vcfpy.Reader.from_path(str(input_vcf))

    # Get all sample names from VCF
    all_samples = reader.header.samples.names

    # Find trios that exist in the VCF
    trios = []
    child_samples = set()  # Track which samples are children
    for sample in all_samples:
        if sample in dict_ped:
            father_name, mother_name, _ = dict_ped[sample]
            if father_name in all_samples and mother_name in all_samples:
                trios.append([sample, father_name, mother_name])
                child_samples.add(sample)
            else:
                log.warning(
                    f"Trio for {sample} incomplete in VCF (missing father or mother)"
                )

    if len(trios) == 0:
        log.warning("No complete trios found in VCF file")
        return

    log.info(f"Found {len(trios)} complete trios in VCF")

    # Create a temporary file for writing

    # Update header with MNDL FORMAT field
    header_line = vcfpy.FormatHeaderLine(
        "FORMAT",
        '<ID=MNDL,Number=1,Type=String,Description="Mendelian consistency"',
        {
            "ID": "MNDL",
            "Number": 1,
            "Type": "String",
            "Description": "Mendelian consistency",
        },
    )
    insertion_idx = []
    for i, line in enumerate(reader.header.lines):
        if type(line) == vcfpy.FormatHeaderLine:
            insertion_idx.append(i)
    # insert after last FORMAT line
    insert_pos = insertion_idx[-1]
    new_lines = (
        reader.header.lines[: insert_pos + 1]
        + [header_line]
        + reader.header.lines[insert_pos + 1 :]
    )
    new_header = vcfpy.Header(lines=new_lines, samples=reader.header.samples)

    try:
        tmp_output = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".vcf")
        tmp_output_writer = vcfpy.Writer.from_path(tmp_output.name, new_header)
        # Process variants
        for variant in reader:
            # Skip if passonly and variant is not PASS
            if passonly and variant.FILTER != ["PASS"]:
                continue

            # Skip SVTYPEs not in the requested set
            if (
                svtypes is not None
                and "SVTYPE" in variant.INFO
                and variant.INFO["SVTYPE"] not in svtypes
            ):
                continue

            # Skip variants outside requested regions
            if regions is not None and not variant_overlaps_regions(variant, regions):
                continue

            # Skip variants smaller than min_size
            if min_size > 0:
                svlen = variant.INFO.get("SVLEN", [0])
                if isinstance(svlen, list):
                    svlen = svlen[0]
                if abs(int(svlen)) < min_size:
                    continue

            # Check Mendelian consistency for each trio
            trio_statuses = {}
            for child, father, mother in trios:
                status = is_variant_inconsistent(
                    variant=variant,
                    names_trio=[child, father, mother],
                    missing_as_reference=missing_as_reference,
                )
                trio_statuses[child] = status.name

            # Update FORMAT field to include MNDL
            if "MNDL" not in variant.FORMAT:
                variant.FORMAT.append("MNDL")

            # Update each sample's call data
            for call in variant.calls:
                if call.sample in trio_statuses:
                    # This sample is a child in a trio
                    call.data["MNDL"] = trio_statuses[call.sample]
                else:
                    # This sample is not a child (parent or unrelated)
                    call.data["MNDL"] = "."

            # Write the variant
            tmp_output_writer.write_record(variant)
        tmp_output_writer.close()

        # Sort and compress if needed
        if str(output_vcf).endswith(".vcf.gz"):
            # Sort to temporary file
            tmp_sorted = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
            cmd_sort = f"bcftools sort {tmp_output.name} -o {tmp_sorted.name}"
            subprocess.check_call(cmd_sort.split())

            # Compress
            with open(output_vcf, "wb") as f:
                cmd_zip = f"bgzip -c {tmp_sorted.name}"
                subprocess.check_call(cmd_zip.split(), stdout=f)

            # Index
            cmd_index = f"tabix -f -p vcf {str(output_vcf)}"
            subprocess.check_call(cmd_index.split())

            # Clean up temporary sorted file
            Path(tmp_sorted.name).unlink()
        else:
            # Just sort
            cmd_sort = f"bcftools sort {tmp_output.name} -o {str(output_vcf)}"
            subprocess.check_call(cmd_sort.split())

        log.info(f"Successfully annotated VCF written to {output_vcf}")

    finally:
        # Clean up temporary file
        Path(tmp_output.name).unlink(missing_ok=True)


def mendelian_consistency(
    input: PosixPath,
    output: PosixPath,
    ped: PosixPath | None = None,
    passonly: bool = False,
    tsv: PosixPath | None = None,
    latex: PosixPath | None = None,
    json_out: PosixPath | None = None,
    bed: PosixPath | None = None,
    vcf: PosixPath | None = None,
    min_size: int = 0,
    regions: PosixPath | None = None,
    svtypes: set[str] | None = None,
    size_stratification: bool = True,
    fig: PosixPath | None = None,
    plot_color: str = "#1f77b4",
    missing_as_reference: bool = False,
) -> None:
    """Run the Mendelian consistency check over *input* and write the reports.

    With *missing_as_reference* the pre-v0.3 reading of a no-call is restored so
    that previously published numbers reproduce exactly; see
    ``--legacy-missing-as-reference``.
    """
    if missing_as_reference:
        log.warning(
            "--legacy-missing-as-reference: missing genotypes (./.) are counted "
            "as observed homozygous reference calls.  The consistency figure "
            "produced under this flag is the pre-v0.3 one and does not "
            "distinguish a technical dropout from an observed reference allele."
        )
    regions_tree: dict[str, IntervalTree] | None = (
        load_regions_from_bed(regions) if regions is not None else None
    )

    # If VCF annotation is requested, do that and return
    if vcf is not None:
        if ped is None:
            raise ValueError("--ped is required when --vcf is specified")
        annotate_vcf_with_mendelian_consistency(
            input_vcf=input,
            output_vcf=vcf,
            ped=ped,
            passonly=passonly,
            min_size=min_size,
            regions=regions_tree,
            svtypes=svtypes,
            missing_as_reference=missing_as_reference,
        )

    variants: list[vcfpy.Record] = parse_variants(path_vcf=input, passonly=passonly)
    if ped:
        trios = get_trios_from_variants(variants, ped)
    else:
        sample = variants[0].calls[0].sample
        father = variants[0].calls[1].sample
        mother = variants[0].calls[2].sample
        trios = [[sample, father, mother]]
    all_dict_stats: dict[str, dict[str, dict[GTInheritanceStatus, int]]] = {}
    all_size_dict_stats: (
        dict[str, dict[str, dict[str, dict[GTInheritanceStatus, int]]]] | None
    ) = {} if size_stratification else None
    with open(output, "w") as f:
        for sample, father, mother in sorted(trios, key=lambda x: x[0]):
            dict_stats, size_stats = variants_consistency_stats(
                names_trio=[sample, father, mother],
                variants=variants,
                output=output,
                bed=bed,
                min_size=min_size,
                regions=regions_tree,
                svtypes=svtypes,
                size_stratification=size_stratification,
                missing_as_reference=missing_as_reference,
            )
            all_dict_stats[sample] = dict_stats
            if all_size_dict_stats is not None and size_stats is not None:
                all_size_dict_stats[sample] = size_stats
            line_consistency = dict_stats_to_line(
                sample, dict_stats, size_stats if size_stratification else None
            )
            print(line_consistency, file=f)
            log.info(line_consistency)
    if tsv:
        all_dict_stats_to_tsv(all_dict_stats, tsv, all_size_dict_stats)
    if latex:
        all_dict_stats_to_latex(all_dict_stats, latex, all_size_dict_stats)
    if json_out:
        all_dict_stats_to_json(all_dict_stats, json_out, all_size_dict_stats)
    if fig:
        plot_mendelian_results(
            all_dict_stats,
            all_size_dict_stats,
            str(fig),
            plot_color=plot_color,
        )


def run(args, **kwargs):
    svtypes = set(args.svtypes.split(",")) if args.svtypes else None
    mendelian_consistency(
        input=args.input,
        output=args.output,
        ped=args.ped,
        passonly=args.passonly,
        tsv=args.tsv,
        latex=args.latex,
        json_out=args.json,
        bed=args.bed,
        vcf=args.vcf,
        min_size=args.min_size,
        regions=args.regions,
        svtypes=svtypes,
        size_stratification=args.size_stratification,
        fig=args.fig,
        plot_color=args.plot_color,
        missing_as_reference=args.missing_as_reference,
    )
    return


def get_parser():
    parser = argparse.ArgumentParser(
        description="Check Mendelian consistency in a vcf file. The first column needs to be the child, the second the father and the third the mother"
    )
    parser.add_argument(
        "-i", "--input", type=PosixPath, required=True, help="Path to vcf(.gz) file"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=PosixPath,
        required=True,
        help="Path to the output summary file",
    )
    parser.add_argument(
        "-p",
        "--ped",
        type=PosixPath,
        required=False,
        help="Path to .ped file. Required when --vcf is used.",
    )
    parser.add_argument(
        "--passonly",
        action="store_true",
        help="If set, only variants that have the FILTER field set to 'PASS' will be considered.",
    )
    parser.add_argument(
        "--tsv",
        required=False,
        default=None,
        help="If a path is provided, a tsv of the results will be saved to this path",
    )
    parser.add_argument(
        "--latex",
        required=False,
        default=None,
        help="If a path is provided, a latex table of the results will be saved to this path",
    )
    parser.add_argument(
        "--json",
        type=PosixPath,
        required=False,
        default=None,
        help="If a path is provided, a JSON file of the results will be saved to this path",
    )
    parser.add_argument(
        "--bed",
        type=PosixPath,
        required=False,
        default=None,
        help="If a path is provided, a BED file with variant features (SVLEN,SVTYPE,PRECISION) and consistency status will be created for each sample",
    )
    parser.add_argument(
        "--vcf",
        type=PosixPath,
        required=False,
        default=None,
        help="Path to an output vcf file that will contain annotated variants with Mendelian consistency status in the MNDL FORMAT field. When specified, the script only performs VCF annotation (no summary statistics).",
    )
    parser.add_argument(
        "--min-size",
        type=int,
        required=False,
        default=0,
        help="Minimum variant size (absolute SVLEN). Variants smaller than this are skipped. Default: 0 (no filtering).",
    )
    parser.add_argument(
        "--regions",
        type=PosixPath,
        required=False,
        default=None,
        help="Path to a BED file. If provided, only variants overlapping at least one region are included in the analysis.",
    )
    parser.add_argument(
        "--svtypes",
        type=str,
        required=False,
        default=None,
        help="Comma-separated list of SVTYPEs to include (e.g. 'DEL,INS'). If not specified, all SVTYPEs present in the VCF are analysed.",
    )
    parser.add_argument(
        "--no-size-stratification",
        action="store_false",
        dest="size_stratification",
        default=True,
        help=(
            "If set, disable size stratification (on by default). "
            "By default, additional sub-tables with counts stratified by variant size are produced. "
            "Size is taken from SVLEN when available, otherwise from REF length (DEL), "
            "ALT length (INS), or max(REF, ALT) length (INV). "
            f"Bins: {', '.join(SIZE_BIN_LABELS)}. BNDs are excluded from size stratification."
        ),
    )
    parser.add_argument(
        "--fig",
        type=PosixPath,
        required=False,
        default=None,
        help=(
            "Path prefix for figure output. Produces SVG and PNG files: "
            "<prefix>_stacked_combined.{svg,png}, "
            "<prefix>_stacked_per_svtype.{svg,png}, and "
            "<prefix>_size_scatter.{svg,png} (when size stratification is on)."
        ),
    )
    parser.add_argument(
        "--legacy-missing-as-reference",
        action="store_true",
        dest="missing_as_reference",
        default=False,
        help=(
            "Restore the pre-v0.3 reading of a missing genotype: count './.' as "
            "an observed homozygous reference call instead of reporting it as a "
            "'no_call'.  Only for reproducing numbers published before the "
            "no-call outcome existed -- under this flag a technical dropout in "
            "the child scores as a correct inheritance and a dropout in a parent "
            "of a carrier child scores as a Mendelian violation."
        ),
    )
    parser.add_argument(
        "--plot-color",
        type=str,
        required=False,
        default="#1f77b4",
        help=(
            "Base color for bars and markers in plots (any matplotlib color string). "
            "Stacked bars use the solid color for 'consistent' and hatched lighter shades "
            "for other statuses. Scatter markers use this color directly. "
            "Default: '#1f77b4' (matplotlib tab:blue)."
        ),
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
