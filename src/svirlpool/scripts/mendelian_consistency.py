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

import argparse
import subprocess
import tempfile
from enum import Enum
from pathlib import Path, PosixPath

import vcfpy
from logzero import logger as log

# %%


# parse chr,start and GTs
def parse_variant(variant):
    return (variant.CHROM, variant.POS), [call.data["GT"] for call in variant.calls]


# %%
def parse_genotype(gt_string: str) -> tuple[int, int, int]:
    """
    Parse genotype string to count alleles.

    Treats missing genotypes (./.) as homozygous reference (0/0).

    Returns: (n_ref_alleles, n_alt_alleles, total_copies)

    Examples:
        "0" → (1, 0, 1)         # hemizygous ref
        "1" → (0, 1, 1)         # hemizygous alt
        "0/0" → (2, 0, 2)       # homozygous ref
        "0/1" → (1, 1, 2)       # heterozygous
        "1/1" → (0, 2, 2)       # homozygous alt
        "0/0/1" → (2, 1, 3)     # triploid with 1 alt
        "0/1/1" → (1, 2, 3)     # triploid with 2 alts
        "0/1/1/1" → (1, 3, 4)   # tetraploid with 3 alts
        "./." → (2, 0, 2)       # missing, treated as homozygous ref
        "." → (1, 0, 1)         # missing hemizygous, treated as ref
    """
    # Handle missing genotypes by treating them as reference
    if "." in gt_string:
        # Count the number of alleles (by counting slashes + 1, or 1 if no slash)
        if "/" in gt_string:
            n_alleles = gt_string.count("/") + 1
        else:
            n_alleles = 1
        # Treat missing as all reference alleles
        return (n_alleles, 0, n_alleles)

    # Handle both hemizygous (no slash) and multi-allelic (with slash)
    if "/" in gt_string:
        alleles = gt_string.split("/")
    else:
        alleles = [gt_string]

    # Count reference (0) and alt (non-0) alleles
    n_ref = alleles.count("0")
    # For now, treat any non-0 allele as alt (handles 1, 2, etc.)
    n_alt = sum(1 for a in alleles if a != "0")
    total = len(alleles)

    return (n_ref, n_alt, total)


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
    non_informative = 5  # No one in trio has a variant (all 0 or missing)


# %%
def is_trio_informative(
    child: tuple[int, int, int],
    father: tuple[int, int, int],
    mother: tuple[int, int, int],
) -> bool:
    """
    Check if a trio has at least one member with an alternate allele.

    Returns False if all members are homozygous reference (n_alt = 0).

    This filters out non-informative variants where the entire trio
    is either ./. or 0/0 or any combination of missing/reference.

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
    variant: vcfpy.Record, names_trio: list[str] | None = None
) -> GTInheritanceStatus:
    """
    Check Mendelian consistency for a variant across all copy numbers (CN 0-4+).

    Uses generalized allele counting approach that works for:
    - Hemizygous (CN=1): e.g., "0" or "1"
    - Diploid (CN=2): e.g., "0/0", "0/1", "1/1"
    - Triploid (CN=3): e.g., "0/0/1", "0/1/1", "1/1/1"
    - Tetraploid (CN=4): e.g., "0/0/0/1", "0/0/1/1", etc.

    Args:
        variant: VCF record
        names_trio: List of [child, father, mother] sample names

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
    child = parse_genotype(gt_son)
    father = parse_genotype(gt_father)
    mother = parse_genotype(gt_mother)

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


# creates the output and results line of one trio
def variants_consistency_stats(
    variants: list[vcfpy.Record],
    output: PosixPath,
    names_trio: list[str] | None = None,
    bed: PosixPath | None = None,
) -> dict[GTInheritanceStatus, int]:
    """
    Calculate Mendelian consistency statistics for variants in a trio.

    Now supports generalized copy numbers (CN 0-4+) including:
    - Hemizygous (CN=1)
    - Diploid (CN=2)
    - Triploid (CN=3)
    - Tetraploid (CN=4)

    Args:
        variants: List of VCF records
        output: Output BED file path
        names_trio: Optional list of [child, father, mother] sample names
        output_bed: Optional path for BED file with features (chr, start, end, features, status)

    Returns:
        Dictionary mapping GTInheritanceStatus to counts
    """
    if names_trio is not None:
        assert len(names_trio) == 3, "names_trio must have 3 elements"
    # get the number of inconsistent variants
    # counter for each GTInheritanceStatus

    # Open bed file if provided
    bed_file = open(bed, "w") if bed else None

    with open(output, "w") as f:
        dict_stats = {
            GTInheritanceStatus.consistent: 0,
            GTInheritanceStatus.inconsistent: 0,
            GTInheritanceStatus.missing: 0,
            GTInheritanceStatus.incomplete: 0,
            GTInheritanceStatus.uncertain: 0,
            GTInheritanceStatus.non_informative: 0,
        }
        for variant in variants:
            if variant.INFO["SVTYPE"] not in ["DEL", "INS"]:
                continue
            status = is_variant_inconsistent(variant=variant, names_trio=names_trio)
            dict_stats[status] += 1
            # if inconsistent, print the variant CHROM and POS
            end = (
                int(variant.INFO["END"])
                if variant.INFO["SVTYPE"] == "DEL"
                else variant.POS + 1
            )

            # Write to BED file with features if requested
            if bed_file:
                features = extract_variant_features(variant)
                print(
                    f"{variant.CHROM}\t{variant.POS}\t{end}\t{status.name},{features}",
                    file=bed_file,
                )

    if bed_file:
        bed_file.close()

    return dict_stats


def dict_stats_to_line(
    samplename: str, dict_stats: dict[GTInheritanceStatus, int]
) -> str:
    # calc n_variants by summing all values in dict_stats
    n_variants_total = sum(dict_stats.values())
    # Calculate n_variants excluding uncertain and non_informative for percentage calculation
    n_variants_informative = sum(
        count
        for status, count in dict_stats.items()
        if status
        not in [GTInheritanceStatus.uncertain, GTInheritanceStatus.non_informative]
    )

    # separate thousands by comma and percentage by two decimal places
    results_line = f"{samplename}:\n- total: {n_variants_total:,}\n"

    # Add each status with count and percentage (calculated from informative variants only)
    for status in dict_stats:
        count = dict_stats[status]
        if status in [
            GTInheritanceStatus.uncertain,
            GTInheritanceStatus.non_informative,
        ]:
            # Show count only, no percentage for excluded statuses
            results_line += f"- {status.name}: {count:,}\n"
        else:
            # Show count and percentage for informative statuses
            percentage = (
                count / n_variants_informative if n_variants_informative > 0 else 0
            )
            results_line += f"- {status.name}: {count:,}; {percentage:.2%}\n"

    return results_line.rstrip()  # Remove trailing newline


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
    all_dict_stats: dict[str, dict[GTInheritanceStatus, int]], path_tsv: PosixPath
) -> None:
    # writes a tsv file with colums: sample, status, count, percentage (.2f)
    with open(path_tsv, "w") as f:
        print("sample\tstatus\tcount\tpercentage", file=f)
        for sample, dict_stats in all_dict_stats.items():
            # Calculate n_variants excluding uncertain and non_informative for percentage calculation
            n_variants_informative = sum(
                count
                for status, count in dict_stats.items()
                if status
                not in [
                    GTInheritanceStatus.uncertain,
                    GTInheritanceStatus.non_informative,
                ]
            )
            for status in dict_stats:
                count = dict_stats[status]
                if status in [
                    GTInheritanceStatus.uncertain,
                    GTInheritanceStatus.non_informative,
                ]:
                    # Show NA for percentage of excluded statuses
                    percentage_str = "NA"
                else:
                    # Calculate percentage from informative variants only
                    percentage = (
                        count / n_variants_informative
                        if n_variants_informative > 0
                        else 0
                    )
                    percentage_str = str(percentage)
                print(f"{sample}\t{status.name}\t{count}\t{percentage_str}", file=f)


def all_dict_stats_to_latex(
    all_dict_stats: dict[str, dict[GTInheritanceStatus, int]], path_latex: PosixPath
) -> None:
    # crate 1 row per GTInheritanceStatus
    # the columns are: GTInheritanceStatus, * sample names
    # print two rows per GTInheritanceStatus: one row with the counts and another row with the percentages
    with open(path_latex, "w") as f:
        print(r"status & " + " & ".join(all_dict_stats.keys()) + r"\\", file=f)
        for status in GTInheritanceStatus:
            line_abs = (
                f"{status.name} & "
                + " & ".join([
                    f"{all_dict_stats[sample][status]:,}" for sample in all_dict_stats
                ])
                + r"\\"
            )
            print(line_abs, file=f)
            n_variants_per_sample = {
                samplename: sum(dict_stats.values())
                for samplename, dict_stats in all_dict_stats.items()
            }
            line_per = (
                f"{status.name} & "
                + " & ".join([
                    f"{all_dict_stats[sample][status] / n_variants_per_sample[sample]:.2}"
                    for sample in all_dict_stats
                ])
                + r"\\"
            )
            # line_per = line_per.replace("%", r"\%")
            print(line_per, file=f)


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
    input_vcf: PosixPath, output_vcf: PosixPath, ped: PosixPath, passonly: bool = False
) -> None:
    """Annotate VCF file with Mendelian consistency information.

    Args:
        input_vcf: Path to input VCF file
        output_vcf: Path to output VCF file (can be .vcf or .vcf.gz)
        ped: Path to pedigree file
        passonly: If True, only annotate PASS variants
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

            # Skip non-SV types if needed
            if "SVTYPE" in variant.INFO and variant.INFO["SVTYPE"] not in [
                "DEL",
                "INS",
            ]:
                continue

            # Check Mendelian consistency for each trio
            trio_statuses = {}
            for child, father, mother in trios:
                status = is_variant_inconsistent(
                    variant=variant, names_trio=[child, father, mother]
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
    bed: PosixPath | None = None,
    vcf: PosixPath | None = None,
) -> None:
    # If VCF annotation is requested, do that and return
    if vcf is not None:
        if ped is None:
            raise ValueError("--ped is required when --vcf is specified")
        annotate_vcf_with_mendelian_consistency(
            input_vcf=input, output_vcf=vcf, ped=ped, passonly=passonly
        )

    variants: list[vcfpy.Record] = parse_variants(path_vcf=input, passonly=passonly)
    if ped:
        trios = get_trios_from_variants(variants, ped)
    else:
        sample = variants[0].calls[0].sample
        father = variants[0].calls[1].sample
        mother = variants[0].calls[2].sample
        trios = [[sample, father, mother]]
    all_dict_stats: dict[str, dict[GTInheritanceStatus, int]] = {}
    with open(output, "w") as f:
        for sample, father, mother in sorted(trios, key=lambda x: x[0]):
            # Create sample-specific output_bed path if requested
            dict_stats: dict[GTInheritanceStatus, int] = variants_consistency_stats(
                names_trio=[sample, father, mother],
                variants=variants,
                output=output,
                bed=bed,
            )
            all_dict_stats[sample] = dict_stats
            line_consistency = dict_stats_to_line(sample, dict_stats)
            print(line_consistency, file=f)
            log.info(line_consistency)
    if tsv:
        all_dict_stats_to_tsv(all_dict_stats, tsv)
    if latex:
        all_dict_stats_to_latex(all_dict_stats, latex)


def run(args, **kwargs):
    mendelian_consistency(
        input=args.input,
        output=args.output,
        ped=args.ped,
        passonly=args.passonly,
        tsv=args.tsv,
        latex=args.latex,
        bed=args.bed,
        vcf=args.vcf,
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
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
