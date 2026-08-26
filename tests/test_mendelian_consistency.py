"""Tests for :mod:`svirlpool.analysis.mendelian_consistency`.

The module had no tests before N9.  The suite below pins three things:

1. every genotype string form the caller can put into a VCF (haploid,
   diploid, triploid, tetraploid, phased, multi-allelic, missing, partially
   missing and absent);
2. the Mendelian rule itself over the classic trio cases, built as real
   ``vcfpy`` records read from a real VCF text rather than as hand-made
   objects;
3. the N9 policy: a trio that contains a no-call is reported as its own
   outcome (``no_call``) instead of being silently read as an observed
   homozygous reference, and the pre-N9 reading stays reachable behind
   ``missing_as_reference=True`` / ``--legacy-missing-as-reference``.
"""

from pathlib import Path

import pytest
import vcfpy

from svirlpool.analysis.barplot_mendelian_consistency import (
    create_plot,
    create_table,
)
from svirlpool.analysis.mendelian_consistency import (
    _STACK_ORDER,
    _STATUS_HATCH,
    _STATUS_LABEL,
    GTInheritanceStatus,
    _consistency_rate,
    _empty_status_dict,
    all_dict_stats_to_json,
    all_dict_stats_to_tsv,
    dict_stats_to_line,
    is_no_call,
    is_trio_informative,
    is_variant_inconsistent,
    parse_genotype,
    parse_variants,
    variants_consistency_stats,
)

TRIO = ["child", "father", "mother"]

_VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DEL,Description="Deletion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother
"""


def write_trio_vcf(
    tmp_path: Path,
    trios: list[tuple[str, str, str]],
    *,
    svtype: str = "INS",
    svlen: int = 500,
    name: str = "trio.vcf",
) -> Path:
    """Write a VCF with one INS record per ``(child, father, mother)`` triple."""
    path = tmp_path / name
    lines = [_VCF_HEADER]
    for i, (child, father, mother) in enumerate(trios):
        pos = 1000 * (i + 1)
        info = f"SVTYPE={svtype};SVLEN={svlen}"
        if svtype == "DEL":
            info += f";END={pos + svlen}"
        lines.append(
            f"chr1\t{pos}\t.\tN\t<{svtype}>\t.\tPASS\t{info}\tGT"
            f"\t{child}\t{father}\t{mother}\n"
        )
    path.write_text("".join(lines))
    return path


def read_records(path: Path) -> list[vcfpy.Record]:
    return parse_variants(path)


def status_of(
    tmp_path: Path, child: str, father: str, mother: str, **kwargs
) -> GTInheritanceStatus:
    """Build a one-record VCF for the trio and return its consistency status."""
    path = write_trio_vcf(
        tmp_path, [(child, father, mother)], name=f"one_{abs(hash((child, father, mother)))}.vcf"
    )
    record = read_records(path)[0]
    return is_variant_inconsistent(variant=record, names_trio=TRIO, **kwargs)


# ---------------------------------------------------------------------------
# parse_genotype: every form the VCF can actually contain
# ---------------------------------------------------------------------------


class TestParseGenotypeObserved:
    """Observed genotypes: (n_ref, n_alt, ploidy), and n_ref + n_alt == ploidy."""

    @pytest.mark.parametrize(
        ("gt", "expected"),
        [
            # haploid -- emitted by the caller at copy number 0 and 1
            ("0", (1, 0, 1)),
            ("1", (0, 1, 1)),
            # diploid -- the common case
            ("0/0", (2, 0, 2)),
            ("0/1", (1, 1, 2)),
            ("1/1", (0, 2, 2)),
            # triploid / tetraploid -- emitted at higher copy number
            ("0/0/0", (3, 0, 3)),
            ("0/0/1", (2, 1, 3)),
            ("0/1/1", (1, 2, 3)),
            ("1/1/1", (0, 3, 3)),
            ("0/0/0/0", (4, 0, 4)),
            ("0/0/0/1", (3, 1, 4)),
            ("0/1/1/1", (1, 3, 4)),
            ("1/1/1/1", (0, 4, 4)),
            # multi-allelic: any non-zero allele index counts as alt
            ("2", (0, 1, 1)),
            ("1/2", (0, 2, 2)),
            ("0/2", (1, 1, 2)),
        ],
    )
    def test_observed_forms(self, gt, expected):
        assert parse_genotype(gt) == expected

    @pytest.mark.parametrize(
        ("gt", "expected"),
        [("0|0", (2, 0, 2)), ("0|1", (1, 1, 2)), ("1|0", (1, 1, 2)), ("1|1", (0, 2, 2))],
    )
    def test_phased_separator_is_understood(self, gt, expected):
        """``|`` is a valid VCF allele separator and must not be read as one allele.

        Before N9 ``"0|1"`` parsed to ``(0, 1, 1)`` -- a hemizygous alt -- because
        only ``/`` was recognised as a separator.
        """
        assert parse_genotype(gt) == expected

    @pytest.mark.parametrize(
        "gt",
        ["0", "1", "0/0", "0/1", "1/1", "0/0/0", "0/1/1", "1/1/1/1", "0|1", "1/2"],
    )
    def test_observed_genotypes_are_not_no_calls(self, gt):
        assert is_no_call(parse_genotype(gt)) is False

    @pytest.mark.parametrize(
        "gt", ["0", "1", "0/0", "0/1", "1/1", "0/0/0", "0/1/1", "0|1", "1/2"]
    )
    def test_allele_counts_sum_to_ploidy(self, gt):
        n_ref, n_alt, ploidy = parse_genotype(gt)
        assert n_ref + n_alt == ploidy


class TestParseGenotypeMissing:
    """Missing genotypes are 'no information', not an observed reference."""

    @pytest.mark.parametrize(
        ("gt", "ploidy"),
        [("./.", 2), (".", 1), ("././.", 3), ("./././.", 4)],
    )
    def test_fully_missing_is_a_no_call_of_the_right_ploidy(self, gt, ploidy):
        parsed = parse_genotype(gt)
        assert is_no_call(parsed) is True
        assert parsed == (0, 0, ploidy)

    @pytest.mark.parametrize("gt", ["0/.", "./0", "1/.", "./1", "0|."])
    def test_partially_missing_is_a_no_call(self, gt):
        """A half-observed genotype cannot support a Mendelian determination."""
        assert is_no_call(parse_genotype(gt)) is True

    def test_absent_genotype_is_a_no_call(self):
        """vcfpy decodes a bare ``.`` GT field as ``None``; it must not raise."""
        assert is_no_call(parse_genotype(None)) is True

    def test_empty_genotype_is_a_no_call(self):
        assert is_no_call(parse_genotype("")) is True

    def test_no_call_carries_no_alt_alleles(self):
        _, n_alt, _ = parse_genotype("./.")
        assert n_alt == 0


class TestParseGenotypeLegacy:
    """The pre-N9 reading stays reachable so published numbers reproduce."""

    @pytest.mark.parametrize(
        ("gt", "expected"),
        [("./.", (2, 0, 2)), (".", (1, 0, 1)), ("././.", (3, 0, 3)), ("0/.", (2, 0, 2))],
    )
    def test_missing_as_reference_reproduces_the_old_mapping(self, gt, expected):
        assert parse_genotype(gt, missing_as_reference=True) == expected

    def test_legacy_missing_is_not_flagged_as_a_no_call(self):
        assert is_no_call(parse_genotype("./.", missing_as_reference=True)) is False

    @pytest.mark.parametrize("gt", ["0/0", "0/1", "1/1", "0", "1", "0/0/1"])
    def test_legacy_flag_does_not_touch_observed_genotypes(self, gt):
        assert parse_genotype(gt, missing_as_reference=True) == parse_genotype(gt)


# ---------------------------------------------------------------------------
# is_trio_informative
# ---------------------------------------------------------------------------


class TestIsTrioInformative:
    def test_all_reference_is_not_informative(self):
        gts = [parse_genotype("0/0")] * 3
        assert is_trio_informative(*gts) is False

    @pytest.mark.parametrize("carrier", [0, 1, 2])
    def test_a_single_carrier_makes_the_trio_informative(self, carrier):
        gts = [parse_genotype("0/0") for _ in range(3)]
        gts[carrier] = parse_genotype("0/1")
        assert is_trio_informative(*gts) is True

    def test_a_no_call_alone_does_not_make_a_trio_informative(self):
        """A no-call carries no observed alt allele, so it cannot be evidence."""
        gts = [parse_genotype("./."), parse_genotype("0/0"), parse_genotype("0/0")]
        assert is_trio_informative(*gts) is False


# ---------------------------------------------------------------------------
# The consistency rule over classic trio cases (real records)
# ---------------------------------------------------------------------------


class TestDiploidRule:
    @pytest.mark.parametrize(
        ("child", "father", "mother"),
        [
            ("0/1", "0/1", "0/0"),  # one transmitting parent
            ("0/1", "0/0", "0/1"),
            ("0/1", "1/1", "0/0"),
            ("1/1", "0/1", "0/1"),  # both parents transmit an alt
            ("1/1", "1/1", "1/1"),
            ("0/0", "1/1", "0/0"),  # child reference: always possible
            ("0/0", "0/1", "0/1"),
        ],
    )
    def test_mendelian_consistent_trios(self, tmp_path, child, father, mother):
        assert (
            status_of(tmp_path, child, father, mother)
            is GTInheritanceStatus.consistent
        )

    @pytest.mark.parametrize(
        ("child", "father", "mother"),
        [
            ("0/1", "0/0", "0/0"),  # de novo alt
            ("1/1", "0/0", "0/0"),
            ("1/1", "0/1", "0/0"),  # only one parent can transmit an alt
            ("1/1", "0/0", "1/1"),
        ],
    )
    def test_mendelian_violating_trios(self, tmp_path, child, father, mother):
        assert (
            status_of(tmp_path, child, father, mother)
            is GTInheritanceStatus.inconsistent
        )

    def test_all_reference_trio_is_non_informative(self, tmp_path):
        assert (
            status_of(tmp_path, "0/0", "0/0", "0/0")
            is GTInheritanceStatus.non_informative
        )


class TestHemizygousAndPolyploidRules:
    def test_hemizygous_alt_from_a_reference_mother_is_inconsistent(self, tmp_path):
        assert status_of(tmp_path, "1", "0/1", "0/0") is GTInheritanceStatus.inconsistent

    def test_hemizygous_alt_from_a_carrier_mother_is_consistent(self, tmp_path):
        assert status_of(tmp_path, "1", "0/0", "0/1") is GTInheritanceStatus.consistent

    def test_haploid_reference_child_is_non_informative_with_reference_parents(
        self, tmp_path
    ):
        assert (
            status_of(tmp_path, "0", "0/0", "0/0")
            is GTInheritanceStatus.non_informative
        )

    def test_triploid_child_needs_enough_parental_alts(self, tmp_path):
        assert (
            status_of(tmp_path, "0/1/1", "0/1", "0/0")
            is GTInheritanceStatus.inconsistent
        )
        assert (
            status_of(tmp_path, "0/1/1", "0/1", "0/1") is GTInheritanceStatus.consistent
        )

    def test_triploid_reference_child_is_non_informative(self, tmp_path):
        """F9 emits ``0/0/0`` at copy number 3; it must read as a reference call."""
        assert (
            status_of(tmp_path, "0/0/0", "0/0", "0/0")
            is GTInheritanceStatus.non_informative
        )


# ---------------------------------------------------------------------------
# N9: no-calls in every position
# ---------------------------------------------------------------------------


class TestNoCallPolicy:
    @pytest.mark.parametrize(
        ("child", "father", "mother"),
        [
            ("./.", "0/1", "0/0"),  # child missing
            ("0/1", "./.", "0/0"),  # father missing
            ("0/1", "0/0", "./."),  # mother missing
            ("./.", "./.", "0/0"),  # two missing
            ("0/1", "./.", "./."),
            ("./.", "./.", "./."),  # all missing
            ("0/0", "./.", "0/0"),  # missing parent, nobody observed carrying
            ("./.", "0/0", "0/0"),
            (".", "0/1", "0/0"),  # haploid missing
            ("0/.", "0/1", "0/0"),  # partially missing
        ],
    )
    def test_any_no_call_in_the_trio_yields_no_call(
        self, tmp_path, child, father, mother
    ):
        assert status_of(tmp_path, child, father, mother) is GTInheritanceStatus.no_call

    def test_no_call_takes_precedence_over_non_informative(self, tmp_path):
        """'nobody was observed to carry it' and 'somebody was not observed'
        are different claims and must not share a bucket."""
        observed = status_of(tmp_path, "0/0", "0/0", "0/0")
        dropped = status_of(tmp_path, "0/0", "./.", "0/0")
        assert observed is GTInheritanceStatus.non_informative
        assert dropped is GTInheritanceStatus.no_call

    # -- the characterisation test: fails before the fix, passes after --

    def test_missing_child_no_longer_manufactures_a_consistent_verdict(self, tmp_path):
        """A no-call child used to be read as ``0/0``, and a ``0/0`` child is
        Mendelian-consistent with *any* parents.  The metric therefore scored a
        technical dropout as a correctly inherited genotype."""
        assert (
            status_of(tmp_path, "./.", "0/1", "0/0", missing_as_reference=True)
            is GTInheritanceStatus.consistent
        )
        assert status_of(tmp_path, "./.", "0/1", "0/0") is GTInheritanceStatus.no_call

    def test_missing_parent_no_longer_manufactures_an_inconsistent_verdict(
        self, tmp_path
    ):
        """The mirror image: a no-call parent used to be read as a
        non-transmitting ``0/0``, turning a dropout into a Mendelian violation."""
        assert (
            status_of(tmp_path, "0/1", "./.", "0/0", missing_as_reference=True)
            is GTInheritanceStatus.inconsistent
        )
        assert status_of(tmp_path, "0/1", "./.", "0/0") is GTInheritanceStatus.no_call

    def test_legacy_mode_reproduces_every_pre_n9_verdict(self, tmp_path):
        cases = {
            ("./.", "0/1", "0/0"): GTInheritanceStatus.consistent,
            ("0/1", "./.", "0/0"): GTInheritanceStatus.inconsistent,
            ("0/1", "0/0", "./."): GTInheritanceStatus.inconsistent,
            ("./.", "./.", "./."): GTInheritanceStatus.non_informative,
            ("0/0", "./.", "0/0"): GTInheritanceStatus.non_informative,
            ("1/1", "./.", "0/1"): GTInheritanceStatus.inconsistent,
        }
        for (child, father, mother), expected in cases.items():
            assert (
                status_of(tmp_path, child, father, mother, missing_as_reference=True)
                is expected
            ), (child, father, mother)


# ---------------------------------------------------------------------------
# Aggregation and reporting
# ---------------------------------------------------------------------------

# 1 consistent, 1 inconsistent, 1 non-informative, 3 no-calls
_MIXED_TRIO = [
    ("0/1", "0/1", "0/0"),  # consistent
    ("0/1", "0/0", "0/0"),  # inconsistent
    ("0/0", "0/0", "0/0"),  # non-informative
    ("./.", "0/1", "0/0"),  # no-call (child)
    ("0/1", "./.", "0/0"),  # no-call (father)
    ("0/1", "0/0", "./."),  # no-call (mother)
]


@pytest.fixture
def mixed_stats(tmp_path):
    path = write_trio_vcf(tmp_path, _MIXED_TRIO, name="mixed.vcf")
    dict_stats, size_stats = variants_consistency_stats(
        variants=read_records(path),
        output=tmp_path / "summary.txt",
        names_trio=TRIO,
    )
    return dict_stats, size_stats


class TestAggregation:
    def test_no_calls_are_counted_as_their_own_outcome(self, mixed_stats):
        dict_stats, _ = mixed_stats
        counts = dict_stats["all"]
        assert counts[GTInheritanceStatus.no_call] == 3
        assert counts[GTInheritanceStatus.consistent] == 1
        assert counts[GTInheritanceStatus.inconsistent] == 1
        assert counts[GTInheritanceStatus.non_informative] == 1

    def test_no_calls_are_out_of_the_headline_denominator(self, mixed_stats):
        """The headline consistency figure is over fully observed, informative
        trios: 1 of 2 here, not 1 of 5."""
        dict_stats, _ = mixed_stats
        assert _consistency_rate(dict_stats["all"]) == pytest.approx(50.0)

    def test_every_record_is_still_accounted_for(self, mixed_stats):
        dict_stats, _ = mixed_stats
        assert sum(dict_stats["all"].values()) == len(_MIXED_TRIO)

    def test_legacy_mode_restores_the_pre_n9_counts(self, tmp_path):
        path = write_trio_vcf(tmp_path, _MIXED_TRIO, name="mixed_legacy.vcf")
        dict_stats, _ = variants_consistency_stats(
            variants=read_records(path),
            output=tmp_path / "summary_legacy.txt",
            names_trio=TRIO,
            missing_as_reference=True,
        )
        counts = dict_stats["all"]
        assert counts[GTInheritanceStatus.no_call] == 0
        assert counts[GTInheritanceStatus.consistent] == 2  # +1 from the ./. child
        assert counts[GTInheritanceStatus.inconsistent] == 3  # +2 from the ./. parents
        assert counts[GTInheritanceStatus.non_informative] == 1

    def test_size_stratified_counts_carry_no_calls_too(self, mixed_stats):
        _, size_stats = mixed_stats
        assert size_stats is not None
        assert size_stats["all"]["350-1000"][GTInheritanceStatus.no_call] == 3

    def test_empty_status_dict_covers_every_status(self):
        assert set(_empty_status_dict()) == set(GTInheritanceStatus)


class TestReportedOutput:
    def test_no_call_rate_appears_in_the_text_summary(self, mixed_stats):
        dict_stats, size_stats = mixed_stats
        text = dict_stats_to_line("child", dict_stats, size_stats)
        assert "no_call" in text
        assert "no_call_rate" in text
        # 3 no-calls out of 6 records
        assert "50.00%" in text

    def test_tsv_carries_the_no_call_count_and_rate(self, mixed_stats, tmp_path):
        dict_stats, size_stats = mixed_stats
        tsv = tmp_path / "stats.tsv"
        all_dict_stats_to_tsv({"child": dict_stats}, tsv, {"child": size_stats})
        lines = tsv.read_text().splitlines()
        header = lines[0].split("\t")
        assert "denominator" in header
        rows = [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]
        no_call_rows = [
            r
            for r in rows
            if r["status"] == "no_call"
            and r["svtype"] == "all"
            and r["size_bin"] == "all"
        ]
        assert len(no_call_rows) == 1
        assert no_call_rows[0]["count"] == "3"
        assert float(no_call_rows[0]["percentage"]) == pytest.approx(0.5)
        assert no_call_rows[0]["denominator"] == "total"

    def test_tsv_headline_rows_are_still_over_informative_trios(
        self, mixed_stats, tmp_path
    ):
        dict_stats, size_stats = mixed_stats
        tsv = tmp_path / "stats2.tsv"
        all_dict_stats_to_tsv({"child": dict_stats}, tsv, {"child": size_stats})
        lines = tsv.read_text().splitlines()
        header = lines[0].split("\t")
        rows = [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]
        consistent = next(
            r
            for r in rows
            if r["status"] == "consistent"
            and r["svtype"] == "all"
            and r["size_bin"] == "all"
        )
        assert float(consistent["percentage"]) == pytest.approx(0.5)
        assert consistent["denominator"] == "informative"

    def test_json_carries_the_no_call_rate(self, mixed_stats, tmp_path):
        import json

        dict_stats, size_stats = mixed_stats
        path = tmp_path / "stats.json"
        all_dict_stats_to_json({"child": dict_stats}, path, {"child": size_stats})
        data = json.loads(path.read_text())
        entry = data["child"]["all"]
        assert entry["no_call"]["count"] == 3
        assert entry["no_call"]["percentage"] == pytest.approx(0.5)
        assert entry["no_call_rate"] == pytest.approx(0.5)


class TestPlottingMetadata:
    def test_every_status_is_plottable(self):
        assert set(_STACK_ORDER) == set(GTInheritanceStatus)
        for status in GTInheritanceStatus:
            assert status in _STATUS_HATCH
            assert status in _STATUS_LABEL


class TestBarplotConsumer:
    def test_barplot_still_reads_the_emitted_tsv(self, mixed_stats, tmp_path):
        import pandas as pd

        dict_stats, size_stats = mixed_stats
        tsv = tmp_path / "stats.tsv"
        all_dict_stats_to_tsv({"child": dict_stats}, tsv, {"child": size_stats})

        df = pd.read_csv(tsv, sep="\t")
        df = df[df["status"].isin(["consistent", "inconsistent"])]
        df = df[(df["svtype"] == "all") & (df["size_bin"] == "all")]
        df = df[["sample", "status", "count", "percentage"]]
        df["file_group"] = "arm"
        df["file_group"] = pd.Categorical(df["file_group"], categories=["arm"])
        fig = create_plot(df.copy(), title="t")
        assert fig is not None
        table = create_table(df.copy())
        assert "consistent" in table.columns
