"""The vertical merge's size test: |a - b| <= fraction * max(|a|, |b|)."""

import pytest
from test_can_merge_svComposites import _make_insertion_composite

from svirlpool.svcalling import svcomposite_merging as m
from svirlpool.svcalling.multisample_sv_calling import get_parser

_REQUIRED = ["--input", "x.db", "--output", "o.vcf", "--reference", "r.fa"]


def test_the_tolerance_is_a_fraction_of_the_larger_size():
    assert m.size_tolerance(100, 110, 0.1) == pytest.approx(11.0)
    assert m.size_tolerance(110, 100, 0.1) == pytest.approx(11.0)
    assert m.size_tolerance(-100, -110, 0.1) == pytest.approx(11.0)


def test_there_is_no_absolute_floor():
    # 10% of ~50 bp is ~5 bp, and nothing lifts it
    assert m.sizes_similar(50, 55, 0.1)
    assert not m.sizes_similar(50, 56, 0.1)
    assert m.sizes_similar(1000, 1100, 0.1)
    assert not m.sizes_similar(1000, 1112, 0.1)


def test_identical_sizes_are_similar_at_zero_tolerance():
    for size in (0, 1, 50, 5000):
        assert m.sizes_similar(size, size, 0.0)


def _pair(size_a, size_b, *, sibling=False):
    a = _make_insertion_composite(size=size_a, samplename="sample1", consensusID="1.0")
    b = _make_insertion_composite(
        size=size_b,
        samplename="sample1" if sibling else "sample2",
        consensusID="1.1" if sibling else "2.0",
    )
    return a, b


def test_sizes_alone_decide():
    """At fraction 0.1, 500 vs 600 is rejected and 500 vs 540 accepted."""
    a, b = _pair(500, 600)
    assert not m._similar_size(a, b, 0.1)
    a, b = _pair(500, 540)
    assert m._similar_size(a, b, 0.1)


def test_detail_reports_the_tolerance():
    a, b = _pair(100, 108)
    assert m._similar_size_detail(a, b, 0.05) == {
        "similar": False,
        "size_tolerance": pytest.approx(5.4),
    }


def test_the_sibling_test_uses_the_same_rule(monkeypatch):
    """Two haplotypes of one sample are compared at SIBLING_SIZE_TOLERANCE."""
    monkeypatch.setattr(m, "HAPLOTYPE_AWARE_MERGE", True)
    monkeypatch.setattr(m, "SIBLING_SIZE_TOLERANCE", 0.1)
    a, b = _pair(50, 57, sibling=True)
    assert m._assembly_relation(a, b) == "sibling_assembly"
    assert m._haplotype_gate(a, b, "INS") is False  # 7 > 0.1 * 57
    a, b = _pair(55, 57, sibling=True)
    assert m._haplotype_gate(a, b, "INS") is None  # left to proximity and k-mers


def test_cli_default_fraction_is_one_tenth():
    args = get_parser().parse_args(_REQUIRED)
    assert args.apriori_size_difference_fraction_tolerance == 0.1
    assert args.max_cohens_d is None


@pytest.mark.parametrize(
    "removed",
    [
        ["--size-tolerance-reference", "hmean"],
        ["--size-tolerance-floor", "5"],
        ["--size-gates", "fraction"],
    ],
)
def test_removed_options_are_refused(removed):
    with pytest.raises(SystemExit):
        get_parser().parse_args(_REQUIRED + removed)
