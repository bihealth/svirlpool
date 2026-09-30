"""The vertical merge's size test: one floored relative bound,
|a - b| <= max(floor, fraction * ref(a, b)), with ref = min, max or hmean."""

import math
import random

import pytest
from test_can_merge_svComposites import _make_insertion_composite

from svirlpool.svcalling import svcomposite_merging as m
from svirlpool.svcalling.multisample_sv_calling import get_parser


@pytest.fixture
def rule():
    """Set the module-level reference and floor, and restore the defaults."""
    saved = (m.SIZE_TOLERANCE_REFERENCE, m.SIZE_TOLERANCE_FLOOR)

    def set_rule(reference="max", floor=0.0):
        m.SIZE_TOLERANCE_REFERENCE = reference
        m.SIZE_TOLERANCE_FLOOR = floor

    yield set_rule
    m.SIZE_TOLERANCE_REFERENCE, m.SIZE_TOLERANCE_FLOOR = saved


def test_defaults_are_the_previous_fraction_gate():
    assert m.SIZE_TOLERANCE_REFERENCE == "max"
    assert m.SIZE_TOLERANCE_FLOOR == 0.0


@pytest.mark.parametrize(
    "reference, expected",
    [("min", 10.0), ("max", 11.0), ("hmean", 0.1 * 2 * 100 * 110 / 210)],
)
def test_the_reference_size(rule, reference, expected):
    rule(reference)
    assert m.size_tolerance(100, 110, 0.1) == pytest.approx(expected)
    assert m.size_tolerance(110, 100, 0.1) == pytest.approx(expected)
    assert m.size_tolerance(-100, -110, 0.1) == pytest.approx(expected)


@pytest.mark.parametrize("reference", m.SIZE_TOLERANCE_REFERENCES)
def test_the_floor_is_absolute(rule, reference):
    rule(reference, floor=5)
    # 6% of ~50 bp is 3 bp: the floor of 5 bp decides
    assert m.sizes_similar(50, 55, 0.06)
    assert not m.sizes_similar(50, 56, 0.06)
    # above the floor the relative bound decides
    assert m.sizes_similar(1000, 1050, 0.06)
    assert not m.sizes_similar(1000, 1100, 0.06)


@pytest.mark.parametrize("reference", m.SIZE_TOLERANCE_REFERENCES)
def test_identical_sizes_are_similar_at_zero_tolerance(rule, reference):
    rule(reference)
    for size in (0, 1, 50, 5000):
        assert m.sizes_similar(size, size, 0.0)


def test_the_three_references_are_one_ratio_test():
    """Without a floor each reference accepts exactly the pairs with max/min <= R.

    R = 1 + f (min), 1 / (1 - f) (max), f + sqrt(f^2 + 1) (hmean), so at
    matched R the three verdicts coincide on every pair.
    """
    ratio = 1.1234567  # no ratio of two integers below 10^4 lands on it
    fractions = {
        "min": ratio - 1,
        "max": 1 - 1 / ratio,
        "hmean": (ratio**2 - 1) / (2 * ratio),
    }
    for reference, f in fractions.items():
        m.SIZE_TOLERANCE_REFERENCE = reference
        rng = random.Random(7)
        for _ in range(20000):
            a, b = rng.randint(1, 10000), rng.randint(1, 10000)
            want = max(a, b) / min(a, b) <= ratio
            assert m.sizes_similar(a, b, f) is want, (reference, a, b)
    m.SIZE_TOLERANCE_REFERENCE = "max"


def test_hmean_fraction_maps_to_its_ratio():
    f = 0.1
    ratio = f + math.sqrt(f * f + 1)
    assert (ratio - 1) == pytest.approx(f * 2 * ratio / (1 + ratio))


def test_an_unknown_reference_is_refused(rule):
    rule("mean")
    with pytest.raises(ValueError, match="SIZE_TOLERANCE_REFERENCE"):
        m.size_tolerance(100, 110, 0.1)


def _pair(size_a, size_b, *, sibling=False, dist_a=None, dist_b=None):
    a = _make_insertion_composite(
        size=size_a, size_distortions=dist_a, samplename="sample1", consensusID="1.0"
    )
    b = _make_insertion_composite(
        size=size_b,
        size_distortions=dist_b,
        samplename="sample1" if sibling else "sample2",
        consensusID="1.1" if sibling else "2.0",
    )
    return a, b


def test_noise_populations_do_not_enter_the_size_test(rule):
    """Wide, overlapping distortion populations no longer rescue a size gap."""
    rule()
    wide = {f"r{i}": v for i, v in enumerate([-80, -50, -20, 0, 20, 50, 80])}
    for dist in (None, wide):
        a, b = _pair(500, 600, dist_a=dist, dist_b=dist)
        assert not m._similar_size(a, b, 0.1)
        a, b = _pair(500, 540, dist_a=dist, dist_b=dist)
        assert m._similar_size(a, b, 0.1)


def test_detail_reports_the_tolerance(rule):
    rule("min", floor=3)
    a, b = _pair(100, 108)
    r = m._similar_size_detail(a, b, 0.05)
    assert r == {
        "similar": False,
        "size_tolerance": 5.0,
        "size_reference": "min",
        "size_floor": 3,
    }


def test_the_sibling_test_uses_the_same_rule(rule, monkeypatch):
    """Two haplotypes of one sample: SIBLING_SIZE_TOLERANCE, same reference and floor."""
    monkeypatch.setattr(m, "HAPLOTYPE_AWARE_MERGE", True)
    monkeypatch.setattr(m, "SIBLING_SIZE_TOLERANCE", 0.1)
    a, b = _pair(50, 57, sibling=True)
    assert m._assembly_relation(a, b) == "sibling_assembly"
    rule("max", floor=0)
    assert m._haplotype_gate(a, b, "INS") is False  # 7 > 0.1 * 57
    rule("max", floor=8)
    assert m._haplotype_gate(a, b, "INS") is None  # left to proximity and k-mers


def test_cli_options():
    args = get_parser().parse_args(
        [
            "--input", "x.db",
            "--output", "o.vcf",
            "--reference", "r.fa",
            "--size-tolerance-reference", "hmean",
            "--size-tolerance-floor", "5",
        ]
    )  # fmt: skip
    assert args.size_tolerance_reference == "hmean"
    assert args.size_tolerance_floor == 5.0
    assert args.max_cohens_d is None
    with pytest.raises(SystemExit):
        get_parser().parse_args(
            ["--input", "x.db", "--output", "o.vcf", "--reference", "r.fa",
             "--size-gates", "fraction"]
        )  # fmt: skip
