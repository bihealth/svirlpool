"""sv-calling --size-gates: the three gates of the vertical merge's size test,
each on its own, none, and the default (all three, the previous behaviour)."""

import numpy as np
import pytest
from test_can_merge_svComposites import _make_insertion_composite, _random_dna

from svirlpool.svcalling import svcomposite_merging as m
from svirlpool.svcalling.multisample_sv_calling import parse_size_gates

ALL = frozenset(m.SIZE_GATE_NAMES)


def _pair(size_a, size_b, dist_a, dist_b, random_sequence):
    kw_a, kw_b = {}, {}
    if random_sequence:
        kw_a["sequence"] = _random_dna(size_a, seed=11)
        kw_b["sequence"] = _random_dna(size_b, seed=12)
    a = _make_insertion_composite(
        size=size_a,
        size_distortions=dist_a,
        samplename="sample1",
        consensusID="1.0",
        **kw_a,
    )
    b = _make_insertion_composite(
        size=size_b,
        size_distortions=dist_b,
        samplename="sample2",
        consensusID="2.0",
        **kw_b,
    )
    return a, b


SPREAD = {f"r{i}": float(v) for i, v in enumerate([0, 2, 4, 6, 8, 10])}


def detail(a, b, gates, tol=0.06, factor=1.0, d=2.0):
    return m._similar_size_detail(a, b, tol, factor, d, gates=frozenset(gates))


def test_no_gate_passes_every_pair():
    a, b = _pair(100, 500, SPREAD, SPREAD, random_sequence=True)
    r = detail(a, b, [])
    assert r["similar"]
    assert not (r["fraction_similar"] or r["complexity_similar"] or r["noise_similar"])


def test_fraction_gate_alone():
    a, b = _pair(100, 104, SPREAD, SPREAD, random_sequence=True)
    assert detail(a, b, ["fraction"])["similar"]
    a, b = _pair(100, 120, SPREAD, SPREAD, random_sequence=True)
    r = detail(a, b, ["fraction"])
    assert not r["similar"]
    assert r["cohens_d_status"] == "not_reached"


def test_complexity_gate_alone_accepts_a_homopolymer_far_beyond_the_fraction():
    """The unbounded allowance: 100 vs 180 bp in a homopolymer passes."""
    a, b = _pair(100, 180, {"r1": 0.0}, {"r1": 0.0}, random_sequence=False)
    r = detail(a, b, ["complexity"])
    assert r["similar"] and r["complexity_similar"]
    assert not detail(a, b, ["fraction"])["similar"]


def test_complexity_gate_alone_never_computes_cohens_d():
    a, b = _pair(100, 300, SPREAD, SPREAD, random_sequence=True)
    r = detail(a, b, ["complexity"])
    assert not r["similar"]
    assert r["cohens_d_status"] == "not_reached"


def test_noise_gate_alone_is_not_shifted_by_the_complexity_allowance():
    """In a homopolymer the allowance would close the gap; alone, the noise gate
    sees the unshifted populations and rejects a separation of 8 sd."""
    a, b = _pair(100, 130, SPREAD, SPREAD, random_sequence=False)
    alone = detail(a, b, ["noise"])
    assert alone["cohens_d_status"] == "computed"
    assert not alone["similar"]
    assert detail(a, b, ["complexity", "noise"])["similar"]


def test_noise_gate_rejects_constant_populations_that_differ():
    """No spread and different values: d is infinite in the limit."""
    a, b = _pair(
        100, 130, {"r1": 0.0, "r2": 0.0}, {"r1": 0.0, "r2": 0.0}, random_sequence=True
    )
    r = detail(a, b, ["noise"])
    assert r["cohens_d"] is None
    assert r["cohens_d_status"] == "constant_different"
    assert not r["similar"]


@pytest.mark.parametrize(
    "dist",
    [
        {"r1": 0.0, "r2": 0.0, "r3": 0.0},
        {"r1": 0.0},  # a single read on each side is constant too
        {f"r{i}": 0.1 for i in range(19)},  # constant, np.std(ddof=1) != 0
    ],
)
def test_noise_gate_merges_constant_populations_that_coincide(dist):
    """No spread and the same value is a perfect match, not a missing one.

    The gate used to abstain here, rejecting about half of all identical-size
    cross-sample pairs of a trio call."""
    a, b = _pair(250, 250, dist, dist, random_sequence=True)
    r = detail(a, b, ["noise"])
    assert r["cohens_d"] is None
    assert r["cohens_d_status"] == "constant_equal"
    assert r["noise_similar"] and r["similar"]


def test_noise_gate_constants_compare_after_the_distortion_shift():
    """Equal sizes but different constant distortions are different values."""
    a, b = _pair(250, 250, {"r1": 0.0, "r2": 0.0}, {"r1": 3.0, "r2": 3.0}, True)
    r = detail(a, b, ["noise"])
    assert r["cohens_d_status"] == "constant_different"
    assert not r["similar"]


@pytest.mark.parametrize("sizes", [(100, 104), (100, 130), (100, 180), (300, 310)])
@pytest.mark.parametrize("random_sequence", [True, False])
@pytest.mark.parametrize("dist", [SPREAD, {"r1": 0.0}])
def test_default_is_all_three_and_the_old_tuple(sizes, random_sequence, dist):
    a, b = _pair(*sizes, dist, dist, random_sequence)
    r = detail(a, b, ALL)
    similar, fraction, population, cohens = m._similar_size(a, b, 0.06, 1.0, 2.0)
    assert m.SIZE_GATES == ALL
    assert similar == r["similar"]
    assert fraction == r["fraction_similar"]
    assert population == r["population_similar"]
    assert (cohens is None and r["cohens_d"] is None) or np.isnan(cohens) == np.isnan(
        r["cohens_d"]
    )


def test_parse_size_gates():
    assert parse_size_gates("none") == frozenset()
    assert parse_size_gates("fraction, noise") == frozenset({"fraction", "noise"})
    assert parse_size_gates("fraction,complexity,noise") == ALL
    with pytest.raises(ValueError):
        parse_size_gates("fraction,cohen")
    with pytest.raises(ValueError):
        parse_size_gates("")
