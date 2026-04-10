"""Tests for pairwise distance/similarity, BND adjustment, and outlier detection."""

import numpy as np

from svirlpool.localassembly.consensus_lib import (
    DirectedSignals,
    _directed_distance,
    _has_interior_bnd,
    detect_outlier_reads,
    importance_densities_from_ava_signals,
    pairwise_similarity_matrix,
    size_damping,
)
from svirlpool.util.datatypes import SVsignal

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _ins(
    ref_start: int,
    ref_end: int,
    size: int,
    read_start: int = 0,
    read_end: int | None = None,
) -> SVsignal:
    return SVsignal(
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=read_start,
        read_end=read_end if read_end is not None else read_start + size,
        size=size,
        sv_type=0,
    )


def _del(
    ref_start: int,
    ref_end: int,
    size: int,
    read_start: int = 0,
    read_end: int | None = None,
) -> SVsignal:
    return SVsignal(
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=read_start,
        read_end=read_end if read_end is not None else read_start,
        size=size,
        sv_type=1,
    )


def _bnd(ref_start: int, size: int, sv_type: int = 3) -> SVsignal:
    """Create a BNDL (type=3) or BNDR (type=4) signal."""
    return SVsignal(
        ref_start=ref_start,
        ref_end=ref_start,
        read_start=0,
        read_end=0,
        size=size,
        sv_type=sv_type,
    )


# ===========================================================================
# size_damping
# ===========================================================================


class TestSizeDamping:
    def test_small_signal_is_strongly_damped(self) -> None:
        assert size_damping(10.0, s0=30.0, alpha=1.5) < 0.20

    def test_large_signal_is_nearly_undamped(self) -> None:
        assert size_damping(200.0, s0=30.0, alpha=1.5) > 0.99

    def test_transition_is_monotonic(self) -> None:
        values = [size_damping(float(s), s0=30.0, alpha=1.5) for s in range(1, 500)]
        assert all(a <= b for a, b in zip(values, values[1:], strict=True))

    def test_zero_size_returns_zero(self) -> None:
        assert size_damping(0.0) == 0.0


# ===========================================================================
# importance_densities_from_ava_signals (updated: BNDs excluded, damping)
# ===========================================================================


class TestImportanceDensitiesUpdated:
    def test_bnd_signals_are_excluded_from_density(self) -> None:
        """A track built only from BND signals should be all zeros."""
        ava_signals = {
            "readA": [
                _bnd(ref_start=500, size=300, sv_type=3),
                _bnd(ref_start=2000, size=500, sv_type=4),
            ]
        }
        densities = importance_densities_from_ava_signals(ava_signals)
        # ref_length = max(ref_end) = 2000, but all BND => density should be all zeros
        # Actually BND has ref_start==ref_end so max ref_end = 2000
        track = densities.get("readA")
        if track is not None:
            assert track.max() == 0.0

    def test_small_signal_damped_versus_large(self) -> None:
        """A 15bp INS should produce a smaller peak than a 200bp INS."""
        ava_signals = {
            "readA": [
                _ins(ref_start=490, ref_end=510, size=15),
                _ins(ref_start=2490, ref_end=2510, size=200),
            ]
        }
        densities = importance_densities_from_ava_signals(
            ava_signals, size_densities=None
        )
        track = densities["readA"]
        peak_small = track[490:520].max()
        peak_large = track[2490:2520].max()
        assert peak_large > peak_small

    def test_backward_compat_no_size_densities(self) -> None:
        """Without size_densities, large signals still get scale=1.0."""
        ava_signals = {"readA": [_ins(ref_start=450, ref_end=550, size=200)]}
        densities = importance_densities_from_ava_signals(
            ava_signals, size_densities=None
        )
        assert densities["readA"].max() > 0


# ===========================================================================
# _has_interior_bnd
# ===========================================================================


class TestHasInteriorBnd:
    def test_bnd_at_start_is_not_interior(self) -> None:
        signals = [_bnd(ref_start=50, size=300, sv_type=3)]
        assert not _has_interior_bnd(signals, ref_length=5000, bnd_margin=100)

    def test_bnd_at_end_is_not_interior(self) -> None:
        signals = [_bnd(ref_start=4950, size=300, sv_type=4)]
        assert not _has_interior_bnd(signals, ref_length=5000, bnd_margin=100)

    def test_bnd_in_middle_is_interior(self) -> None:
        signals = [_bnd(ref_start=2500, size=300, sv_type=3)]
        assert _has_interior_bnd(signals, ref_length=5000, bnd_margin=100)

    def test_bnd_exactly_at_margin_is_not_interior(self) -> None:
        signals = [_bnd(ref_start=100, size=300, sv_type=3)]
        assert not _has_interior_bnd(signals, ref_length=5000, bnd_margin=100)

    def test_empty_signals(self) -> None:
        assert not _has_interior_bnd([], ref_length=5000)

    def test_non_bnd_signals_are_ignored(self) -> None:
        signals = [_ins(ref_start=2500, ref_end=2600, size=100)]
        assert not _has_interior_bnd(signals, ref_length=5000)


# ===========================================================================
# _directed_distance
# ===========================================================================


class TestDirectedDistance:
    def test_bnd_signals_are_excluded(self) -> None:
        signals = [_bnd(ref_start=500, size=1000, sv_type=3)]
        assert _directed_distance(signals, None, None, 30.0, 1.5) == 0.0

    def test_basic_ins_contributes(self) -> None:
        signals = [
            _ins(ref_start=500, ref_end=600, size=100, read_start=50, read_end=150)
        ]
        d = _directed_distance(signals, None, None, 30.0, 1.5)
        assert d > 0

    def test_importance_weighting_amplifies(self) -> None:
        signals = [
            _ins(ref_start=500, ref_end=600, size=100, read_start=50, read_end=150)
        ]
        # With uniform density of 2.0 at the signal positions
        ref_density = np.full(1000, 2.0)
        query_density = np.full(200, 2.0)
        d_weighted = _directed_distance(signals, ref_density, query_density, 30.0, 1.5)
        d_unweighted = _directed_distance(signals, None, None, 30.0, 1.5)
        # 2.0 * 2.0 = 4x amplification
        np.testing.assert_allclose(d_weighted, d_unweighted * 4.0, rtol=1e-6)

    def test_small_signal_contributes_less(self) -> None:
        small = [_ins(ref_start=500, ref_end=510, size=12)]
        large = [_ins(ref_start=500, ref_end=600, size=200)]
        d_small = _directed_distance(small, None, None, 30.0, 1.5)
        d_large = _directed_distance(large, None, None, 30.0, 1.5)
        assert d_small < d_large


# ===========================================================================
# pairwise_similarity_matrix
# ===========================================================================


class TestPairwiseSimilarityMatrix:
    def _make_three_read_scenario(self):
        """Three reads: A and B are similar, C has an interior BND to A."""
        read_lengths = {"A": 5000, "B": 5000, "C": 5000}
        # A->B and B->A: small INS signals (similar reads)
        ds_ab = DirectedSignals(
            query_name="A",
            ref_name="B",
            signals=[
                _ins(
                    ref_start=1000,
                    ref_end=1100,
                    size=50,
                    read_start=1000,
                    read_end=1050,
                )
            ],
            ref_length=5000,
        )
        ds_ba = DirectedSignals(
            query_name="B",
            ref_name="A",
            signals=[
                _ins(
                    ref_start=1000,
                    ref_end=1100,
                    size=50,
                    read_start=1000,
                    read_end=1050,
                )
            ],
            ref_length=5000,
        )
        # A->C: interior BND (reads are very different)
        ds_ac = DirectedSignals(
            query_name="A",
            ref_name="C",
            signals=[_bnd(ref_start=2500, size=1000, sv_type=3)],
            ref_length=5000,
        )
        ds_ca = DirectedSignals(
            query_name="C",
            ref_name="A",
            signals=[
                _ins(
                    ref_start=1000,
                    ref_end=1100,
                    size=50,
                    read_start=1000,
                    read_end=1050,
                )
            ],
            ref_length=5000,
        )
        return [ds_ab, ds_ba, ds_ac, ds_ca], read_lengths

    def test_diagonal_is_one(self) -> None:
        ds, rl = self._make_three_read_scenario()
        sim, names = pairwise_similarity_matrix(
            ds, {}, rl, all_read_names=["A", "B", "C"]
        )
        for i in range(len(names)):
            assert sim[i, i] == 1.0

    def test_symmetry(self) -> None:
        ds, rl = self._make_three_read_scenario()
        sim, names = pairwise_similarity_matrix(
            ds, {}, rl, all_read_names=["A", "B", "C"]
        )
        np.testing.assert_array_equal(sim, sim.T)

    def test_interior_bnd_zeroes_similarity(self) -> None:
        """A<->C should be 0 because A->C has an interior BND."""
        ds, rl = self._make_three_read_scenario()
        sim, names = pairwise_similarity_matrix(
            ds, {}, rl, all_read_names=["A", "B", "C"]
        )
        ai, ci = names.index("A"), names.index("C")
        assert sim[ai, ci] == 0.0

    def test_no_alignment_gives_zero_similarity(self) -> None:
        """B and C have no alignments between them."""
        ds, rl = self._make_three_read_scenario()
        sim, names = pairwise_similarity_matrix(
            ds, {}, rl, all_read_names=["A", "B", "C"]
        )
        bi, ci = names.index("B"), names.index("C")
        assert sim[bi, ci] == 0.0

    def test_similar_reads_have_high_similarity(self) -> None:
        """A and B have small symmetric signals -> positive similarity."""
        ds, rl = self._make_three_read_scenario()
        sim, names = pairwise_similarity_matrix(
            ds, {}, rl, all_read_names=["A", "B", "C"]
        )
        ai, bi = names.index("A"), names.index("B")
        assert sim[ai, bi] > 0.3

    def test_length_difference_reduces_similarity(self) -> None:
        """Two reads with identical signals but different lengths should have
        lower similarity than two with the same length."""
        sig = _ins(
            ref_start=1000, ref_end=1100, size=100, read_start=1000, read_end=1100
        )
        ds_equal = [
            DirectedSignals("X", "Y", [sig], ref_length=5000),
            DirectedSignals("Y", "X", [sig], ref_length=5000),
        ]
        ds_unequal = [
            DirectedSignals("X", "Y", [sig], ref_length=5000),
            DirectedSignals("Y", "X", [sig], ref_length=10000),
        ]
        rl_equal = {"X": 5000, "Y": 5000}
        rl_unequal = {"X": 10000, "Y": 5000}

        sim_eq, _ = pairwise_similarity_matrix(
            ds_equal, {}, rl_equal, all_read_names=["X", "Y"]
        )
        sim_uneq, _ = pairwise_similarity_matrix(
            ds_unequal, {}, rl_unequal, all_read_names=["X", "Y"]
        )
        assert sim_eq[0, 1] > sim_uneq[0, 1]

    def test_bnd_at_read_edge_does_not_zero_similarity(self) -> None:
        """A BND at position 50 on a 5000bp ref (within margin) should NOT
        trigger the zeroing."""
        ds = [
            DirectedSignals(
                "A",
                "B",
                [
                    _bnd(ref_start=50, size=300, sv_type=3),
                    _ins(
                        ref_start=1000,
                        ref_end=1100,
                        size=100,
                        read_start=1000,
                        read_end=1100,
                    ),
                ],
                ref_length=5000,
            ),
            DirectedSignals(
                "B",
                "A",
                [
                    _ins(
                        ref_start=1000,
                        ref_end=1100,
                        size=100,
                        read_start=1000,
                        read_end=1100,
                    )
                ],
                ref_length=5000,
            ),
        ]
        rl = {"A": 5000, "B": 5000}
        sim, names = pairwise_similarity_matrix(ds, {}, rl, all_read_names=["A", "B"])
        assert sim[0, 1] > 0.0

    def test_empty_directed_signals(self) -> None:
        rl = {"A": 5000, "B": 5000}
        sim, names = pairwise_similarity_matrix([], {}, rl, all_read_names=["A", "B"])
        assert sim[0, 1] == 0.0
        assert sim[0, 0] == 1.0


# ===========================================================================
# detect_outlier_reads
# ===========================================================================


class TestDetectOutlierReads:
    def test_uniform_reads_no_outliers(self) -> None:
        """All reads have similar lengths and are well-connected."""
        names = [f"r{i}" for i in range(10)]
        rl = dict.fromkeys(names, 5000)
        sim = np.ones((10, 10), dtype=np.float64) * 0.5
        np.fill_diagonal(sim, 1.0)
        outliers = detect_outlier_reads(sim, names, rl, n_clusters=2)
        assert outliers == set()

    def test_length_only_outlier_is_flagged(self) -> None:
        """A read with strongly outlier length is flagged even when well-connected,
        because detect_outlier_reads uses OR: failing either criterion is sufficient."""
        names = ["r0", "r1", "r2", "r3", "r4"]
        rl = {"r0": 5000, "r1": 5000, "r2": 5000, "r3": 5000, "r4": 50000}
        sim = np.ones((5, 5), dtype=np.float64) * 0.5
        np.fill_diagonal(sim, 1.0)
        outliers = detect_outlier_reads(sim, names, rl, n_clusters=2)
        assert "r4" in outliers

    def test_both_criteria_flags_outlier(self) -> None:
        """A read with outlier length OR low connectivity IS flagged."""
        names = ["r0", "r1", "r2", "r3", "r4"]
        rl = {"r0": 5000, "r1": 5000, "r2": 5000, "r3": 5000, "r4": 50000}
        sim = np.ones((5, 5), dtype=np.float64) * 0.5
        np.fill_diagonal(sim, 1.0)
        # Disconnect r4
        sim[4, :] = 0.0
        sim[:, 4] = 0.0
        sim[4, 4] = 1.0
        outliers = detect_outlier_reads(sim, names, rl, n_clusters=2)
        assert "r4" in outliers

    def test_empty_input(self) -> None:
        assert detect_outlier_reads(np.zeros((0, 0)), [], {}, n_clusters=1) == set()

    def test_coherent_length_groups_not_flagged(self) -> None:
        """Two groups of reads with distinct but internally consistent lengths
        should NOT be flagged as outliers — only truly isolated reads should."""
        # Group A: lengths ~1200 (18 reads), Group B: lengths ~1660 (10 reads),
        # plus one singleton outlier at 162.
        names_a = [f"a{i}" for i in range(18)]
        names_b = [f"b{i}" for i in range(10)]
        singleton = "outlier"
        names = [singleton] + names_a + names_b
        rl = {singleton: 162}
        rl.update({r: 1200 + i * 2 for i, r in enumerate(names_a)})
        rl.update({r: 1660 + i * 2 for i, r in enumerate(names_b)})

        n = len(names)
        sim = np.ones((n, n), dtype=np.float64) * 0.5
        np.fill_diagonal(sim, 1.0)

        outliers = detect_outlier_reads(sim, names, rl, n_clusters=1)
        # The singleton should be flagged
        assert singleton in outliers
        # Neither coherent group should be flagged
        for r in names_a + names_b:
            assert r not in outliers, f"{r} should not be an outlier"

    def test_low_connectivity_outlier_median_based(self) -> None:
        """A read with very few connections relative to the median is flagged."""
        names = [f"r{i}" for i in range(10)]
        rl = dict.fromkeys(names, 5000)
        sim = np.ones((10, 10), dtype=np.float64) * 0.5
        np.fill_diagonal(sim, 1.0)
        # Disconnect r9 from almost everyone
        sim[9, :] = 0.0
        sim[:, 9] = 0.0
        sim[9, 9] = 1.0
        outliers = detect_outlier_reads(sim, names, rl, n_clusters=2)
        assert "r9" in outliers
