"""Tests for importance_densities_from_ava_signals in ava_visualization."""

import importlib.util
import sys
from pathlib import Path

import numpy as np

from svirlpool.localassembly.consensus_lib import (
    importance_densities_from_ava_signals,
    sv_size_densities_from_ava_signals,
)
from svirlpool.util.datatypes import SVsignal

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _ins(ref_start: int, ref_end: int, size: int) -> SVsignal:
    return SVsignal(
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=0,
        read_end=size,
        size=size,
        sv_type=0,  # INS
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_importance_densities_without_size_scaling_uses_unit_scale() -> None:
    """Without size_densities every signal contributes with scale=1.0.

    Place two INS signals at widely separated positions with one at each location.
    Without size scaling both bumps should have equal peak contributions (scale=1.0).
    """
    ava_signals = {
        "readA": [
            _ins(ref_start=450, ref_end=550, size=500),  # position ~500
            _ins(ref_start=2450, ref_end=2550, size=500),  # position ~2500
        ]
    }

    densities = importance_densities_from_ava_signals(ava_signals, size_densities=None)
    track = densities["readA"]

    peak_at_500 = track[480:520].max()
    peak_at_2500 = track[2480:2520].max()

    # Both bumps are identical shape (scale=1.0 for both) so their peaks must be equal
    np.testing.assert_allclose(peak_at_500, peak_at_2500, rtol=1e-6)


def test_importance_densities_size_scaling_reduces_rare_size_contribution() -> None:
    """With size_densities, a rare SV size contributes less than a common size.

    Three signals cluster near size=100 (common) and one signal at size=500 (rare),
    all at different positions so positional bumps do not overlap.

    After size-density scaling the density peak at the rare-size signal's position
    must be strictly less than the same position's peak computed without scaling.
    """
    ava_signals = {
        "readA": [
            # three common-size signals far from the rare one
            _ins(ref_start=450, ref_end=550, size=100),
            _ins(ref_start=650, ref_end=750, size=100),
            _ins(ref_start=850, ref_end=950, size=100),
            # one rare-size signal at a distinct position
            _ins(ref_start=2450, ref_end=2550, size=500),
        ]
    }

    size_densities = sv_size_densities_from_ava_signals(ava_signals)

    unscaled = importance_densities_from_ava_signals(ava_signals, size_densities=None)
    scaled = importance_densities_from_ava_signals(
        ava_signals, size_densities=size_densities
    )

    # The rare-size position (2500) should be reduced after scaling
    rare_peak_unscaled = unscaled["readA"][2480:2520].max()
    rare_peak_scaled = scaled["readA"][2480:2520].max()

    assert rare_peak_scaled < rare_peak_unscaled, (
        f"Expected scaled peak ({rare_peak_scaled:.4f}) < unscaled peak "
        f"({rare_peak_unscaled:.4f}) for the rare-size signal"
    )


def test_importance_densities_size_scaling_preserves_common_size_contribution() -> None:
    """Common-size signals (size=100) keep scale≈1.0 and are not reduced.

    The size density at size=100 is the maximum (three signals), so the
    normalised scale there is 1.0 — scaled and unscaled peaks must be equal.
    """
    ava_signals = {
        "readA": [
            _ins(ref_start=450, ref_end=550, size=100),
            _ins(ref_start=650, ref_end=750, size=100),
            _ins(ref_start=850, ref_end=950, size=100),
            _ins(ref_start=2450, ref_end=2550, size=500),
        ]
    }

    size_densities = sv_size_densities_from_ava_signals(ava_signals)

    unscaled = importance_densities_from_ava_signals(ava_signals, size_densities=None)
    scaled = importance_densities_from_ava_signals(
        ava_signals, size_densities=size_densities
    )

    # At the common-size position the scale should be 1.0, so peaks are equal
    common_peak_unscaled = unscaled["readA"][480:520].max()
    common_peak_scaled = scaled["readA"][480:520].max()

    np.testing.assert_allclose(common_peak_unscaled, common_peak_scaled, rtol=1e-6)


def test_importance_densities_returns_empty_for_empty_input() -> None:
    result = importance_densities_from_ava_signals({}, size_densities=None)
    assert result == {}


def test_importance_densities_bnd_signals_always_use_unit_scale() -> None:
    """BND signals (sv_type=3/4) are not in size_densities and must keep scale=1.0.

    The peak of a BND signal's positional bump must be the same with and
    without size_densities.
    """
    bnd_signal = SVsignal(
        ref_start=1000,
        ref_end=1000,
        read_start=0,
        read_end=0,
        size=300,
        sv_type=3,  # BNDL
    )
    ins_signal = _ins(ref_start=200, ref_end=300, size=100)
    ava_signals = {"readA": [ins_signal, bnd_signal]}

    size_densities = sv_size_densities_from_ava_signals(ava_signals)

    unscaled = importance_densities_from_ava_signals(ava_signals, size_densities=None)
    scaled = importance_densities_from_ava_signals(
        ava_signals, size_densities=size_densities
    )

    bnd_peak_unscaled = unscaled["readA"][980:1020].max()
    bnd_peak_scaled = scaled["readA"][980:1020].max()

    np.testing.assert_allclose(bnd_peak_unscaled, bnd_peak_scaled, rtol=1e-6)
