"""Shared helper utilities for SVcomposite objects.

Kept in a separate module to avoid circular imports between
multisample_sv_calling and svcomposite_merging.
"""

import numpy as np

from .SVcomposite import SVcomposite


def _crIDs_from_svcomposite(svc: SVcomposite) -> set[int]:
    """Extract crIDs from an SVcomposite. The crID is the integer prefix of the consensusID (format: crID.subID)."""
    crIDs: set[int] = set()
    for svp in svc.svPatterns:
        try:
            crIDs.add(int(svp.consensusID.split(".")[0]))
        except (ValueError, IndexError):
            pass
    return crIDs


def _regions_str_from_svcomposite(svc: SVcomposite) -> str:
    """Get a compact string representation of the regions in an SVcomposite."""
    try:
        regions = svc.get_regions(tolerance_radius=1)
        if not regions:
            return "no_regions"
        unique_regions = list({(r[0], r[1], r[2]) for r in regions})
        unique_regions.sort()
        return ";".join(f"{r[0]}:{r[1]}-{r[2]}" for r in unique_regions[:5])
    except Exception:
        return "region_error"


def _svcomposite_log_id(svc: SVcomposite) -> str:
    """Get a concise identifier string for an SVcomposite for logging.
    Format: sv_type|size=N|crIDs={...}|regions=chr:start-end;...|consensusIDs=[...]|representative=...
    """
    sv_type = svc.sv_type.get_sv_type() if svc.sv_type else "UNKNOWN"
    crIDs = sorted(_crIDs_from_svcomposite(svc))
    consensusIDs = sorted({svp.samplenamed_consensusID for svp in svc.svPatterns})
    regions_str = _regions_str_from_svcomposite(svc)
    return f"{sv_type}|size={svc.get_size()}|crIDs={{{','.join(map(str, crIDs))}}}|regions={regions_str}|consensusIDs={consensusIDs}|representative={svc.get_representative_SVpattern()._log_id()}"


def _svcomposite_short_id(svc: SVcomposite) -> str:
    """Shorter identifier for pairwise log messages."""
    crIDs = sorted(_crIDs_from_svcomposite(svc))
    cids = [svp.samplenamed_consensusID for svp in svc.svPatterns]
    return f"crIDs={{{','.join(map(str, crIDs))}}}|cIDs={cids}"


def cohens_d(x: list | np.ndarray, y: list | np.ndarray) -> float | None:
    """
    Calculate Cohen's d effect size between two samples.
    Cohen's d = (mean1 - mean2) / pooled_standard_deviation
    Where pooled_standard_deviation = sqrt(((n1-1)*s1² + (n2-1)*s2²) / (n1+n2-2))

    Returns ``None`` when the effect size is **not estimable**, i.e. when both
    samples are constant -- a single observation counts as constant, and that
    case additionally has a pooled variance of 0/0 that cannot be formed at all.

    Cohen's *d* expresses a difference of means *in units of the within-group
    spread*. With no observed spread there are no such units, so there is no
    effect size to report, whatever the means do. This function used to
    fabricate one anyway: ``float("inf")`` when the means differed and ``0.0``
    when they did not. Both are wrong in the same way and only differ in which
    direction they mislead -- ``inf`` reads as "maximally separated" and ``0.0``
    as "indistinguishable", and neither is supported by the data. This branch is
    what let a genome-wide failure of the read-derived noise model (every
    distortion value pinned at 0.0, hence every population constant) pass
    through an entire benchmark campaign unnoticed -- and on real data it did
    not even produce the `inf` that would at least have looked odd, but a finite
    value of order 1e16, for the reason given below.

    ``None`` rather than ``nan`` is deliberate. ``nan`` compares False against
    every threshold, so a caller that forgets to check keeps running on a
    plausible-looking answer -- the same failure mode in a new disguise; and
    ``nan`` is already in use by ``svcomposite_merging._similar_size`` to mean
    "not computed", which is a different statement. ``None`` cannot be compared
    with ``<=`` at all, so no caller can ignore it by accident.

    What a non-estimable effect size *means* is left to the caller, and the two
    callers in this codebase answer differently -- see
    ``svcomposite_merging._similar_size`` (refuse to merge on it) and
    ``candidateregions.signalstrength_to_crs`` (fall back to whether the two
    constant groups coincide). That divergence is the reason the policy is not
    baked in here.
    """
    if len(x) == 0 or len(y) == 0:
        raise ValueError("Both samples must contain at least one value")

    nx = len(x)
    ny = len(y)

    x_arr = np.array(x) if not isinstance(x, np.ndarray) else x
    y_arr = np.array(y) if not isinstance(y, np.ndarray) else y

    # Degeneracy is a property of the INPUTS and has to be decided on the inputs.
    # Testing the computed pooled standard deviation against 0.0 does not work,
    # which is how the old `inf` guard came to be dead code in production: the
    # mean of an array of *identical* elements is not exactly that element,
    # because numpy sums pairwise and the value need not be representable. A
    # genuinely constant sample of 19 copies of 6164.5678 has ptp exactly 0 but
    # `np.std(..., ddof=1) == 9.3e-13`, so `pooled_std == 0` is False and the
    # quotient comes out around 1e16 -- a finite, plausible-looking effect size
    # rather than the `inf` the guard was watching for. `np.ptp` is exact for
    # every magnitude and length, so it is what decides here.
    #
    # A single observation is a constant sample too, which folds in the old
    # `nx == 1 and ny == 1` early return; that one has to stay guarded anyway,
    # since the pooled variance below would be 0/0.
    x_constant = nx == 1 or bool(np.ptp(x_arr) == 0)
    y_constant = ny == 1 or bool(np.ptp(y_arr) == 0)
    if x_constant and y_constant:
        return None

    # Calculate means
    mean_x = np.mean(x_arr)
    mean_y = np.mean(y_arr)

    # Calculate sample standard deviations (with Bessel's correction, ddof=1)
    std_x = np.std(x_arr, ddof=1) if nx > 1 else 0.0
    std_y = np.std(y_arr, ddof=1) if ny > 1 else 0.0

    # Calculate pooled standard deviation
    pooled_var = ((nx - 1) * std_x**2 + (ny - 1) * std_y**2) / (nx + ny - 2)
    pooled_std = np.sqrt(pooled_var)

    # Belt and braces: only reachable if a non-constant sample's variance
    # underflows, which needs a spread at the very bottom of the subnormal
    # range. Still not an effect size.
    if pooled_std == 0:
        return None

    # Calculate Cohen's d
    cohens_d_value = (mean_x - mean_y) / pooled_std

    return float(cohens_d_value)
