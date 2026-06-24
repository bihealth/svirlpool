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


def cohens_d(x: list | np.ndarray, y: list | np.ndarray) -> float:
    """
    Calculate Cohen's d effect size between two samples.
    Cohen's d = (mean1 - mean2) / pooled_standard_deviation
    Where pooled_standard_deviation = sqrt(((n1-1)*s1² + (n2-1)*s2²) / (n1+n2-2))
    """
    if len(x) == 0 or len(y) == 0:
        raise ValueError("Both samples must contain at least one value")

    nx = len(x)
    ny = len(y)

    x_arr = np.array(x) if not isinstance(x, np.ndarray) else x
    y_arr = np.array(y) if not isinstance(y, np.ndarray) else y

    # Calculate means
    mean_x = np.mean(x_arr)
    mean_y = np.mean(y_arr)

    # Handle edge cases
    if nx == 1 and ny == 1:
        # Can't calculate pooled standard deviation with only one observation each
        # Return a large effect size if means differ, 0 if they're the same
        return float("inf") if mean_x != mean_y else 0.0

    # Calculate sample standard deviations (with Bessel's correction, ddof=1)
    std_x = np.std(x_arr, ddof=1) if nx > 1 else 0.0
    std_y = np.std(y_arr, ddof=1) if ny > 1 else 0.0

    # Calculate pooled standard deviation
    pooled_var = ((nx - 1) * std_x**2 + (ny - 1) * std_y**2) / (nx + ny - 2)
    pooled_std = np.sqrt(pooled_var)

    # Handle case where pooled standard deviation is 0 (all values identical)
    if pooled_std == 0:
        return 0.0 if mean_x == mean_y else float("inf")

    # Calculate Cohen's d
    cohens_d_value = (mean_x - mean_y) / pooled_std

    return float(cohens_d_value)
