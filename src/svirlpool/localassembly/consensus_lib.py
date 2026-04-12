from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import matplotlib
import numpy as np
import pysam
from Bio import SeqIO

matplotlib.use("Agg")  # non-interactive backend
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from ..signalprocessing import alignments_to_rafs
from ..util.datatypes import SVsignal

logging.getLogger("matplotlib").setLevel(logging.WARNING)
log = logging.getLogger(__name__)


def size_damping(size: float, s0: float = 30.0, alpha: float = 1.5) -> float:
    """Monotonic size-dependent damping that suppresses small indel signals.

    Returns a value in [0, 1). Signals much larger than *s0* contribute ~1.0;
    signals near or below *s0* are damped toward 0.

    f(s) = 1 - exp(-(s / s0) ** alpha)

    Args:
        size: Absolute signal size in bp.
        s0: Transition size (bp) below which signals are strongly damped.
        alpha: Sharpness of the transition (higher = sharper).
    """
    return 1.0 - np.exp(-((size / s0) ** alpha))


def parse_sv_signals_from_ava_alignments(
    ava_alignments: Path,
    min_signal_size: int,
    min_bnd_size: int,
) -> Dict[str, list[SVsignal]]:
    """Parse SV signals from an all-vs-all alignment SAM file.

    For each alignment where read A is the query and read B is the reference,
    SV signals are extracted and collected under the reference read (B).

    Returns a dict mapping each reference read name to all SVsignals observed
    across every query read that aligned to it.
    """
    result: Dict[str, list[SVsignal]] = {}

    with pysam.AlignmentFile(str(ava_alignments), "r", check_sq=False) as sam:
        for aln in sam.fetch(until_eof=True):
            if aln.is_unmapped or aln.reference_name is None:
                continue
            if aln.cigartuples is None:
                continue

            ref_name = aln.reference_name
            ref_start = aln.reference_start
            ref_end = aln.reference_end
            read_start = aln.query_alignment_start
            read_end = aln.query_alignment_end

            if (
                ref_start is None
                or ref_end is None
                or read_start is None
                or read_end is None
            ):
                raise ValueError(
                    f"parse_directed_signals_from_ava::Alignment missing required positions: ref_start={ref_start}, ref_end={ref_end}, read_start={read_start}, read_end={read_end}"
                )

            signals = alignments_to_rafs.parse_SVsignals_from_alignment(
                alignment=aln,
                ref_start=ref_start,
                ref_end=ref_end,
                read_start=read_start,
                read_end=read_end,
                min_signal_size=min_signal_size,
                min_bnd_size=min_bnd_size,
            )

            if signals:
                result.setdefault(ref_name, []).extend(signals)
    return result


def importance_densities_from_ava_signals(
    ava_signals: Dict[str, list[SVsignal]],
    size_densities: Dict[int, Dict[str, np.ndarray]] | None = None,
    window_size: int = 100,
    falloff: float = 0.5,
    damping_s0: float = 30.0,
    damping_alpha: float = 1.5,
    density_power: float = 1.5,
) -> Dict[str, np.ndarray]:
    """Build position-space importance density tracks, optionally scaled by SV size density.

    Each signal contributes a sigmoid-shaped plateau centred at its reference
    position.  Contributions are multiplicatively scaled by:

    1. ``size_damping(abs(signal.size), damping_s0, damping_alpha)`` — a
       monotonic function that suppresses small indel signals typical of
       Nanopore sequencing noise.
    2. When *size_densities* is provided, the normalised size-density at the
       signal's absolute size, raised to *density_power* (>1 makes the
       damping of rare sizes more aggressive).

    BND signals (sv_type 3/4) are excluded from the density tracks because
    their sizes reflect clipping length rather than SV content.

    Args:
        ava_signals: Output of ``parse_sv_signals_from_ava_alignments``.
        size_densities: Optional output of ``sv_size_densities_from_ava_signals``.
        window_size: Half-plateau width (bp) around each signal centre.
        falloff: Sigmoid steepness at plateau edges (higher = sharper).
        damping_s0: Transition size (bp) for ``size_damping``.
        damping_alpha: Sharpness exponent for ``size_damping``.
        density_power: Exponent applied to normalised size-density scaling.

    Returns:
        Dict mapping reference read name to its positional importance density array.
    """
    densities: Dict[str, np.ndarray] = {}

    for ref_name, signals in ava_signals.items():
        if not signals:
            continue

        ref_length = max(s.ref_end for s in signals)
        if ref_length == 0:
            continue

        # Pre-compute normalised size-density lookup tables for this read
        size_scale_cache: Dict[int, np.ndarray] = {}  # sv_type -> normalised array
        if size_densities is not None:
            for sv_type in (0, 1):
                arr = size_densities[sv_type].get(ref_name)
                if arr is not None and arr.max() > 0:
                    size_scale_cache[sv_type] = (arr / arr.max()) ** density_power

        density = np.zeros(ref_length, dtype=np.float64)
        half_win = window_size / 2.0

        for sig in signals:
            # Skip BND signals — they don't contribute to positional importance
            if sig.sv_type in (3, 4):
                continue

            center = (sig.ref_start + sig.ref_end) / 2.0
            lo = int(max(0, center - half_win - 6.0 / falloff))
            hi = int(min(ref_length, center + half_win + 6.0 / falloff + 1))
            if lo >= hi:
                continue
            x = np.arange(lo, hi, dtype=np.float64)
            left = 1.0 / (1.0 + np.exp(-falloff * (x - (center - half_win))))
            right = 1.0 / (1.0 + np.exp(falloff * (x - (center + half_win))))
            bump = left * right

            # Size-damping: suppress small indels
            damp = size_damping(abs(sig.size), s0=damping_s0, alpha=damping_alpha)

            # Scale by size density for INS/DEL
            scale = 1.0
            if sig.sv_type in (0, 1) and sig.sv_type in size_scale_cache:
                scale_arr = size_scale_cache[sig.sv_type]
                size_idx = abs(sig.size)
                if size_idx < len(scale_arr):
                    scale = float(scale_arr[size_idx])

            density[lo:hi] += damp * scale * bump

        densities[ref_name] = density

    return densities


def sv_size_densities_from_ava_signals(
    ava_signals: Dict[str, list[SVsignal]],
    falloff: float = 0.5,
    figures_dir: Path | None = None,
) -> Dict[int, Dict[str, np.ndarray]]:
    """Build size-space density tracks for INS and DEL signals.

    For each reference read, creates a 1-D density array indexed by SV size.
    Each INS (sv_type=0) or DEL (sv_type=1) signal contributes a smooth plateau
    bump centred at abs(signal.size), with a half-window of 20% of that size
    (i.e., the plateau spans from 0.8*size to 1.2*size).
    BND signals are excluded because their sizes are not informative.

    The sigmoid transition sharpness at the plateau edges is controlled by
    `falloff` (higher = sharper).

    Returns:
        {sv_type: {ref_read_name: density_array}}
        where sv_type is 0 (INS) or 1 (DEL).
        density_array[k] is the accumulated density at SV size k.
        Query with: size_densities[sv_signal.sv_type][abs(sv_signal.size)]
    """
    result: Dict[int, Dict[str, np.ndarray]] = {0: {}, 1: {}}

    for ref_name, signals in ava_signals.items():
        for sv_type in (0, 1):
            type_signals = [s for s in signals if s.sv_type == sv_type]
            if not type_signals:
                continue

            max_size = max(abs(s.size) for s in type_signals)
            if max_size == 0:
                continue

            density = np.zeros(max_size + 1, dtype=np.float64)

            for sig in sorted(type_signals, key=lambda s: abs(s.size)):
                center = float(abs(sig.size))
                half_win = 0.2 * center
                lo = int(max(0, center - half_win - 6.0 / falloff))
                hi = int(min(max_size + 1, center + half_win + 6.0 / falloff + 1))
                if lo >= hi:
                    continue
                x = np.arange(lo, hi, dtype=np.float64)
                # smooth plateau: rising sigmoid at (center - half_win), falling at (center + half_win)
                left = 1.0 / (1.0 + np.exp(-falloff * (x - (center - half_win))))
                right = 1.0 / (1.0 + np.exp(falloff * (x - (center + half_win))))
                density[lo:hi] += left * right

            result[sv_type][ref_name] = density

    if figures_dir is not None:
        sv_type_names = {0: "INS", 1: "DEL"}
        for sv_type, label in sv_type_names.items():
            type_tracks = result[sv_type]
            if not type_tracks:
                continue
            read_names = sorted(type_tracks.keys())
            fig, ax = plt.subplots(
                figsize=(
                    max(10, max(len(d) for d in type_tracks.values()) // 100),
                    max(4, len(read_names) * 0.6),
                )
            )
            for name in read_names:
                track = type_tracks[name]
                ax.plot(np.arange(len(track)), track, label=name, linewidth=1)
            ax.set_xlabel("SV size (bp)")
            ax.set_ylabel("Size density")
            ax.set_xscale("log")
            ax.set_title(f"{label} size density per reference read")
            ax.legend(fontsize=7, loc="upper right")
            fig.tight_layout()
            out = figures_dir / f"sv_size_densities_{label.lower()}.png"
            fig.savefig(str(out), dpi=150)
            plt.close(fig)
            log.info(f"Saved {label} size-density figure to {out}")

    return result


# ---------------------------------------------------------------------------
# Bidirectional AVA signal parsing (keeps query→ref direction)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DirectedSignals:
    """SV signals from one alignment direction (query aligned to reference)."""

    query_name: str
    ref_name: str
    signals: list[SVsignal]
    ref_length: int


def parse_directed_signals_from_ava(
    ava_alignments: Path,
    min_signal_size: int,
    min_bnd_size: int,
) -> list[DirectedSignals]:
    """Parse AVA alignments keeping the query→ref direction for each alignment.

    Unlike ``parse_sv_signals_from_ava_alignments`` which groups signals only
    by reference read, this function preserves which query produced which
    signals on which reference, needed for asymmetric pairwise distance
    computation.

    Returns a list of ``DirectedSignals`` objects, one per mapped alignment.
    """
    result: list[DirectedSignals] = []

    with pysam.AlignmentFile(str(ava_alignments), "r", check_sq=False) as sam:
        for aln in sam.fetch(until_eof=True):
            if aln.is_unmapped or aln.reference_name is None:
                continue
            if aln.cigartuples is None:
                continue

            ref_name = aln.reference_name
            query_name = aln.query_name
            ref_start = aln.reference_start
            ref_end = aln.reference_end
            read_start = aln.query_alignment_start
            read_end = aln.query_alignment_end

            if query_name is None:
                raise ValueError(
                    "parse_directed_signals_from_ava::Alignment query name is missing"
                )
            if (
                ref_end is None
                or ref_start is None
                or read_start is None
                or read_end is None
            ):
                raise ValueError(
                    f"parse_directed_signals_from_ava::Alignment missing required positions: ref_start={ref_start}, ref_end={ref_end}, read_start={read_start}, read_end={read_end}"
                )

            # Use inferred reference length from the alignment.
            # For the ref read in AVA this equals the read's own length.
            ref_length = aln.reference_length if aln.reference_length else 0

            signals = alignments_to_rafs.parse_SVsignals_from_alignment(
                alignment=aln,
                ref_start=ref_start,
                ref_end=ref_end,
                read_start=read_start,
                read_end=read_end,
                min_signal_size=min_signal_size,
                min_bnd_size=min_bnd_size,
            )

            result.append(
                DirectedSignals(
                    query_name=query_name,
                    ref_name=ref_name,
                    signals=signals,
                    ref_length=ref_length,
                )
            )
    return result


# ---------------------------------------------------------------------------
# Pairwise similarity matrix
# ---------------------------------------------------------------------------


def _directed_distance(
    signals: list[SVsignal],
    ref_density: np.ndarray | None,
    query_density: np.ndarray | None,
    damping_s0: float,
    damping_alpha: float,
    densities_weight: float = 1.0,
) -> float:
    """Compute the directed distance contribution from one alignment's INS/DEL signals.

    Each signal is weighted by:
      w_s = imp_ref[ref_center] * imp_query[query_center] * size_damping(|size|) * |size|
    The absolute signal size is included so that the distance scales with
    the amount of structural variation (in bp), not just the signal count.
    BND signals are excluded (handled separately as a posterior adjustment).

    ``densities_weight`` controls how much the importance densities influence
    the weights: 0.0 → flat weights of 1.0 (densities ignored); 1.0 → current
    behaviour; >1.0 → densities amplified beyond their raw value.
    The weight is applied as a linear blend:
      w = 1.0 + densities_weight * (density_value - 1.0)
    """
    total = 0.0
    for sig in signals:
        if sig.sv_type in (3, 4):
            continue

        damp = size_damping(abs(sig.size), s0=damping_s0, alpha=damping_alpha)

        ref_center = (sig.ref_start + sig.ref_end) // 2
        query_center = (sig.read_start + sig.read_end) // 2

        w_ref = 1.0
        if ref_density is not None and 0 <= ref_center < len(ref_density):
            w_ref = 1.0 + densities_weight * (float(ref_density[ref_center]) - 1.0)

        w_query = 1.0
        if query_density is not None and 0 <= query_center < len(query_density):
            w_query = 1.0 + densities_weight * (
                float(query_density[query_center]) - 1.0
            )

        total += w_ref * w_query * damp * abs(sig.size)
    return total


def _has_interior_bnd(
    signals: list[SVsignal],
    ref_length: int,
    bnd_margin: int = 100,
) -> bool:
    """Return True if any BND signal is more than *bnd_margin* bp from the
    start or end of the reference read, indicating a true structural breakpoint.
    """
    for sig in signals:
        if sig.sv_type not in (3, 4):
            continue
        pos = sig.ref_start  # BNDs have ref_start == ref_end
        if bnd_margin < pos < ref_length - bnd_margin:
            return True
    return False


def pairwise_similarity_matrix(
    directed_signals: list[DirectedSignals],
    densities: Dict[str, np.ndarray],
    read_lengths: Dict[str, int],
    all_read_names: list[str] | None = None,
    length_penalty_lambda: float | None = None,
    damping_s0: float = 30.0,
    damping_alpha: float = 1.5,
    bnd_margin: int = 100,
    densities_weight: float = 1.0,
) -> Tuple[np.ndarray, list[str]]:
    """Compute a symmetric pairwise similarity matrix from AVA signals.

    The underlying distance for a pair (A, B) is::

        D(A,B) = max(d(A→B), d(B→A)) + λ * (1 − size_sim(A,B))

    where d(A→B) sums importance-weighted, size-damped, size-scaled signal
    contributions.  The *max* (rather than average) of the two directions is
    used so that a signal-free reverse alignment never dilutes a clear signal
    from the forward direction.
    Similarity is then computed via a Gaussian kernel::

        sim(A,B) = exp(-D² / (2σ²))

    where σ is auto-calibrated to the median of finite non-zero distances.

    Posterior adjustments:
    - If any alignment between A and B contains a BND more than *bnd_margin*
      bp from either end of the reference read, the similarity is set to 0.
    - If no alignment exists between A and B in either direction, the
      similarity is 0.
    - The diagonal is 1.

    Args:
        directed_signals: Output of ``parse_directed_signals_from_ava``.
        densities: Importance density tracks from ``importance_densities_from_ava_signals``.
        read_lengths: Read name → length mapping.
        all_read_names: Ordered list of read names to include.  When *None*,
            derived from union of all reads seen in *directed_signals* and
            *read_lengths*.
        length_penalty_lambda: Additive penalty weight for length dissimilarity.
            When *None* (default), auto-calibrated to the median of the
            signal-based directed distances.
        damping_s0: Transition size for ``size_damping``.
        damping_alpha: Sharpness exponent for ``size_damping``.
        bnd_margin: Minimum distance (bp) from reference start/end for a BND
            to be considered interior (and trigger similarity = 0).
        densities_weight: Scaling factor for the importance density weights.
            0.0 → densities have no effect (all weights = 1.0, plain signals);
            1.0 → current behaviour (raw density values used as weights);
            >1.0 → density influence amplified beyond raw values.
        max_intra_distance: Hard upper bound on the signal-only pairwise distance.
            When > 0, any pair whose signal-only distance (before the length
            penalty) exceeds this value has its similarity forced to 0.
            -1.0 (default) disables the threshold.

    Returns:
        (similarity_matrix, read_names) — a symmetric numpy array of shape
        (n, n) with values in [0, 1], and the ordered list of read names.
    """
    # Determine the set of reads
    if all_read_names is not None:
        read_names = list(all_read_names)
    else:
        names: set[str] = set(read_lengths.keys())
        for ds in directed_signals:
            names.add(ds.query_name)
            names.add(ds.ref_name)
        read_names = sorted(names)

    n = len(read_names)
    idx = {name: i for i, name in enumerate(read_names)}

    # Accumulate directed distances and track alignment presence + interior BNDs
    # d_forward[i][j] collects distances from alignment query=i, ref=j
    d_sum = np.zeros((n, n), dtype=np.float64)
    d_count = np.zeros((n, n), dtype=np.int32)
    has_alignment = np.zeros((n, n), dtype=bool)
    has_interior_bnd = np.zeros((n, n), dtype=bool)

    for ds in directed_signals:
        qi = idx.get(ds.query_name)
        ri = idx.get(ds.ref_name)
        if qi is None or ri is None:
            continue
        has_alignment[qi, ri] = True

        ref_density = densities.get(ds.ref_name)
        query_density = densities.get(ds.query_name)

        dd = _directed_distance(
            ds.signals,
            ref_density,
            query_density,
            damping_s0=damping_s0,
            damping_alpha=damping_alpha,
            densities_weight=densities_weight,
        )
        d_sum[qi, ri] += dd
        d_count[qi, ri] += 1

        if _has_interior_bnd(ds.signals, ds.ref_length, bnd_margin=bnd_margin):
            has_interior_bnd[qi, ri] = True

    # Build symmetric distance matrix
    distance = np.full((n, n), np.inf, dtype=np.float64)
    for i in range(n):
        distance[i, i] = 0.0
        for j in range(i + 1, n):
            # Check alignment existence in either direction
            has_ij = has_alignment[i, j] or has_alignment[j, i]
            if not has_ij:
                # No alignment → infinite distance (similarity will be 0)
                continue

            # Directed distances: take the max so that a signal-free
            # reverse alignment never dilutes a clear forward signal.
            d_ij = d_sum[i, j] / max(d_count[i, j], 1) if has_alignment[i, j] else 0.0
            d_ji = d_sum[j, i] / max(d_count[j, i], 1) if has_alignment[j, i] else 0.0
            pair_d = max(d_ij, d_ji)

            distance[i, j] = pair_d
            distance[j, i] = pair_d

    # Compute sigma from signal-only distances (before length penalty).
    # This ensures that the length penalty shifts similarities relative to
    # a fixed kernel width derived from the SV-signal content alone.
    finite_signal_dists = distance[np.isfinite(distance) & (distance > 0)]
    sigma = (
        float(np.median(finite_signal_dists)) if finite_signal_dists.size > 0 else 1.0
    )
    if sigma == 0:
        sigma = 1.0

    # Add length-dissimilarity penalty
    lengths_arr = np.array(
        [float(read_lengths.get(r, 0)) for r in read_names], dtype=np.float64
    )
    # Auto-calibrate lambda if not provided: scale to the median signal distance
    if length_penalty_lambda is None:
        length_penalty_lambda = sigma

    for i in range(n):
        for j in range(i + 1, n):
            if not np.isfinite(distance[i, j]):
                continue
            lo = min(lengths_arr[i], lengths_arr[j])
            hi = max(lengths_arr[i], lengths_arr[j])
            size_sim = lo / hi if hi > 0 else 0.0
            penalty = length_penalty_lambda * (1.0 - size_sim)
            distance[i, j] += penalty
            distance[j, i] += penalty

    # Convert distance to similarity using a Gaussian kernel:
    #   sim = exp(-D^2 / (2 * sigma^2))
    # sigma was fixed above from signal-only distances so that the length
    # penalty genuinely reduces similarity for length-mismatched pairs.

    similarity = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        similarity[i, i] = 1.0
        for j in range(i + 1, n):
            if np.isfinite(distance[i, j]):
                sim = float(np.exp(-(distance[i, j] ** 2) / (2.0 * sigma**2)))
                similarity[i, j] = sim
                similarity[j, i] = sim
            # else: stays 0 (no alignment)

    # Posterior BND adjustment: if interior BND in either direction → similarity = 0
    for i in range(n):
        for j in range(i + 1, n):
            if has_interior_bnd[i, j] or has_interior_bnd[j, i]:
                similarity[i, j] = 0.0
                similarity[j, i] = 0.0

    return similarity, read_names


def visualize_raw_signal_matrix(
    directed_signals: list["DirectedSignals"],
    all_read_names: list[str],
    output: Path,
) -> None:
    """Visualise raw pairwise signal sums from all-vs-all alignments.

    Each cell (i, j) shows the combined INS + DEL size sum between reads i and
    j (both alignment directions summed).  Cells where any BND signal is
    present (in either direction) are drawn in grey instead of using colour
    scaling.

    Args:
        directed_signals: Output of ``parse_directed_signals_from_ava``.
        all_read_names: Ordered list of read names.
        output: Path where the PNG figure is saved.
    """
    read_names = list(all_read_names)
    n = len(read_names)
    if n == 0:
        log.warning("No reads – nothing to visualise.")
        return

    idx = {name: i for i, name in enumerate(read_names)}

    signal_sums = np.zeros((n, n), dtype=np.float64)
    has_bnd = np.zeros((n, n), dtype=bool)

    for ds in directed_signals:
        qi = idx.get(ds.query_name)
        ri = idx.get(ds.ref_name)
        if qi is None or ri is None:
            continue
        for sig in ds.signals:
            if sig.sv_type in (3, 4):  # BNDL / BNDR
                has_bnd[qi, ri] = True
            elif sig.sv_type in (0, 1, 2):  # INS, DEL, DELR
                signal_sums[qi, ri] += abs(sig.size)

    # Make symmetric: [i,j] = sum from i→j direction + sum from j→i direction
    signal_sums = signal_sums + signal_sums.T
    has_bnd = has_bnd | has_bnd.T
    np.fill_diagonal(signal_sums, 0.0)
    np.fill_diagonal(has_bnd, False)

    fig, ax = plt.subplots(figsize=(max(8, n), max(7, n - 1)))

    # Build RGBA image: YlOrRd for signal values, grey for BND cells
    cmap = plt.get_cmap("YlOrRd")
    non_bnd_values = signal_sums[~has_bnd & (np.arange(n)[:, None] != np.arange(n))]
    vmax = (
        float(non_bnd_values.max())
        if non_bnd_values.size > 0 and non_bnd_values.max() > 0
        else 1.0
    )
    norm = mcolors.Normalize(vmin=0, vmax=vmax)

    rgba = cmap(norm(signal_sums))
    rgba[has_bnd] = [0.65, 0.65, 0.65, 1.0]  # grey for BND cells

    ax.imshow(rgba, aspect="auto", interpolation="nearest")

    fontsize = max(5, 9 - n // 5)
    for i in range(n):
        for j in range(n):
            val = signal_sums[i, j]
            if has_bnd[i, j]:
                label = f"{val:.0f}\nBND"
                ax.text(
                    j,
                    i,
                    label,
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    color="white",
                    linespacing=1.2,
                )
            else:
                text_color = "black" if norm(val) < 0.65 else "white"
                ax.text(
                    j,
                    i,
                    f"{val:.0f}",
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    color=text_color,
                )

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(read_names, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(read_names, fontsize=8)
    ax.set_xlabel("Read B")
    ax.set_ylabel("Read A")
    ax.set_title("Raw pairwise signal sums (INS + DEL sizes; grey = BND present)")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)  # type: ignore[attr-defined]
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Sum of INS + DEL sizes")

    fig.tight_layout()
    fig.savefig(str(output), dpi=150)
    log.info(f"Saved raw-signal matrix figure to {output}")
    plt.close(fig)


def visualize_normalized_similarity_matrix(
    similarity_matrix: np.ndarray,
    output: Path,
    read_names: list[str] | None = None,
) -> None:
    """Visualise the pairwise similarity matrix as an annotated heatmap.

    Values are expected to be in [0, 1] (already normalised).

    Args:
        similarity_matrix: Symmetric (n, n) array with values in [0, 1].
        output: Path where the PNG figure is saved.
        read_names: Ordered list of read names for axis labels.
            When None, integer indices are used.
    """
    n = similarity_matrix.shape[0]
    if n == 0:
        log.warning("Empty similarity matrix – nothing to visualise.")
        return

    labels = read_names if read_names is not None else [str(i) for i in range(n)]

    fig, ax = plt.subplots(figsize=(max(8, n), max(7, n - 1)))
    im = ax.imshow(
        similarity_matrix,
        cmap="viridis",
        vmin=0,
        vmax=1,
        aspect="auto",
        interpolation="nearest",
    )

    for i in range(n):
        for j in range(n):
            val = similarity_matrix[i, j]
            text_color = "white" if val < 0.6 else "black"
            ax.text(
                j,
                i,
                f"{val:.2f}",
                ha="center",
                va="center",
                fontsize=max(6, 10 - n // 5),
                color=text_color,
            )

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Read B")
    ax.set_ylabel("Read A")
    ax.set_title("Pairwise similarity matrix")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Similarity")

    fig.tight_layout()
    fig.savefig(str(output), dpi=150)
    log.info(f"Saved similarity-matrix figure to {output}")


# Palette of visually distinct colours for cluster nodes
_CLUSTER_COLORS = [
    "#e41a1c",
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#ff7f00",
    "#a65628",
    "#f781bf",
    "#999999",
]


def visualize_similarity_graph(
    similarity_matrix: np.ndarray,
    read_names: list[str],
    output: Path,
    node_colors: dict[str, str] | None = None,
    edge_threshold: float = 0.1,
) -> None:
    """Visualise the pairwise similarity as a weighted graph.

    Nodes represent reads; edges are drawn for pairs with similarity at or
    above *edge_threshold*.  Edge colour **and** width encode the similarity
    weight (thicker / warmer-coloured = more similar).

    Two files are written:

    * A static PNG at *output* (matplotlib).
    * An interactive HTML at ``output.with_suffix(".html")`` (gravis / d3.js).
      Nodes are draggable and the layout uses a force simulation.  Clicking a
      node shows its full read name in the details panel where it can be
      selected and copied into IGV.

    Layout is initialised via MDS on the dissimilarity matrix so that closely
    related reads appear near each other.  Node positions can be freely
    adjusted in the browser by dragging.

    Args:
        similarity_matrix: Symmetric (n, n) array with values in [0, 1].
        read_names: Ordered list of read names matching matrix rows/columns.
        output: Path for the PNG file.
        node_colors: Optional mapping of read name → CSS colour string.  When
            provided each node is drawn in the given colour instead of the
            default steelblue.  Outliers should be passed as ``"#aaaaaa"``.
        edge_threshold: Minimum similarity to draw an edge (default 0.1).
    """
    import gravis as gv
    from sklearn.manifold import MDS

    n = len(read_names)
    if n == 0:
        log.warning("Empty read list – nothing to visualise in similarity graph.")
        return

    # ------------------------------------------------------------------
    # 2-D layout via MDS on the dissimilarity matrix
    # ------------------------------------------------------------------
    dissimilarity = np.clip(1.0 - similarity_matrix, 0.0, 1.0)
    np.fill_diagonal(dissimilarity, 0.0)

    if n == 1:
        pos = np.array([[0.0, 0.0]])
    elif n == 2:
        pos = np.array([[-0.5, 0.0], [0.5, 0.0]])
    else:
        mds = MDS(
            n_components=2,
            dissimilarity="precomputed",
            random_state=42,
            normalized_stress="auto",
        )
        pos = mds.fit_transform(dissimilarity)

    # ------------------------------------------------------------------
    # Collect edges above threshold
    # ------------------------------------------------------------------
    edge_pairs: list[tuple[int, int, float]] = []
    for i in range(n):
        for j in range(i + 1, n):
            w = float(similarity_matrix[i, j])
            if w >= edge_threshold:
                edge_pairs.append((i, j, w))

    # ------------------------------------------------------------------
    # matplotlib PNG
    # ------------------------------------------------------------------
    cmap = plt.get_cmap("plasma")
    norm = mcolors.Normalize(vmin=0.0, vmax=1.0)

    fig_mpl, ax = plt.subplots(figsize=(max(8, n), max(8, n)))

    for i, j, w in edge_pairs:
        ax.plot(
            [pos[i, 0], pos[j, 0]],
            [pos[i, 1], pos[j, 1]],
            color=cmap(norm(w)),
            linewidth=1.0 + w * 5.0,
            alpha=0.7,
            zorder=1,
        )
        mx = (pos[i, 0] + pos[j, 0]) / 2
        my = (pos[i, 1] + pos[j, 1]) / 2
        ax.text(
            mx,
            my,
            f"{w:.2f}",
            fontsize=6,
            ha="center",
            va="center",
            color="dimgrey",
            zorder=2,
        )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig_mpl.colorbar(sm, ax=ax, label="Similarity", fraction=0.03, pad=0.02)

    for idx, name in enumerate(read_names):
        node_color = node_colors.get(name, "steelblue") if node_colors else "steelblue"
        ax.scatter(
            pos[idx, 0],
            pos[idx, 1],
            s=200,
            zorder=3,
            color=node_color,
            edgecolors="white",
            linewidths=1,
        )
        ax.text(
            pos[idx, 0],
            pos[idx, 1],
            name,
            fontsize=7,
            ha="center",
            va="bottom",
            zorder=4,
            bbox={
                "boxstyle": "round,pad=0.1",
                "fc": "white",
                "alpha": 0.6,
                "ec": "none",
            },
        )

    ax.set_title("Read similarity graph")
    ax.axis("off")
    fig_mpl.tight_layout()
    fig_mpl.savefig(str(output), dpi=150)
    plt.close(fig_mpl)
    log.info(f"Saved similarity-graph PNG to {output}")

    # ------------------------------------------------------------------
    # gravis interactive HTML via d3.js
    # No fixed positions: the d3 force layout runs freely so nodes can be
    # dragged into any arrangement.  Clicking a node opens a details panel
    # showing the full read name for easy copy into IGV.
    # ------------------------------------------------------------------

    # Build gJGF nodes – read name is used as node id and label
    # No fixed x/y so the d3 force layout runs freely and nodes can be dragged
    gjgf_nodes: dict = {}
    for name in read_names:
        node_color = node_colors.get(name, "steelblue") if node_colors else "steelblue"
        gjgf_nodes[name] = {
            "metadata": {
                "size": 18,
                "color": node_color,
                "label_color": "black",
                "hover": f"<b>{name}</b>",
                "click": name,  # shown in details panel → easy copy
            }
        }

    # Build gJGF edges
    gjgf_edges: list = []
    for i, j, w in edge_pairs:
        rgba = cmap(norm(w))
        hex_color = "#{:02x}{:02x}{:02x}".format(
            int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255)
        )
        gjgf_edges.append({
            "source": read_names[i],
            "target": read_names[j],
            "label": f"{w:.2f}",
            "metadata": {
                "size": float(max(0.5, w * 8.0)),
                "color": hex_color,
                "hover": f"similarity: {w:.2f}",
            },
        })

    gjgf = {
        "graph": {
            "directed": False,
            "metadata": {
                "background_color": "#f9f9f9",
                "node_label_size": 8,
                "edge_label_size": 7,
            },
            "nodes": gjgf_nodes,
            "edges": gjgf_edges,
        }
    }

    fig_gv = gv.d3(
        gjgf,
        graph_height=700,
        details_height=60,
        show_details=True,
        show_menu=True,
        show_node_label=True,
        show_node_label_border=False,
        node_label_data_source="id",
        show_edge_label=True,
        show_edge_label_border=False,
        edge_label_data_source="label",
        use_edge_size_normalization=True,
        edge_size_normalization_min=0.5,
        edge_size_normalization_max=10.0,
        edge_size_data_source="size",
        node_hover_tooltip=True,
        edge_hover_tooltip=True,
        node_drag_fix=True,
        zoom_factor=0.9,
        links_force_distance=150,
    )

    html_output = output.with_suffix(".html")
    fig_gv.export_html(str(html_output))
    log.info(f"Saved interactive similarity-graph HTML to {html_output}")


# ---------------------------------------------------------------------------
# Outlier detection
# ---------------------------------------------------------------------------


def detect_outlier_reads(
    similarity: np.ndarray,
    read_names: list[str],
    read_lengths: Dict[str, int],
    n_clusters: int,
    length_factor: float = 1.2,
    min_connections: int | None = None,
    connectivity_threshold: float = 0.1,
) -> set[str]:
    n = len(read_names)
    if n == 0:
        return set()
    log.debug(
        f"detect_outlier_reads: n_reads={n}, n_clusters={n_clusters}, "
        f"length_factor={length_factor}, connectivity_threshold={connectivity_threshold}, "
        f"min_connections={min_connections}"
    )
    by_length = _outlier_reads_by_length(
        read_names=read_names,
        read_lengths=read_lengths,
        n_clusters=n_clusters,
        length_factor=length_factor,
    )
    by_similarity = _outlier_reads_by_similarity(
        similarity=similarity,
        read_names=read_names,
        n_clusters=n_clusters,
        min_connections=min_connections,
        connectivity_threshold=connectivity_threshold,
    )
    log.debug(
        f"detect_outlier_reads: length outliers={by_length}, similarity outliers={by_similarity}"
    )
    return by_length | by_similarity


def _outlier_reads_by_length(
    read_names: list[str],
    read_lengths: Dict[str, int],
    n_clusters: int,
    length_factor: float = 1.2,
) -> set[str]:
    """Flag reads whose length makes them true outliers.

    Reads are sorted by length and clustered by detecting gaps where the
    ratio between consecutive lengths exceeds *length_factor*.  Only reads
    in clusters that are too small to represent a genuine allele/haplotype
    are flagged as outliers.
    """
    n = len(read_names)
    if n == 0:
        return set()

    sorted_reads = sorted(read_names, key=lambda r: read_lengths.get(r, 0))

    # Build clusters by detecting gaps in sorted lengths
    clusters: list[list[str]] = [[sorted_reads[0]]]
    for i in range(1, n):
        prev_len = read_lengths.get(sorted_reads[i - 1], 0)
        curr_len = read_lengths.get(sorted_reads[i], 0)
        if prev_len > 0 and curr_len / prev_len > length_factor:
            clusters.append([])
        elif prev_len == 0 and curr_len > 0:
            clusters.append([])
        clusters[-1].append(sorted_reads[i])

    min_cluster_size = max(2, int(np.sqrt(n / max(n_clusters, 1))))
    log.debug(
        f"_outlier_reads_by_length: n={n}, n_clusters={n_clusters}, "
        f"length_factor={length_factor}, {len(clusters)} length-clusters found, "
        f"min_cluster_size={min_cluster_size}"
    )

    length_outliers: set[str] = set()
    for idx, cluster in enumerate(clusters):
        cluster_lengths = [read_lengths.get(r, 0) for r in cluster]
        log.debug(f"  cluster {idx}: size={len(cluster)}, lengths={cluster_lengths}")
        if len(cluster) < min_cluster_size:
            log.debug(
                f"    -> flagging {len(cluster)} read(s) as length outlier(s) "
                f"(cluster too small: {len(cluster)} < {min_cluster_size})"
            )
            length_outliers.update(cluster)

    return length_outliers


def _outlier_reads_by_similarity(
    similarity: np.ndarray,
    read_names: list[str],
    n_clusters: int,
    min_connections: int | None = None,
    connectivity_threshold: float = 0.1,
) -> set[str]:
    n = len(read_names)
    if min_connections is None:
        expected_cluster_size = int(np.ceil(n / max(n_clusters, 1)))
        cluster_based = max(1, int(0.5 * expected_cluster_size))
        # Also adapt to actual observed connectivity via the median
        conn_counts = np.array(
            [
                int(np.sum(similarity[i, :] > connectivity_threshold))
                - 1  # -1 for self
                for i in range(n)
            ],
            dtype=np.int64,
        )
        median_conn = float(np.median(conn_counts)) if n > 0 else 0.0
        median_based = max(1, int(0.5 * median_conn))
        min_connections = max(1, min(cluster_based, median_based))
    log.debug(
        f"_outlier_reads_by_similarity: n={n}, n_clusters={n_clusters}, "
        f"min_connections={min_connections}, connectivity_threshold={connectivity_threshold}"
    )
    idx = {name: i for i, name in enumerate(read_names)}
    connectivity_outliers: set[str] = set()
    for r in read_names:
        i = idx[r]
        connections = [
            read_names[j]
            for j in range(n)
            if j != i and similarity[i, j] > connectivity_threshold
        ]
        n_conn = len(connections)
        log.debug(
            f"  '{r}': {n_conn} connections >= {connectivity_threshold} "
            f"(need {min_connections}): {connections}"
        )
        if n_conn < min_connections:
            log.debug(f"    -> flagging '{r}' as similarity outlier")
            connectivity_outliers.add(r)
    return connectivity_outliers


def visualize_ava_alignments(ava_alignments: Path, output: Path) -> None:
    """Parse a SAM AVA file and save an alignment-length matrix plot."""
    lengths: dict[tuple[str, str], int] = {}
    reads: set[str] = set()

    with pysam.AlignmentFile(str(ava_alignments), "r", check_sq=False) as sam:
        for aln in sam.fetch(until_eof=True):
            if aln.is_unmapped or aln.reference_name is None:
                continue
            query = aln.query_name
            ref = aln.reference_name
            alen = aln.query_alignment_length or 0
            if query is None:
                raise ValueError(
                    "visualize_alignment_presence::Alignment query name is missing"
                )
            reads.update([query, ref])
            key = (query, ref)
            lengths[key] = lengths.get(key, 0) + alen

    read_names = sorted(reads)
    n = len(read_names)
    idx = {name: i for i, name in enumerate(read_names)}

    matrix = np.zeros((n, n), dtype=np.int64)
    for (q, r), alen in lengths.items():
        if q in idx and r in idx:
            matrix[idx[q], idx[r]] = alen

    max_len = matrix.max() if matrix.max() > 0 else 1
    norm_matrix = matrix / max_len

    fig, ax = plt.subplots(figsize=(max(8, n), max(7, n - 1)))
    im = ax.imshow(norm_matrix, cmap="viridis", vmin=0, vmax=1, aspect="auto")

    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            norm_val = norm_matrix[i, j]
            text_color = "white" if norm_val < 0.6 else "black"
            ax.text(
                j,
                i,
                str(val),
                ha="center",
                va="center",
                fontsize=max(6, 10 - n // 5),
                color=text_color,
            )

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(read_names, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(read_names, fontsize=8)
    ax.set_xlabel("Reference read (B)")
    ax.set_ylabel("Query read (A)")
    ax.set_title("All-vs-all alignment lengths")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Normalised alignment length")

    fig.tight_layout()
    fig.savefig(str(output), dpi=150)
    log.info(f"Saved figure to {output}")


def visualize_importance_densities(
    densities: Dict[str, np.ndarray],
    output: Path,
) -> None:
    """Visualise per-read importance density tracks as a heatmap.

    Each row is one reference read. The x-axis is the position along that read.
    Because track lengths differ, shorter tracks are zero-padded on the right.
    A single global colour scale (vmin=0, vmax=global_max) is used so that
    colours are comparable across reads.
    """
    if not densities:
        log.warning("No density tracks to visualise.")
        return

    read_names = sorted(densities.keys())
    max_len = max(d.size for d in densities.values())
    global_max = max(d.max() for d in densities.values() if d.size > 0)
    if global_max == 0:
        global_max = 1.0

    # Build a 2-D matrix (reads × positions), zero-padded
    matrix = np.zeros((len(read_names), max_len), dtype=np.float64)
    for row, name in enumerate(read_names):
        track = densities[name]
        matrix[row, : track.size] = track

    n_reads = len(read_names)
    fig_h = max(3, n_reads * 0.5 + 1.5)
    fig_w = max(10, max_len / 500)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(
        matrix,
        cmap="viridis",
        vmin=0,
        vmax=global_max,
        aspect="auto",
        interpolation="nearest",
    )

    ax.set_yticks(range(n_reads))
    ax.set_yticklabels(read_names, fontsize=8)
    ax.set_xlabel("Position along reference read (bp)")
    ax.set_title("Importance density tracks (unified scale)")

    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label(f"Importance density (max={global_max:.2f})")

    fig.tight_layout()
    fig.savefig(str(output), dpi=150)
    log.info(f"Saved importance-density figure to {output}")


def read_lengths_from_fasta(fasta: Path) -> Dict[str, int]:
    """Return a dict mapping each record ID to its sequence length."""
    return {record.id: len(record.seq) for record in SeqIO.parse(str(fasta), "fasta")}


def visualize_size_similarity(read_lengths: Dict[str, int], output: Path) -> None:
    """Visualise the reciprocal size similarity matrix of all reads.

    Cell (A, B) = min(len_A, len_B) / max(len_A, len_B), ranging 0–1.
    Values are annotated with the absolute shorter/longer ratio as a percentage.
    """
    read_names = sorted(read_lengths.keys())
    n = len(read_names)
    lengths = np.array([read_lengths[r] for r in read_names], dtype=np.float64)

    matrix = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            lo = min(lengths[i], lengths[j])
            hi = max(lengths[i], lengths[j])
            matrix[i, j] = lo / hi if hi > 0 else 0.0

    fig, ax = plt.subplots(figsize=(max(8, n), max(7, n - 1)))
    im = ax.imshow(matrix, cmap="viridis", vmin=0, vmax=1, aspect="auto")

    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            text_color = "white" if val < 0.6 else "black"
            ax.text(
                j,
                i,
                f"{val:.2f}",
                ha="center",
                va="center",
                fontsize=max(6, 10 - n // 5),
                color=text_color,
            )

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(read_names, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(read_names, fontsize=8)
    ax.set_xlabel("Read B")
    ax.set_ylabel("Read A")
    ax.set_title("Reciprocal size similarity (min/max read length)")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Size similarity")

    fig.tight_layout()
    fig.savefig(str(output), dpi=150)
    log.info(f"Saved size-similarity figure to {output}")


def visualize_alignment_presence(ava_alignments: Path, output: Path) -> None:
    """Visualise which reads are aligned to which as a boolean matrix.

    Cell (A, B) is blue when read A appears as a query aligned to read B as
    reference; grey otherwise. Row and column labels include the total number
    of alignments as query and as reference respectively.
    """
    aligned: set[tuple[str, str]] = set()
    reads: set[str] = set()
    aln_count: dict[str, int] = {}  # total alignments per read (as query or ref)

    with pysam.AlignmentFile(str(ava_alignments), "r", check_sq=False) as sam:
        for aln in sam.fetch(until_eof=True):
            if aln.is_unmapped or aln.reference_name is None:
                continue
            query = aln.query_name
            if query is None:
                raise (
                    ValueError(
                        "visualize_alignment_presence::Alignment query name is missing"
                    )
                )
            ref = aln.reference_name
            reads.update([query, ref])
            aligned.add((query, ref))
            aln_count[query] = aln_count.get(query, 0) + 1
            aln_count[ref] = aln_count.get(ref, 0) + 1

    read_names = sorted(reads)
    n = len(read_names)
    idx = {name: i for i, name in enumerate(read_names)}

    matrix = np.zeros((n, n), dtype=np.float64)
    for q, r in aligned:
        if q in idx and r in idx:
            matrix[idx[q], idx[r]] = 1.0

    # Labels include per-read total alignment count
    labels = [f"{name} ({aln_count.get(name, 0)})" for name in read_names]

    cmap = mcolors.ListedColormap(["#aaaaaa", "#1f77b4"])  # grey, blue

    fig, ax = plt.subplots(figsize=(max(8, n), max(7, n - 1)))
    ax.imshow(matrix, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Reference read (B)")
    ax.set_ylabel("Query read (A)")
    ax.set_title("Alignment presence (blue = aligned, grey = no alignment)")

    from matplotlib.patches import Patch

    ax.legend(
        handles=[
            Patch(facecolor="#1f77b4", label="Aligned"),
            Patch(facecolor="#aaaaaa", label="Not aligned"),
        ],
        loc="upper right",
        fontsize=8,
    )

    fig.tight_layout()
    fig.savefig(str(output), dpi=150)
    log.info(f"Saved alignment-presence figure to {output}")
