# Consensus stage: where the compute goes, and what was cut

Branch `perf/consensus-profiling`, on `exp/ava-read-phasing` d00f5ec. HG002 ONT
20x, the svp_improvements 5% genome (2004 containers, 30 batches of <= 100).

## Profile (before)

`profile_batch.py` runs `crs_containers_to_consensus` on real batches
under cProfile. It also attributes the time of every external tool, taken
from `os.wait4` rusage (which includes the tool's own children, e.g. lastal
under lamassemble), to the pipeline stage that launched it. Batch 25 was the
slowest batch of ava_phase2; batches 0 and 12 agree:

| stage (batch 25, 100 containers) | legacy | phased |
|---|---|---|
| total wall | 126 s | 338 s |
| all-vs-all minimap2 (`read_phasing.run_ava`) | -- | 190 s (56%) |
| lamassemble (`assemble_consensus`) | 95 s (76%) | 113 s (33%) |
| phasing in Python (`snv_sites`, `parse_pair`, `recurrent_sites`) | -- | 21 s (6%) |
| KMeans (`consensus_while_clustering_with_kmeans`) | 17 s wall, **361 s CPU** | 0.3 s |
| rest (read fetch, cutting, final_consensus, padding, ...) | ~13 s | ~14 s |

Over batches 0, 12 and 25: AVA ~60%, lamassemble ~28%, phasing Python ~7%.

* **All-vs-all alignment.** It costs ~1-2 s per container of 20-30 reads,
  and the pair count grows with the square of the read count; containers with
  more than 60 reads hit the 20 s timeout and phase nothing. Every read pair
  was aligned twice (q on t and t on q, `-D` without `-X`).
* **lamassemble** is ~40% mafft, ~27% lastal (`-m 50` from d35b23d) and the
  rest its own Python waiting on their pipes. Its cost is inherent to the
  method and nearly the same in both modes.
* **KMeans** on a few dozen points started OpenMP threads on all 24 cores:
  361 s CPU for 17 s wall per batch, which oversubscribes the machine when
  batches run in parallel.

## Changes

1. **Each read pair is aligned once** (`PhasingParams.align_pairs_once`,
   default on): `minimap2 --dual=no`, and `parse_pair_both` derives the other
   direction by inverting the alignment (insertions <-> deletions, mismatch
   bases from t, coordinates mirrored for reverse pairs). Exact: across 5.9 M
   `=` columns no base disagrees under the mapping (`check_inversion_exact.py`).
2. **Reads are put on the reference strand before phasing**
   (`consensus.orient_reads_to_reference`). Without this, align-once lost
   accuracy (pair acc. 0.974 -> 0.964). For a reverse-strand pair the
   inverted alignment places ambiguous gaps (homopolymers) at the other end,
   so all opposite-strand reads share "differences" on a target and get split
   by strand. With one strand, the inverted alignments equal minimap2's own
   other direction (mismatch Jaccard 0.99-1.00, base agreement 0.999;
   `check_inversion.py`).
3. **At most 50 reads are phased** (`PhasingParams.max_reads`, the longest;
   the others are re-added by alignment as unassigned reads always were).
   Bounds the quadratic tail; affects 13 of 2004 containers at 20x, but far
   more at higher coverage.
4. **Phasing Python on numpy**: CIGAR parsed from the string, alignments held
   as arrays, the SNV scan as a read x candidate-column matrix, recurrence as
   matrix products, net indels by prefix sums. Same algorithm: identical
   partitions on all 2004 containers except the timeout-bound ones, in both
   alignment modes. 0.39 -> 0.145 s CPU per container.
5. **KMeans / spectral clustering limited to one thread** (`threadpoolctl`).
   Legacy batch 25: Python CPU 395 s -> 8 s.

## Phasing accuracy and time (`phasing_eval.py`, all 2004 containers)

Trio read labels as truth (`../ava_phasing/data/trio_read_labels.tsv.gz`);
pair accuracy over the trio-labelled reads that the phasing assigns.

| variant | phasing time (sum) | pair acc. | containers perfect | trio reads assigned |
|---|---|---|---|---|
| before (d00f5ec) | 5044 s | 0.9737 | 0.944 | 0.827 |
| oriented, both directions | 4529 s | 0.9736 | 0.942 | 0.828 |
| **oriented, align once (new default)** | **2657 s** | **0.9734** | **0.939** | **0.826** |
| + `-k15 -w10` | 2239 s | 0.9714 | 0.938 | 0.827 |
| + `-k19 -w10` | 1947 s | 0.9708 | 0.937 | 0.827 |

Timed at 12 processes in parallel (the first row with other jobs on the
machine, so it reads a little high). Align-once
changes partitions about as much as orienting the reads does (ARI 0.981 vs
0.983); per container 37 worse / 26 better (orientation alone: 23 / 23).
Sparser minimizers (`minimap_params`) save another 20-30% but lose
systematically (16 worse / 8 better vs align-once): left as an option,
not the default.

## End to end

E2E_PLACEHOLDER

## Not done / next

* lamassemble is now the largest cost in both modes. Cutting it means a
  different consensus method (see `feature/poa-consensus`) or fewer/smaller
  assemblies. Small wins: lamassemble checks `mafft --version` on every call
  (~7% of its time).
* `summed_indel_distribution` (needed only by the KMeans fallback) and
  `get_read_alignment_intervals_in_cr` are ~2% of a batch.
* Some containers take the full 20 s timeout with only ~27 reads (e.g.
  crID 1480): repeats, where minimap2 reports many secondary alignments per
  pair.

## Files

`profile_batch.py` (+ `run_profile.sh`) stage/tool profile of real batches,
`dump_phasing_reads.py` the phasing read sets per container,
`phasing_eval.py` time + trio accuracy of `phase_reads` variants,
`profile_phasing.py` CPU profile of the phasing Python with replayed
alignments, `check_inversion*.py` inversion checks, `svp_variants.yaml` the
svp_improvements variants `perf_legacy` / `perf_phased`, `py.sh` runner.
