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
| oriented, align once, numpy | 2657 s | 0.9734 | 0.939 | 0.826 |
| **+ max 50 reads (new default)** | **2692 s** | **0.9731** | **0.939** | **0.826** |
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

svp_improvements, `svp_variants.yaml` here; truvari refined F1; consensus =
sum of the 30 batch wall times (snakemake benchmarks, one job per batch).

| variant | V5 all | V5 non-TRF | T2TQ100 all | T2TQ100 non-TRF | consensus | longest batch |
|---|---|---|---|---|---|---|
| flank200_subset (legacy, d01a165 line) | 0.8552 | 0.9287 | 0.8585 | 0.9249 | 2046 s | 141 s |
| perf_legacy (this branch, legacy) | 0.8552 | 0.9287 | 0.8585 | 0.9249 | 2067 s | 132 s |
| ava_phase2 (phased, ac17af6) | 0.8619 | 0.9503 | 0.8685 | 0.9498 | 6744 s | 436 s |
| perf_phased (this branch, phased) | 0.8558 | 0.9279 | 0.8606 | 0.9274 | **4617 s** | **290 s** |

* Legacy: identical calls. The thread limit saves CPU, not wall time, for a
  batch running alone.
* Phased: the consensus stage is 32% faster (the phasing overhead over legacy
  4700 s -> 2570 s). The rest is lamassemble (~1900 s) and the remaining
  all-vs-all.
* **The phased F1 drop is not a phasing loss.** 15 of the 19 extra V5
  non-TRF FPs come from two containers whose all-vs-all hit the 20 s timeout
  in ava_phase2 (so they fell back to one consensus) and now finish (`fp_containers.py`):
  container 1303 (chr4:7.86 Mb, 13 over-split DELs; phased perfectly, trio
  pair accuracy 1.0 on 23 reads, 1522 SNV sites) and 1741 (2). Over all
  regions, containers that switched from timeout to phased account for +17
  FPs net; the rest of the changes roughly balance (+14 / -9).
  So part of ava_phase2's gain came from the **wall-clock timeout** sending
  paralog-like containers to a single consensus. That depends on the
  machine's speed. It also corrects the repeat-gate analysis in
  `../ava_phasing/README.md`: container 1303's 13 FPs went away through the
  timeout, not the phasing.
* SNV-site density does not separate the containers that timed out (15-157
  sites per read; 163 containers have > 40, 6 of them timed out), so there
  is no deterministic stand-in for the timeout here. The FPs of 1303 come
  from calling on two correct haplotype consensuses.

## Timeout escalation

The workflow meant to run timed-out loci again with more cores (threads
`[1, 2, 4, 12]`, timeouts `[20, 60, 60, 120]` s per snakemake attempt), but
it never did. The nested snakemake had no `--retries`. A tool timeout was
caught inside the batch and never failed the job. `process_consensus_container`
was called with `threads=1` hard-coded. And a retry would have rerun all
100 containers of the batch.

Now the escalation happens per container, inside the batch
(`consensus --escalation`, `svirlpool run --consensus-escalation`, default
`1:20,4:60,12:120` = THREADS:SECONDS). A container starts at the first
level. If the phasing or spectral all-vs-all or lamassemble timed out
(`tool_timeouts.record`), it is processed again at the next level. The last
level is the hard ceiling: after it, the degraded result is kept. Threads
are capped at the snakemake cores and the process's CPU affinity. Memory
cannot be handed to a running process, so it stays at the job level: a
failed batch job (e.g. OOM-killed) is retried by snakemake with double
`mem_mb` (2 GB up to `--consensus-max-mem-mb`, 16 GB) and runtime.

The 7 containers that timed out in perf_phased: 2 finished at 1 thread this
time, 4 at 4 threads / 60 s (1275: 54 s), 1502 at 12 threads / 120 s; none
stayed unresolved. The phased mode's result therefore no longer depends on
the machine's speed, except for containers that exceed the last level.

End to end (`esc_legacy` / `esc_phased` in `svp_variants.yaml`, commit
2647461; `escalations.py`, `counts.py`):

| variant | V5 all | V5 non-TRF | T2TQ100 all | T2TQ100 non-TRF | consensus | longest batch |
|---|---|---|---|---|---|---|
| perf_legacy | 0.8552 | 0.9287 | 0.8585 | 0.9249 | 2067 s | 132 s |
| esc_legacy | 0.8544 | 0.9275 | 0.8584 | 0.9237 | 2826 s | 419 s |
| perf_phased | 0.8558 | 0.9279 | 0.8606 | 0.9274 | 4617 s | 290 s |
| esc_phased | 0.8536 | 0.9279 | 0.8577 | 0.9274 | 5322 s | 416 s |

* Phased: 10 containers escalated, 8 finished at 4 threads, 2 at 12
  (1275, 1502), none unresolved. Outside TRF the calls are identical. In TRF
  the escalated containers now phase and are called on their haplotype
  consensuses: V5 all FP 140 -> 149, TP 1426 -> 1424. The difference is almost
  all container 1275 (chr3:195.48 Mb, 3 alleles, the MUC4 VNTR) and 218
  (chr10:126.9 Mb): the same kind of over-split repeat calls as container
  1303.
* Legacy: 8 containers escalated (the spectral all-vs-all and
  lamassemble timed out), 1059 is still unresolved at 12 threads / 120 s. It
  costs more there: at every level the spectral all-vs-all retries itself up
  to 4 times on subsampled reads, each with the level's timeout (568 and
  1076: ~180 s, 1059: 206 s). V5 all FP +3.
* Time: the escalated containers take 789 s (phased) and 897 s (legacy) of
  the batch wall time. In perf_* they cost ~20 s each and fell back. Wall
  times are from two variants running at once on 24 cores, so +-10%.
* So the escalation buys machine independence, not accuracy: the loci it
  resolves are repeat containers where the resolved result calls slightly
  worse than the fallback did.

## Repeat seeds in the phasing all-vs-all

Why the escalated containers were slow (`ava_modes.py`, `ava_params.py`):
not the base-level alignment bandwidth (`-r500` is slower than `-r2k`), and
not redundant secondary alignments. Aligning every pair separately
(one minimap2 call per target read, `-N 0`) is slower still: 218 takes 58 s
against 37 s, because the joint call skips about half the pairs (117 of 231).
It is chaining. In a tandem repeat every copy of a minimizer in every read
matches, and the anchors explode. Masking them fixes it: with `-U a,b`,
minimap2 does not seed with minimizers occurring more than `b` times in the
pooled reads. The bounds must differ; `-U x,x` falls back to the default.

All 2004 containers, reads cut to CR +- 10 kb, 12 processes, 300 s timeout
(`phasing_eval.py --timeout 300`), trio read labels as truth:

| minimap2 seeds | phasing time (sum) | slowest | pair acc. | containers perfect |
|---|---|---|---|---|
| default (before) | 2959 s | 274 s | 0.9727 | 0.937 |
| **`-U15,20` (new default)** | **1703 s** | **6.4 s** | **0.9753** | **0.945** |
| `-U10,15` | 2262 s | 7.0 s | 0.9742 | 0.942 |
| threshold 1.0 x reads | 1690 s | 12.5 s | 0.9735 | 0.940 |
| threshold 0.75 x reads | 1618 s | 6.9 s | 0.9732 | 0.941 |
| threshold 0.5 x reads | 3455 s | 8.4 s | 0.9664 | 0.925 |

By read count: up to 40 reads `-U15,20` matches or beats the default. At
40-50 reads (23 containers; phasing is capped at 50) the default reaches
pair accuracy 0.74, `-U15,20` 0.90, in a seventh of the time. Thresholds that
drop below about 10 (0.5 x reads in small containers) leave too few seeds;
the result is worse and even slower. What matters is an absolute floor, not
scaling with the read count.

Also, below the last escalation level a container is now abandoned at its
first timeout (`tool_timeouts.Escalate`). It no longer finishes a degraded
result first: the phasing fallback, or legacy's all-vs-all retries on
subsampled reads.

End to end (`esc2_*`, c1c3f4c):

| variant | V5 all | V5 non-TRF | T2TQ100 all | T2TQ100 non-TRF | consensus | longest batch | escalated |
|---|---|---|---|---|---|---|---|
| ava_phase2 | 0.8619 | 0.9503 | 0.8685 | 0.9498 | 6744 s | 436 s | (17 timeouts) |
| perf_phased | 0.8558 | 0.9279 | 0.8606 | 0.9274 | 4617 s | 290 s | (7 timeouts) |
| esc_phased | 0.8536 | 0.9279 | 0.8577 | 0.9274 | 5322 s | 416 s | 10 |
| **esc2_phased** | **0.8577** | **0.9301** | 0.8582 | **0.9296** | **3598 s** | **214 s** | **2** |
| esc_legacy | 0.8544 | 0.9275 | 0.8584 | 0.9237 | 2826 s | 419 s | 8 |
| esc2_legacy | 0.8547 | 0.9275 | 0.8584 | 0.9237 | 2640 s | 413 s | 7 |

Phased mode now has no phasing timeouts at all. The two escalations are
lamassemble (965, 1931), both resolved at 4 threads. Phased mode costs 1.74x
legacy (ava_phase2: 3.3x), with every locus resolved at F1 at or above
perf_phased (V5 FP 140 -> 144, FN 341 -> 327). Legacy mode keeps its
spectral all-vs-all timeouts (`-U 25,35`, 1059 unresolved); the early abort
saves ~190 s.

With the spectral all-vs-all at `-U 15,20` too (`esc3_legacy`, f70cdbc; reverted: legacy mode is to be removed, so it stays as it is):

| variant | V5 all | V5 non-TRF | T2TQ100 all | T2TQ100 non-TRF | consensus | longest batch | escalated |
|---|---|---|---|---|---|---|---|
| perf_legacy | 0.8552 | 0.9287 | 0.8585 | 0.9249 | 2067 s | 132 s | (timeouts degraded) |
| esc3_legacy | 0.8524 | 0.9287 | 0.8631 | 0.9249 | 2481 s | 354 s | 6, all resolved |

Every locus is resolved, including 1059. Outside TRF the calls are
unchanged. In TRF, V5 loses 0.3 points and T2TQ100 gains 0.5. 568, 1059 and
1076 still need 12 threads (150-170 s each): the spectral path has costs
beyond the seeds that were not looked into.

## Not done / next

* With the escalation, 1303-like containers (paralog-rich, phased correctly,
  over-split calls) are phased everywhere. Decide what they should give:
  a question for calling, not for runtime.
* On SLURM the batch job reserves one CPU, and the affinity cap keeps the
  escalation to that CPU (only the timeout grows). Escalating cores there
  needs the job to request them, or a follow-up job for the escalated
  containers.
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
svp_improvements variants `perf_legacy` / `perf_phased`, `e2e_table.py` F1 and
consensus time per variant, `compare_consensus.py` / `fp_containers.py`
consensus and FP differences between two variants per container, `py.sh`
runner.
