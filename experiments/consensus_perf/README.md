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

## Whole-genome extrapolation (crude)

`wg_estimate.py`: esc2_legacy and esc2_phased (same code, c1c3f4c), HG002
20x. The 22 blocks cover 151.7 Mb, 5% of hg38 and ~5% of the truth SVs, so
per-region costs are scaled x20.4. Run CPU is from `/usr/bin/time` (user +
sys, all tools).

| | legacy | phased | phased / legacy |
|---|---|---|---|
| blocks: `svirlpool run` CPU | 3397 s | 4293 s | 1.26x |
| blocks: consensus batches | 2640 s | 3598 s | 1.36x |
| blocks: everything else (wall) | ~134 s | ~145 s | |
| **whole genome: CPU** | **~19 h** | **~24 h** | 1.26x |
| whole genome: consensus CPU | ~15 h | ~20 h | |
| whole genome: consensus batches | ~610 | ~610 | |
| **whole genome: wall, 16 cores** | **~1.5-2 h** | **~2-2.5 h** | |
| whole genome: wall, 64 cores | ~1 h | ~1 h | |
| peak memory | ~13 GB (one step) + ~1-1.5 GB per running batch | same | |

* Wall = consensus CPU / cores, plus the non-consensus stages scaled as they
  ran (~45-50 min if their parallelism does not grow with the genome; less
  if it does). At 64 cores those stages dominate, so both modes take about
  the same time.
* The 12.5 GB peak is `consensus_align_to_initial_reference` (one job),
  most likely minimap2's GRCh38 index, which would not scale; the batches
  need <= 1.5 GB each.
* Legacy here includes its escalations (the spectral all-vs-all timeouts). Without
  them (perf_legacy: degraded, 2067 s) its consensus would be ~12 h, and
  phased would cost 1.7x.
* Likely low: the blocks are autosomal and hold ~5% of the SVs. The rest of
  the genome also has centromere-adjacent sequence, segdups and chrX/Y,
  where repeat containers (the expensive ones) are denser. Containers with
  CN > 4 are skipped in both modes. Cost grows with coverage: the all-vs-all
  grows with the square of the reads up to the 50-read cap, lamassemble
  roughly linearly. Read at 20x as +-50%.

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

## Ablation: phasing sites from all-vs-all vs reference alignments

Branch `exp/ref-vs-ava-phasing` (baba158, on a4a8558). Does aligning the reads
to each other beat reading the same variants off their reference alignments?
`ref_read_phasing.phase_reads_reference` (arm B, `--phasing-sites reference`)
builds the same per-target-read SNV and SV sites from the BAM alignments and
the reference FASTA; everything after the sites (`read_phasing.phase_from_sites`:
weights, correlation clustering, refinement), the read set (`select_reads`, 50
longest), the window (CR +- 10 kb) and the low-quality rule (read divergence to
the reference put on the pairwise scale) are shared with arm A (`phase_reads`).

`ref_vs_ava_eval.py` (lam_orient workdir, 8 processes, 217 s wall) +
`ref_vs_ava_summary.py`; per-container table `results/ref_vs_ava_all.tsv.gz`.
All 2004 HG002 containers:

| | A: all-vs-all | B: reference |
|---|---|---|
| phasing time (sum) | 1547 s | **48 s** |
| status phased / single / no info | 1792 / 124 / 88 | 1787 / 91 / 126 |
| pair accuracy, assigned trio reads (1768 paired containers) | **0.9838** | 0.9707 |
| containers perfect | **0.958** | 0.926 |
| trio reads assigned | 0.844 | 0.855 |
| containers better than the other arm | **93** | 27 (sign test p = 1e-9) |
| pair errors: split / joined | 641 / 901 | 1632 / 2382 |
| T2TQ100 allele set recovered (1747) | 0.974 | 0.969 |
| ... recovered by this arm only | 13 | 5 (p = 0.10) |
| purity | 0.9947 | 0.9923 |

By trio-label coverage (labelled / all reads of the container; the labels
come from reference alignments, so low coverage marks reads the reference
handles badly), pair accuracy A / B: < 0.6: 0.868 / 0.773 (148 containers);
0.6-0.8: 0.978 / 0.942 (145); >= 0.8: 0.993 / 0.987 (672). TR 0.987 / 0.972,
non-TR 0.977 / 0.968. At the allele level, the difference is in compound
hets in TRF (A only 9, B only 4) and hets in TRF (3 / 0).

* Most containers are partitioned identically (median difference 0); the gap
  is a tail of containers where B collapses, almost all with low label
  coverage. Two kinds: B sees far more SNV sites than A (crIDs 258, 2002,
  22: 844 vs 170, 980 vs 83, 2672 vs 1329; mismapped / paralogous reads
  differ from the reference, not from each other), or far fewer (995, 508,
  1614: 0 vs 13, 0 vs 4, 4 vs 59; het SNVs the reference alignments miss).
* The read-level truth is itself reference-based (trio SNVs in the reads'
  reference alignments), so it is evaluated only where reference alignments
  work, which favours B; A wins anyway. The allele-level difference is not
  significant; the end-to-end trio benchmark decides
  (`svp_variants_phasing_sites.yaml`: ps_ava / ps_ref).
* Both arms merge a het SV without a het SNV within the window: one SV is one
  discriminating position, `min_discriminating` is 2
  (`tests/test_ref_read_phasing.py::test_a_deletion_alone_separates_the_haplotypes`).

### Larger set (30%) and a tiered phasing

5% was too small to fit routing rules safely, so a second, disjoint set was
built: `make_regions30.py` (seed 30) draws 5 Mb tiles over 31.6% of the
autosomes' non-N sequence (871 Mb, 128 blocks, `results/regions30.bed`),
avoiding the 5% blocks +- 1 Mb. Upstream stages only (snakemake targets
`crs_containers.db consensus_batches.tsv copy_number_tracks.bed.gz`, 12
threads, in `~/development/phasing_ablation_30pct/`): 8957 containers. Trio
read labels with `../ava_phasing/trio_read_labels.py` (208,785 informative
SNVs; 61.7% of reads labelled vs 77.3% on the 5% set, minority votes 0.85%
vs 0.77% -- random tiles include more repeat-rich and duplicated sequence;
two containers spanning two chromosomes have no labels), truth columns with
`../ava_phasing/container_truth_light.py` (same categories as
container_truth.tsv on the 5% set for all 2004 containers).

A vs B on the 30% set (`results/ref_vs_ava_30pct.tsv.gz`):

| | A: all-vs-all | B: reference |
|---|---|---|
| phasing time (sum) | 10368 s | 1062 s |
| pair accuracy (7142 paired containers) | **0.9688** | 0.9572 |
| containers perfect | **0.921** | 0.897 |
| containers better than the other arm | **412** | 242 |
| allele set recovered (6615) | **0.965** | 0.961 |
| ... by this arm only | **50** | 25 (p = 0.005) |
| containers with 3-4 alleles | 391 | 726 |

The 30% set confirms the 5% result and makes the allele-level difference
significant (mostly compound hets in TRF, 29 vs 11). B over-splits more.

**Routing** (`ref_route_features.py`, `ref_route_rules.py`): when can B stand
in for A? B's own partition is the best guide -- its *discordance* (share of
the grouped reads' bases at the kept SNV columns that differ from their
group's majority; AUC 0.81 for B failing on the 5% set), then allele
imbalance at those columns (0.78), read divergence (0.70); MAPQ is weak
(0.65; most failures have MAPQ 60). Rule grid over allele count, discordance,
imbalance, smallest group share and divergence; a rule must not lose alleles
on balance and not lose read-level containers beyond a budget.

* Rules fitted on the 5% set and tested on the disjoint 30% set: the simple
  rule **B finds 2 alleles, discordance <= 0.02, smallest group >= 30%**
  (also the rule leave-one-chromosome-out picks in 13 of 22 folds on the 5%)
  holds: 30% set routed 42.9%, 37.4% of A's time avoided, pair accuracy
  0.9694 vs 0.9688 always-A (95% bootstrap CI of the difference
  [-0.00002, +0.0012]), containers worse / better 15 / 32, alleles lost /
  gained 2 / 5. Per chromosome the accuracy change is within +-0.003, no
  chromosome loses more than 1 allele set; routed share 13% (chr21) to 56%.
* Looser rules that the 5% set favoured do not transfer: at a 0.2% budget
  the 5%-fitted rule loses 26 alleles and gains 6 on the 30% set. Leave-one-
  chromosome-out on the 30% set (and on both, 10961 containers) picks
  "2 alleles, smallest group >= 30%" without the discordance limit: ~70%
  routed, ~60% of A's time avoided, accuracy -0.001 to -0.002, alleles
  neutral (8 / 15, 12 / 16) -- more churn; not adopted.
* Net saving: B runs for every container in tiered mode. Its mean cost on
  the 30% set was 0.12 s (median 21 ms; a tail of high-read containers took
  up to 10 s), so the simple rule saved ~27% of the phasing time net on the
  30% set and ~47% on the 5% set. The tail is gone since b411f8c (below):
  >= ~31% net on the 30% set.

`--phasing-sites tiered` implements the simple rule
(`ref_read_phasing.accept_reference`, `TIERED_MAX_DISCORDANCE` /
`TIERED_MIN_GROUP_FRACTION`): the reference-site phasing first, kept when it
passes, the all-vs-all phasing otherwise. Its decisions match the rule of the
feature table on 2002 of 2004 containers of the 5% set (the two differ by the
de-duplication of alignments listed under several CRs). End to end:
`svp_variants_phasing_sites.yaml` (ps_ava / ps_ref / ps_tiered).

### B's slow tail (b411f8c)

Profile of the 10 slowest containers (57 s): 81% in `recurrent_sites`, which
multiplied int32 site x read matrices -- NumPy has no BLAS path for integer
products. Then `discriminating_positions` (11%, set intersections per site
for every cluster pair) and `pair_weights` (4%). The reference arm made it
worse by giving every target read its own copy of a column's read lists.
Changes (shared with arm A):

* float32 products (exact for counts < 2^24), the ratio in float64 as before;
* the site x read memberships built once, a list shared by several sites
  indexed once; the reference arm shares a column's lists between targets;
* `discriminating_positions` as matrix products over those memberships;
* `pair_weights` sums with `bincount`, which adds in the loop's order.

`ref_b_bench.py` freezes B's inputs (150 slowest + 350 random containers of
the 30% set) and compares old and new code: **every result identical**
(status, groups, low-quality reads, discordance, digests of the SNV / SV
sites), total 480 -> 121 s (4.0x), slowest 9.6 -> 1.8 s. Arm A on 120 of
them: identical (its time is the minimap2 all-vs-all).

### End to end on a representative 15% (svp_tiered15)

`make_regions15.py`: 5 Mb tiles, 15% of every autosome (410 Mb, 72 blocks,
`results/regions15.bed`), disjoint from the 5% blocks +- 1 Mb where the
tiered rule was fitted; of 5000 random draws the one closest to the whole
autosomal genome (largest relative deviation 0.9%):

| | genome | 15% set |
|---|---|---|
| TRF share | 0.0624 | 0.0623 |
| GC | 0.4095 | 0.4060 |
| T2TQ100 SVs / Mb | 15.37 | 15.42 |
| V5 SVs / Mb | 13.54 | 13.67 |
| benchmark-region share | 0.942 | 0.937 |

`~/development/svp_tiered15` is a copy of svp_improvements (own results/,
shared resources and truvari env) with these regions, the trio and truvari
`--pctsize 0.9`; variants ps_ava, ps_tiered, ps_ref at b411f8c.
`consensus_q100.py` scores every HG002 core consensus against the Q100 v1.1
assembly (located with minimap2 asm20 via the padded consensus, the core
alone aligned with edlib in infix mode to each haplotype);
`tiered15_report.py` collects SV benchmarks, Mendelian consistency, run time
and the consensus identities, paired by container.

Run (2026-10-02): the three variants of a sample together, 8 threads each
(`run_grouped.sh` / `run_rest.sh` in svp_tiered15), so every variant sees
the same load -- the consensus timeouts are wall-clock, and the machine ran
at load 50-100 throughout. One run at a time on 24 threads (tried first)
left 22 cores idle for an hour per run: the 15% set includes the chr19
pericentromere (modelled alpha-satellite sequence passes the non-N filter),
whose containers escalate to the last level and time out by the dozen (one
batch of 100 containers ran > 2 h). The last escalation level is capped at
8 threads instead of 12 by this setup, for all variants alike. Report:
`results/tiered15_report.txt`.

SV benchmark (HG002, truvari --pctsize 0.9, refined F1) and trio Mendelian
consistency:

| | V5 all | V5 non-TRF | T2TQ100 all | T2TQ100 non-TRF | MC |
|---|---|---|---|---|---|
| ps_ava | 0.8487 | 0.9342 | 0.8520 | 0.9295 | 93.36% |
| **ps_tiered** | **0.8505** | **0.9359** | **0.8541** | **0.9311** | **93.52%** |
| ps_ref | 0.8462 | 0.9319 | 0.8467 | 0.9267 | 93.06% |

(T2TQ100 all: tiered TP 2704 / FP 285 / FN 613 vs ava 2697 / 293 / 617; ref
loses recall, 2639 TP.)

HG002 core consensuses vs Q100 (per container `repr` = mean over the two
haplotypes of the best consensus' identity, recall-like; `prec` = mean
identity of the container's consensuses, precision-like):

* ps_tiered routes 1961 of 4528 phased containers (43%) to the reference
  sites, as offline (42.9%). Paired vs ps_ava: repr +0.00013 (10 better / 7
  worse by > 0.005), prec +0.00003. Without the 57 containers with a
  last-level timeout in either run, every container phased on the all-vs-all
  sites is identical in both runs; all differences are in the routed ones
  (repr 7 better / 2 worse).
* ps_ref: per consensus worse (identity >= 0.99: 85.4% vs 90.9%; 2.01
  consensuses per container vs 1.87), per container apparently better (repr
  +0.005) -- but all of that is the chr19 pericentromere: there the
  all-vs-all phasing fails in the satellite arrays (no_information -> one
  consensus at ~0.51 identity to Q100) and the reference sites split the
  reads into 4 groups at ~0.97. Outside chr19 24-30 Mb ps_ref is worse (repr
  -0.00076, 58 better / 81 worse; prec -0.00019), ps_tiered is +0.00008
  (7 / 3). `tiered15_q100_diff.py`.

Time (sum over the trio; wall clock per run is not comparable, the variants
shared the machine):

| | all-vs-all phasing | consensus batches | CPU of `svirlpool run` |
|---|---|---|---|
| ps_ava | 16.74 h (14253 phasings) | 58.6 h | 62.9 h |
| ps_tiered | 12.78 h (8046) | 55.0 h (-6%) | 63.0 h |
| ps_ref | -- | 43.2 h (-26%) | 48.4 h |

The tiered mode avoids 24% of the all-vs-all phasing time (`phase_time.py`,
from the log timestamps, measured under load), but the all-vs-all phasing is
only ~28% of the consensus time on this set: -6% consensus time, and no
difference in total CPU within its noise (escalations, contention). Most of
the consensus time here is lamassemble, much of it timing out in satellites
and repeats.

Conclusions: (1) the tiered mode is safe -- SV F1, Mendelian consistency
and consensus accuracy equal to or slightly better than all-vs-all only;
(2) its saving is small end to end, because the phasing is a minority of the
consensus cost on a genome-representative set; (3) reference-site phasing
alone is worse (recall, MC, consensus precision) except in centromeric
satellites, where the all-vs-all phasing fails -- a candidate for a
satellite-specific route (tiered accepts only 2-allele results, so it does
not take these).
