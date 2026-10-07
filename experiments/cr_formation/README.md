# Candidate-region formation: fewer long and deep regions

Branch `perf/assembly-tail`. Assembly (lamassemble) takes most of the
consensus time, and most of that goes to a few long or deep candidate regions
(CRs) in repeats, many of which time out. In a real run, the only evidence
about a locus is its SV signals, so this asks whether CR formation can avoid
those regions from the signals alone. T2TQ100 is used only to score a CR set,
never to form one.

**State (2026-10-07): `cr_s1_sat` is the candidate worth keeping.** Its rules
are opt-in options, off by default. Next step (not run): the same rules with
the merge buffer left at 1200; see "Next".

## Why long and deep CRs form

* The signal strength measures density, not agreement between reads. It sums
  the scores of other reads' signals within 1 kb (a flat 10 above 200
  signals), so dense noise in segdups and satellites passes `filter_absolute`
  and `filter_normalized`, and `cr_merge_buffer` (1200) chains it together.
* Every TRF interval goes into the `bedtools merge` together with the
  signals. The TRF bed holds whole satellite arrays of 10-62 kb, so a few
  signals in one stretch the CR over the whole array.
* The read-count filter (more than 6x the median reads) drops only the largest
  of these (20 proto-CRs of 0.1-2 Mb in centromere models). Arrays just
  below that limit make up the expensive tail.
* Removing or capping the TRF bridging made things worse: heavy cost went from
  50 to 92-106 Mbp, because the arrays split into pieces that each pass the
  read-count filter. The neighbour merge (`has_similar_sv_signals`) almost
  never fires, so turning it off changes nothing.
* A depth cap alone loses true positives (8-63). Depth is better handled by
  the assembly read cap (`--assembly-max-reads`) than by dropping CRs.
* 82% of the CRs that timed out lie completely outside the T2TQ100 benchmark
  bed, so the absence of a truth SV there says nothing. Inside the bed, long
  and timed-out CRs mostly hold real SVs, mostly large deletions.

## The rules (`signalstrength_to_crs`, also options of `svirlpool run`)

| option | effect |
|---|---|
| `--min-signal-support N` | drop signals that fewer than N other reads confirm: same type (both deletion ends count as one type, as do both BND types), within max(100, size/2) bp or in the same tandem repeat, size ratio >= 0.5. Support is counted over all signals, before the strength filters. |
| `--satellite-depth-factor k` | drop CRs that lie >= 90% in tandem repeats of >= `--satellite-min-tr-length` (10 kb) and are deeper than k x the median CR depth (median over the proto-CRs that pass the read-count filter). They go to the `--dropped` db. |
| `--cr-merge-buffer` | existing option, 1200 by default |

Config keys in `main.smk`: `min_signal_support`, `satellite_depth_factor`.
Both default to 0 (off), and with them off the output is unchanged (same
4,601 CRs on HG002 15%). Support for 511k signals takes 14 s: inside a tandem
repeat it is computed per read rather than per signal pair, which gives
identical values to the quadratic version (82 s).

## Offline: CR sets scored without running consensus

Cost proxy per CR = median depth x (span + 400). "Heavy cost" is the sum
over CRs with a cost >= 100 kbp, which was the best predictor of lamassemble
timeouts.

**Tuning set, HG002 15%** (truvari base set of T2TQ100 'all', 2,773 TPs, all
covered at baseline):

| rule | heavy cost (Mbp) | TPs lost |
|---|---|---|
| current (buffer 1200) | 49.9 | 0 |
| buffer 600 | 40.3 | 2 |
| support >= 2, buffer 600 | 36.1 | 0 |
| support >= 2, buffer 600, satellite 1.5 | 29.5 | 0 |
| support >= 2, buffer 600, satellite 1.25 | 24.7 | 1 |

**Validation, outside the tuning set** (`cr_validate.py`, `cr_lost30.py`).
HG003 and HG004 have no truth set, so "lost" there counts their own trio
calls that no CR covers any more. On HG002 30% (`phasing_ablation_30pct`, disjoint
regions) it counts the T2TQ100 SVs >= 50 bp in the benchmark bed that lose
their CR.

| rule | HG002 15% | HG003 15% | HG004 15% | HG002 30% |
|---|---|---|---|---|
| current | 49.9 | 50.7 | 41.6 | 98.5 |
| buffer 600 | 40.3 / 5 | 40.4 / 17 | 32.1 / 17 | 87.0 / 1 |
| support 1, buffer 600 | 36.7 / 5 | 36.2 / 23 | 30.3 / 18 | 84.1 / 2 |
| support 1, buffer 600, satellite 1.5 | 30.1 / 6 | 30.1 / 35 | 26.2 / 28 | 65.9 / 12 |
| support 2, buffer 600 | 36.1 / 3 | 36.6 / 15 | 30.9 / 9 | 84.8 / 38 |

(heavy cost in Mbp / lost)

* Support >= 2 does not generalise. On the 30% set it loses 38 SVs, at a median
  read depth of 2 near chromosome ends, where a heterozygous SV cannot have
  two confirming reads. Support >= 1 keeps nearly all of the saving.
* Most of the call losses on HG003 and HG004 come from buffer 600 alone, and
  almost all of them are at Mendel-consistent sites.
* On the 30% set the satellite rule's losses are mostly SVs >= 1 kb in long
  VNTRs.

## End to end: trio over the 15% set (svp_tiered15)

`run_grouped_cr.sh` ran the three variants of one sample side by side
(8 threads each) on the same code, followed by truvari, Mendelian consistency
and `consensus_q100.py`. The report is from `tiered15_report.py --variants
cr_base,cr_s1,cr_s1_sat --ref cr_base`. Variants are in
`../consensus_perf/svp_variants_cr_formation.yaml`.

| | cr_base | cr_s1 | **cr_s1_sat** |
|---|---|---|---|
| flags | (defaults) | `--min-signal-support 1 --cr-merge-buffer 600` | cr_s1 + `--satellite-depth-factor 1.5` |
| CRs, HG002 | 4,601 | 4,787 | 4,780 |
| CPU, trio (h) | 58.6 | 48.7 (-17%) | **41.0 (-30%)** |
| wall, trio (min) | 434 | 416 (-4%) | 420 (-3%) |
| summed consensus batch time (h) | 44.5 | 43.9 | 41.6 |
| containers escalated / timed out at the last level | 724 / 238 | 760 / 200 | 778 / **155** |
| F1 V5 all / non_trf | 0.8459 / 0.9279 | 0.8495 / 0.9262 | 0.8492 / 0.9268 |
| F1 T2TQ100 all / non_trf | 0.8482 / 0.9209 | 0.8510 / 0.9198 | 0.8514 / 0.9203 |
| Mendelian consistency (inconsistent) | 0.9347 (338) | 0.9286 (370) | 0.9318 (351) |
| Q100: consensus identity mean / share >= 0.99 | 0.9897 / 0.901 | 0.9896 / 0.890 | 0.9895 / 0.890 |
| Q100: container representation vs cr_base | | -0.0023 (584 worse, 483 better) | -0.0022 (574 worse, 492 better) |

* `cr_s1_sat` saves 30% of the CPU and a third of the last-level timeouts.
  Most of the saving is futile timeout work, so wall time barely moves: each
  run's wall time is set by its few slowest containers.
* SV calls: F1 on all regions is +0.003 (T2TQ100: 16 fewer FNs, 6 fewer FPs).
  Outside tandem repeats, TPs are unchanged and there are 1-2 more FPs.
* Cost: Mendelian consistency drops 0.3 points (0.6 for `cr_s1`), and the
  per-container Q100 representation drops 0.2 points. The pairing by container
  is rough, because the variants form about 190 more CRs with other boundaries.
* `cr_s1_sat` is at least as good as `cr_s1` on every measure, so the
  satellite rule costs nothing extra here. Both variants share buffer 600,
  which is the suspected source of the quality cost (split CRs, alleles
  assembled in pieces). This has not been checked.

## Per-job time limits (`time_limits.py`)

What a wall-clock limit per consensus job would cost: a container (one job
of the batch loop) that needs longer than the limit is dropped with all its
CRs, and so is every call whose consensuses all come from dropped
containers. The times are from the trio runs above (8 threads per run, three
runs side by side), with escalation retries included. TP and FP are counted
on truvari bench for HG002 vs T2TQ100 'all', before refine: 2,549 TPs and
457 FPs. Calls are counted in the family VCF: 6,472 for cr_base and 6,602 for
cr_s1_sat. "Consensus h" is the summed container time with the limit
applied.

| variant | limit | containers lost | CRs lost | calls lost | HG002 TP lost | HG002 FP lost | consensus h |
|---|---|---|---|---|---|---|---|
| cr_base | none | 0 of 13,592 | 0 of 13,676 | 0 | 0 | 0 | 43.5 |
| cr_base | 30 s | 751 (5.5%) | 807 | 850 | 52 | 33 | 20.5 |
| cr_base | 60 s | 356 (2.6%) | 401 | 486 | 6 | 15 | 24.8 |
| cr_base | 120 s | 245 (1.8%) | 290 | 303 | 2 | 10 | 29.4 |
| cr_base | 180 s | 196 (1.4%) | 236 | 229 | 1 | 10 | 33.1 |
| cr_base | 300 s | 114 (0.8%) | 144 | 139 | 0 | 0 | 38.1 |
| cr_s1_sat | none | 0 of 14,149 | 0 of 14,242 | 0 | 0 | 0 | 40.4 |
| cr_s1_sat | 30 s | 794 (5.6%) | 858 | 878 | 38 | 17 | 22.6 |
| cr_s1_sat | 60 s | 318 (2.2%) | 369 | 482 | 4 | 11 | 26.8 |
| cr_s1_sat | 120 s | 199 (1.4%) | 242 | 269 | 1 | 10 | 30.6 |
| cr_s1_sat | 180 s | 152 (1.1%) | 184 | 179 | 0 | 0 | 33.5 |
| cr_s1_sat | 300 s | 67 (0.5%) | 89 | 46 | 0 | 0 | 36.8 |

* The containers lost hold slightly more than one CR each on average (1.1-1.3).
  Over the whole run a container holds 1.006 CRs, so job and locus are nearly
  the same thing here.
* The calls lost are mostly outside the T2TQ100 benchmark: 303 calls but 2
  TPs at 120 s. Most of the slow containers lie outside the benchmark
  bed, so truth cannot tell whether those calls were real.
* Is this an artefact of the region subset? Containers join CRs only through
  split reads, and a subset cuts links to partners outside its regions.
  `container_sizes.py` compares the whole-genome HG002 20x run of the paper
  rerun (hg38, `svirltiles/fix/giab/hg38/20x/HG002`). The whole genome gives
  31,288 containers for 31,430 CRs (1.00 per container). 0.34% of containers
  hold more than one CR, the largest holds 7, and 24 span more than one
  chromosome. Restricted to the containers that touch the 15% regions, the
  whole genome gives 1.01 CRs per container, with 1.6% of CRs in
  multi-CR containers (subset: 1.0%). So containers are not larger
  genome-wide. The consensus batch jobs are, though: batches take up to 100
  containers within 20 Mb (`crs_to_batches`). In the subset, the gaps between
  regions cut them to a median of 52 containers (76 batches); in the whole
  genome they are full (median 100, 357 batches). A time limit belongs on the
  container, not on the batch job.
* A limit is a hard drop. The escalation timeouts inside a job (20/60/120 s per
  tool call) instead keep a degraded result. Moving those per-call limits
  changes quality, not the call count, and needs a rerun to measure.
* The times depend on the load: three runs shared 24 cores.

## Next (not run)

* The same rules with buffer 1200 (`--min-signal-support 1
  --satellite-depth-factor 1.5`). Offline (`cr_b1200.py`) this keeps about half
  of the heavy-cost cut: HG002 49.9 -> 37.2, HG003 50.7 -> 38.6, HG004
  41.6 -> 33.7 Mbp (buffer 600: 30.1 / 30.1 / 26.2). End to end it should show
  whether the Mendel and Q100 cost comes from buffer 600.
* The extra inconsistent sites of `cr_s1` / `cr_s1_sat` versus `cr_base`:
  are they at split CRs?
* The rules were tuned and validated on 20x ONT over the 15% and 30% sets.
  Whole-genome runs and other coverages are untested.
* They combine with the assembly-tail options (`--heavy-container-bp`,
  `--assembly-max-reads`), but that combination has not been run.

## Files

`crtools.py` loaders and the scorer (cost proxy, truth coverage);
`cr_reform.py` re-forms CRs with `signalstrength_to_crs` and caches them in
`$CR_STUDY_OUT` (default `./cr_study_out`); `cr_validate.py` the validation
tables; `cr_lost30.py` which rule loses the 30% SVs, and their depth;
`cr_b1200.py` the buffer-1200 estimate; `time_limits.py` the per-job
time-limit table; `q100_cr.sh` the Q100 comparison of
the three variants. The svirlpool dev env has no `edlib`, so this script uses
the svp_tiered15 truvari env's python. Run the Python scripts with the dev
env, from this directory. The paths to the svp_tiered15 and
phasing_ablation_30pct runs are hard-coded. The original tuning scripts are
not kept; the rules they found are now in `signalstrength_to_crs`.
