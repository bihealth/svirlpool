# Slow consensus containers: time limit, depth gate, read selection

Branch perf/assembly-tail. The data are the HG002 run of the 15% trio set (svp_tiered15, variant
cr_base) unless noted.

## 1. Time is spent on loci that resolve badly (`container_value.py`)

| container wall time | containers | share of consensus time | Q100 representation ≥ 0.99 | TP / FP calls |
|---|---|---|---|---|
| < 20 s | 4,253 | 29% | 92% | 2,459 / 414 |
| > 120 s | 84 (1.8%) | 53% | 22% | 2 / 11 |

Representation: for each Q100 haplotype, the best identity of any of the container's
consensuses, averaged over the two haplotypes.

## 2. Container time limit (`--container-time-limit`, 8d411c3)

This is a wall-clock budget per container over all escalation levels. It is implemented with
SIGALRM; the read-cache fetch is shielded, and the minimap2 helpers kill their process group.
A container still unfinished at the limit is dropped (empty result), and the batch carries on.

`svp_variants_time_limit.yaml` and `run_tl10.sh` set up the HG002 test on a disjoint 10% set
(svp_tl10). That run was stopped on request: the time limit will be tested together with read
selection (section 4).

## 3. A gate on signal depth is weaker (`coverage_gate.py`)

Depth here is the median `ExtendedSVsignal.coverage` of a CR, the maximum over the container's
CRs, relative to the median container.

| rule | containers excluded | share of consensus time | of the 84 slow containers | TP lost | FP lost |
|---|---|---|---|---|---|
| depth > 3× | 41 | 17% | 22 | 5 | 16 |
| time limit 120 s | 84 | 53% (5.2 h saved) | 84 | 2 | 11 |

Most slow loci have normal depth. Signal coverage also does not see piles of supplementary
alignments: chr22:18.89 Mb has signal depth 22, but its fetch returns 21,736 alignments.

## 4. Reads that resolve an allele (`locus_reads.py`, `locus_truth.py`, `readsel_consensus.py`)

Per CR, each read is classified by how it crosses the CR:

- **single:** one alignment covers the CR;
- **split:** collinear split alignments, from the SA tag, cover it;
- **anchor:** how far the read runs on beyond the CR, on its shorter side.

The truth: each read's segment at the CR is mapped (full read, map-ont) to Q100 v1.1. It is
judged against haplotype windows bounded by the nearest unique 5 kbp anchors (≤ 300 kbp away),
or, without a window, at chromosome level.

Locus kinds among the 84 slow containers:

| kind | examples | signature |
|---|---|---|
| anchored repeats (chr19 25–27.4 Mb) | 1888, 1907 | ~1 alignment per read, 30–90% of reads crossing, 2–5 kb anchors |
| unspanned arrays | chr12 34 Mb, chr3 90.6 Mb, chr9 64.9 Mb | CR span 30–68 kb > read alignments; no read crosses (9.6% of slow CRs) |
| fragment pile-ups | 2399 (10% set), 4100 | 4–5 alignments per read, ~80% supplementary, < 1% of reads cross |
| collapsed copies | 1758/1759 | ~1,050 reads (17× depth) on 1 kb, nearly all crossing, short anchors |

Truth, with k = 2 × the controls' median read count = 42:

- **Controls:** 98% of the crossing reads resolve at 0.988 identity; both haplotypes are
  covered in 99% of CRs.
- **Fragment pile-up 2399:** none of the fragment reads lie at the locus (100% map elsewhere).
  12 of its 17 crossing reads resolve, covering both haplotypes.
- **Collapsed copies and slow CRs with a window (20 of 94):** the top-42 crossers resolve in
  50% of cases, the remaining crossers in 1% (99% map to other copies). Per CR, the top-k
  resolve 87% vs 58% of the other reads. Both haplotypes are covered by the top-k in 80% of
  CRs, against 95% with all reads.
- **Truth coverage:** 69 of 94 slow CRs (the chr19 pericentromere) have no unique anchors within
  300 kbp. At chromosome level, the top-k cover both haplotypes in 67% of slow CRs against 95%
  with all reads.

Conclusions:

- Crossing and anchoring separate true-locus reads from paralog reads in pile-ups and collapsed
  copies.
- In spanned satellites the crossing reads of one haplotype are often missing, so a top-k cut
  can lose a haplotype.
- CRs without crossing reads cannot be resolved by spanning reads.
- Non-crossing reads do carry alleles that are longer than the reads (large insertions); a
  selection rule should only bite when a CR has more than k candidate reads.

`readsel_consensus.py` (not yet run) runs the pipeline's consensus on only the top-k reads.

## 5. Trio test: read selection + 240 s container limit (svp_rs10, 2026-10-07)

`--read-selection-factor 2 --container-time-limit 240` (rs_sel240) against the defaults (rs_base).
The data are the trio on the 10% region set (disjoint from the 15% tuning set), 8 threads per run,
three runs at a time. `report.py` is the cluster harness's (evaluation/tuning/workflow/scripts).

| | rs_base | rs_sel240 |
|---|---|---|
| CPU, trio (h) | 62.0 | 34.2 (−45%) |
| consensus batch time (h) | 45.7 | 31.4 (−31%) |
| wall, sum of runs (min; depends on load) | 456 | 379 |
| containers dropped at 240 s / CRs read-selected / CRs dropped (no crossing read) | – | 87 / 1,638 / 31 |
| V5 all P / R / F1 | 0.8752 / 0.7785 / 0.8240 | 0.8788 / 0.7788 / 0.8258 |
| V5 non_trf F1 | 0.9243 | 0.9263 |
| T2TQ100 all P / R / F1 | 0.8782 / 0.7841 / 0.8285 | 0.8827 / 0.7834 / 0.8301 |
| T2TQ100 non_trf F1 | 0.9176 | 0.9196 |
| Mendelian consistency | 0.9321 (5,622 sites) | 0.9336 (5,508 sites) |
| Q100 identity: mean / aggregate / share ≥ 0.99 | 0.9861 / 0.9682 / 0.869 | 0.9879 / 0.9756 / 0.880 |
| Q100 representation, paired by container | – | +0.00013 (68 better, 49 worse; 40 containers without consensus) |

The effect of each option alone is not yet separated.

## 6. Why the time limit works, and what could replace it (`container_cost.py`, `profile_containers.py`, 2026-10-08)

Data:

- `container_cost.py`: per container, the stage times parsed from the consensus batch logs of the
  cluster tuning runs (evaluation/tuning; regions10 = verify10 at 64 threads per run, regions15 =
  rs_tl_grid at 32 threads), joined with features of the container (crs_containers.db), of the
  reference (TRF, k-mer self-repetitiveness) and, for HG002, its value (Q100 representation,
  truvari T2TQ100 TP / FP).
- `profile_containers.py`: HG002 containers rerun locally (1 thread, no timeouts, temporary
  files on local disk) with timers around every step of lamassemble and of the phasing.

### 6.1 On the cluster, a container's wall clock is mostly waiting

Same settings (v_rs3_tl120), one node per sample:

| | HG002 | HG003 | HG004 |
|---|---|---|---|
| consensus batches ran | 10:04–11:16, with ~20 other runs | 11:08–11:48, most others done | 10:03–11:17 |
| container wall, summed | 54.9 h | 6.6 h | 53.7 h |
| median container | 41 s | 2.9 s | 40 s |
| dropped at 120 s | 88 | 19 | 77 |
| CPU of the whole run | 7.3 h | 10.1 h | 8.1 h |

- minimap2's own real time (its stderr) over the ~11,900 read-to-consensus alignments of a
  sample is 0.07 h. The log stage around those calls takes 14 h.
- The HG002 job used 7.3 CPU h in 85 min on 64 CPUs, i.e. its CPUs were 8% busy. So it is not
  oversubscribed; its processes sleep, waiting for I/O.
- Locally, typical containers take 0.4–1.5 s; the same containers took 32–272 s on the cluster
  (table in 6.3).
- Every consensus process of every concurrent run creates and removes its temporary
  directories (lamassemble's LAST database and MAFFT files, minimap2 inputs, ~10 per container)
  in one cephfs directory, `$TMPDIR` = /data/cephfs-1/scratch/.../tmp, and starts its tools
  from the pixi environment on cephfs. Metadata operations on cephfs are the likely bottleneck;
  it is not yet tested directly.

So the effective threshold of a wall-clock limit is set by the cluster load: the same 120 s
dropped 19 containers on a quiet node and 77–88 on busy ones. Because the lost containers differ
by sample, they cost Mendelian consistency too. Of HG002's 88 drops, 54 loci were lost in no
other trio member; with a `cut_bp > 300 kb` gate it is 11 of 49 (23 lost in all three).

### 6.2 What the limit should drop: assemblies that are expensive in themselves

In the quiet run, time is predictable from the reads (AUC for containers over 30 s / 100 s):

| feature | known | AUC > 30 s | AUC > 100 s |
|---|---|---|---|
| reads lamassemble gets × their bp (phased: the groups, ≤ 50 reads; unphased: all) | after phasing | 0.984 | 0.998 |
| cut-read bp | after the read fetch | 0.972 | 0.998 |
| reads | after the read fetch | 0.972 | 0.974 |
| CR span | before (containers DB) | 0.862 | 0.990 |
| k-mer self-repetitiveness of the reference at the CRs | before | 0.817 | 0.957 |
| SV signals | before | 0.661 | 0.949 |
| signal depth | before | 0.876 | 0.893 |
| TRF share of the CRs | before | 0.322 | 0.433 |

In the contended runs, every feature drops to AUC ~0.55 at 30 s and ~0.75 at 100 s: load noise.

- All 19 containers dropped in the quiet run have > 270 kb of cut reads. They lie in satellite
  arrays (chr5 46–49.6 Mb, chr22 12–16 Mb, chr18 15.4 Mb), and 105–116 of their 120 s are
  lamassemble.
- When the phasing finds no two alleles (`single` / `no_information`), the fallback assembles
  all of the container's reads in one lamassemble. The phased path assembles at most the 50
  phased reads. At equal size, unphased containers take 2–4.5× longer: at 100–200 kb of cut
  reads, a median of 41 s vs 9 s, and all 11 unphased ones over 400 kb hit the limit.
- These containers are worth little. Of HG002's 88 drops, the 30 with > 200 kb of cut reads held
  3 TP and 0 FP calls in the uncapped run. The other 58 were ordinary containers dropped only
  because of the load, with 47 TP, 20 FP and a Q100 representation of 0.984.
- They are expensive in CPU because they climb the escalation ladder: they time out at 1 thread,
  then at 4, and run at 12 threads. For HG002 (regions15), the run used 15.5 CPU h without a
  limit, 10.9 h at 240 s (8 dropped) and 7.6 h at 120 s (27 dropped). That is ~17 CPU minutes
  per dropped container. Levels that ended in a timeout make up 40–64% of the wall time of the
  escalated containers.
- lamassemble's own timeout is checked between its steps. It fires within 0–3 s on a quiet node
  but up to 56 s late under load, so that effect is minor.

### 6.3 Where a heavy container spends its compute (local profile)

HG002 containers of v_base (no read selection, no limit), rerun at 1 thread with no timeouts on
the workstation, temporary files on local disk:

| container | locus | reads / cut kb | cluster s (busy, escalations) | local s | of which lastal |
|---|---|---|---|---|---|
| 2561 | chr3 54.5 Mb | 24 / 23 | 32 | 1.5 | 0.1 |
| 3464 | chr5 179.9 Mb | 21 / 19 | 137 | 1.1 | 0.1 |
| 1068 | chr12 133.3 Mb (short TR) | 17 / 16 | 272 | 1.1 | 0.4 |
| 2389 | chr22 18.4 Mb | 1,305 / 3,647 | 125 | 16 | 0.3 |
| 2081 | chr21 10.0 Mb | 83 / 126 | 177 | 15 | 9.4 |
| 597 | chr4 190.1 / chr10 133.7 Mb | 69 / 550 | 260 | 48 | 38 |
| 3061 | chr5 48.0 Mb (unphased) | 80 / 287 | 142 | 68 | 59 |
| 1527 | chr18 15.4 Mb | 54 / 363 | 259 | 113 | 108 |
| 2211 | chr22 11.4 Mb | 91 / 886 | 373 | 343 | 335 |
| 3302 | chr5 49.4 Mb (unphased) | 149 / 792 | 400 | 615 | 597 |
| 3074 | chr5 48.1 Mb (unphased) | 119 / 884 | 469 | 695 | 681 |
| 2249 | chr22 12.4 Mb | 124 / 1,569 | 539 | 1,405 | 1,383 |

- **lastal does the work.** Over 28 lamassemble calls it is 98.5% of their time. The phasing's
  all-vs-all takes ≤ 7 s, and MAFFT, the MAF parsing, the layout and the anchors (Python) take
  ≤ 4 s each.
- **A call's time is quadratic in its input.** Time ∝ (assembled bp)^1.93 (log-log r = 0.97).
  At 1 thread that is ~9 s for 100 kb, ~77 s for 300 kb and ~790 s for 1 Mb. LAST's aligned
  columns predict it almost exactly (r = 0.99, exponent 0.94). The k-mer collision mass of the
  reads (Σ count², the seed hits of an all-vs-all) gives r = 0.95, but only after the short-TR
  calls are left out. In 1068 every k-mer recurs 13–17 times, yet LAST seeds nothing there.
- **2389 is the counter-example to "many reads = slow".** It has 1,305 reads, but the phasing
  hands lamassemble ≤ 50 of them, so it takes 16 s. What counts is what lamassemble gets.
- **LAST's `-m 50` multiplies the work in satellite arrays.** On one cluster each of 2211 and
  2249 (same oriented reads, LAST run as lamassemble runs it):

  | input | -m | alignments | read pairs linked | seconds |
  |---|---|---|---|---|
  | 2211, 14 reads, 256 kb | 5 | 1,278 | 87 / 91 | 12 |
  | | 10 | 4,850 | 91 / 91 | 32 |
  | | 50 | 32,116 | 91 / 91 | 204 |
  | 2249, 26 reads, 559 kb | 5 | 2,392 | 300 / 325 | 14 |
  | | 10 | 3,824 | 316 / 325 | 19 |
  | | 50 | 281,618 | 325 / 325 | 853 |

- **The strand fallback triples it.** In 2249, the one-strand layout left reads unlinked, so
  lamassemble redid LAST on both strands (6 lastal calls for 2 clusters).
- **The cost is CPU, not I/O.** These containers run minutes even on a quiet machine.

### 6.4 A gate on the assembly size against the realized time limits (HG002)

Value of the containers each rule removes: their calls in the uncapped run (truvari T2TQ100
all, before refine) and their Q100 representation.

| base run | rule | containers | TP lost | FP removed |
|---|---|---|---|---|
| regions15 rs3 (32 threads) | 120 s limit (rs3_tl120) | 27 | 5 | 15 |
| | 240 s limit (rs3_tl240) | 8 | 0 | 0 |
| | cut-read bp > 200 kb | 78 | 14 | 22 |
| | > 300 kb | 46 | 10 | 18 |
| | > 400 kb | 31 | 7 | 17 |
| regions10 v_base (64 threads, no read selection) | 120 s limit (v_rs3_tl120, busy) | 88 | 50 | 20 |
| | cut-read bp > 300 kb | 70 | 21 | 0 |
| | > 400 kb | 57 | 11 | 0 |

A bp gate removes about the set a 120 s limit removes on a quiet node (all 19 quiet-run drops
had > 270 kb). It removes it on every node, and it removes it before any work is spent.

### 6.5 Options, from the environment to the algorithm

1. **Make the wall clock reflect the work.** Give each consensus process a temporary directory
   on node-local storage. main.smk does not pass the consensus CLI's `--tmp-dir-path`, so every
   temporary file follows `$TMPDIR` (cephfs scratch on the cluster). Exporting a job-local
   `TMPDIR` in the consensus rule would need no code change. Expected: container times like the
   quiet run's, much shorter consensus wall times, and a time limit that hits only the
   intrinsically heavy containers. Test: one trio member with `TMPDIR` on local disk, against
   the same member on cephfs.
2. **Replace the wall-clock limit with a deterministic budget.**
   - **Assembly-size gate, decided before any LAST call.** Skip the container, or degrade it,
     when the bp lamassemble would get exceed a threshold: phased, the groups' reads; unphased,
     all reads. Its time grows with the square of that (6.3). In practice, cut-read bp after
     read selection > 300–400 kb is nearly as good (6.2, 6.4). It is the same on every machine
     and for every load.
   - **LAST work budget.** Count the alignments lastal streams out and stop above N. It is
     deterministic, but it spends the work up to the budget.
   - **CPU time instead of wall clock.** Use the RUSAGE of the process and its children. It
     ignores I/O waits and oversubscription but still depends on the CPU speed.
3. **Make the heavy containers cheaper instead of dropping them.**
   - **LAST `-m`.** 50 was chosen so that reads in short tandem repeats get linked. In satellite
     arrays it multiplies LAST's work (6.3). Try `-m 5` or `-m 10` first and redo with 50 only
     when the layout leaves reads unlinked, as is already done for the strand. The consensus
     quality of that change must be checked against T2TQ100.
   - **Read cap in the unphased fallback.** It assembles all reads, while the phased path
     assembles at most 50. `--assembly-max-reads` exists and is off by default.
   - **No escalation for predicted-heavy containers.** `--heavy-container-bp` exists and is
     off. Levels that end in a timeout are discarded: 40–64% of the escalated containers' wall.

## 7. Deterministic budgets on the trio (cluster, budget10, 2026-10-08)

Implemented in 425ba7c (all off by default):

- `--consensus-tmp-dir` (run) / `--tmp-dir` (consensus): each consensus process keeps its
  temporary files, and its tools', in its own directory under this base, removed at its end.
- `--max-assembly-bp`: a container about to assemble more bp of reads in one assembly (after
  `--assembly-max-reads`) is dropped before that assembly starts.
- `--lamassemble-max-initial-matches 10,50`: LAST -m 10, then 50 only when the layout leaves reads
  unlinked (after the one-strand → both-strands retry of each -m).

Setup: evaluation/tuning experiment budget10. The data are the trio on the 10% set (verify10),
64 threads per run. Every variant uses the defaults (accurate + tiered phasing, read selection
factor 3, escalation 1:20,4:60,12:120). `/tmp` is the node-local per-job directory of SLURM.
The report is evaluation/tuning/reports/regions10/budget10.txt; `budget_analysis.py` gives the
per-container part. Earlier replicates put the noise at about ±0.0005 F1 and ±0.0007 MC.

| variant | temporary files | container bound | CPU h (trio) | run wall, sum (min) | T2TQ100 F1 | V5 F1 | MC | dropped |
|---|---|---|---|---|---|---|---|---|
| b_default | `$TMPDIR` (cephfs) | 180 s | 36.9 | 171 | 0.8291 | 0.8254 | 0.9285 | 38 (time) |
| b_private | private dirs on cephfs | 180 s | 38.3 | 169 | 0.8287 | 0.8248 | 0.9285 | 39 (time) |
| b_tmp | /tmp | 180 s | 41.9 | 169 | 0.8288 | 0.8246 | 0.9285 | 38 (time) |
| b_tmp_nolimit | /tmp | none | 65.3 | 202 | 0.8303 | 0.8259 | 0.9279 | 0 |
| b_gate | /tmp | 300 kb gate | 27.2 | 142 | 0.8305 | 0.8254 | 0.9266 | 48 (size) |
| b_m10 | /tmp | none | 25.5 | 123 | 0.8301 | 0.8254 | 0.9288 | 0 |
| b_m10_cap | /tmp, read cap 30 | none | 24.0 | 114 | 0.8297 | 0.8254 | 0.9284 | 0 |
| b_m10_gate | /tmp | 300 kb gate | 19.8 | 107 | 0.8293 | 0.8259 | 0.9272 | 47 (size) |

Non-TRF F1 is identical in every variant (V5 0.9316, T2TQ100 0.9276). Q100 representation of the
HG002 containers vs b_default: b_m10 +0.00007 (21 better, 13 worse), b_m10_cap +0.00028 (38
better, 14 worse; aggregate identity 0.9694 → 0.9729).

- **`-m 10,50` is the lever.** Without dropping anything, CPU is −31% and the summed run wall
  −28% against the current defaults, at equal F1, MC and Q100. Per container (HG002
  thread-seconds), the containers whose largest assembly is > 200 kb cost 5–10× less, ordinary
  ones ~20% less. About 1,160 containers per trio retry at -m 50 (484 retry only the strand at
  -m 50).
- **The gate is deterministic.** Its decision depends on the assembly inputs, not the clock. b_gate
  and b_m10_gate dropped the same containers in HG003 (15) and HG004 (22). In HG002 (11 / 10),
  the one difference is a wall-clock effect: under -m 50, both phased assemblies of container 2954
  (247 / 259 kb) ran into the last level's 120 s tool timeout, and the KMeans fallback then
  wanted to assemble all its reads (570 kb).
- **The gated containers are worth little.** HG002's 11 held 2 TP and 0 FP in the uncapped run
  (Q100 representation 0.867). They cost 12.3 of the 24.2 thread-hours of b_tmp_nolimit, but only
  1.8 of 7.7 after `-m 10,50`. The 300 kb gate costs ~0.002 MC: the trio members lose 11 / 15 / 22
  containers, i.e. different loci.
- **Gate threshold sweep** on the largest assembly per container (b_tmp_nolimit inputs, HG002
  truth): 300 kb loses 2 TP, 400 kb 0 TP, with a more even trio (10 / 11 / 10) and 11.5
  thread-hours saved uncapped (1.7 with -m 10,50).
- **Temporary files on node-local /tmp** cut the summed consensus batch time from 35.5 to 22.6 h
  (median container 6.5 → 2.4 s; the alignment stage 1.1 → 0.1 s). A private directory per
  process on cephfs does not help (4.5 s). The cluster was quiet in this run (median 6.5 s vs
  41 s in verify10). The run wall times did not change, and escalations rose from 276 to 373,
  mostly the phasing all-vs-all at 1 thread / 20 s. Their AVA stage took 16–18 s on cephfs and
  29–31 s on /tmp, and minimap2 itself took < 2 s of CPU per finished call. That points to CPU
  contention within the job once the I/O waits are gone (likely, not proven).
- **The escalation's per-tool timeouts are the remaining wall-clock switch.** They decide when a
  container is redone with more threads, and at the last level they drop an assembly. -m 10,50
  halves the escalations (466 → 204 without a limit).

### 7.1 Confirmation: -m 10,50 + read cap 30, with and without a 400 kb gate (budget10b, budget15)

The same settings on the 10% set (budget10b) and on the disjoint 15% set (budget15, against the
current defaults), at 64 threads, with node-local /tmp and no container time limit:

| set | variant | CPU h | run wall, sum (min) | T2TQ100 F1 | V5 F1 | MC | Q100 aggregate identity | dropped |
|---|---|---|---|---|---|---|---|---|
| 15% | defaults (180 s, `$TMPDIR`) | 38.9 | 170 | 0.8570 | 0.8535 | 0.9329 | 0.9797 | 45 (time) |
| 15% | m10 + cap 30 | 21.4 | 116 | 0.8569 | 0.8535 | 0.9338 | 0.9818 | 0 |
| 15% | m10 + cap 30 + gate 400 kb | 18.1 | 104 | 0.8565 | 0.8535 | 0.9333 | 0.9819 | 21 (3 / 9 / 9) |
| 10% | defaults (180 s, `$TMPDIR`) | 36.9 | 171 | 0.8291 | 0.8254 | 0.9285 | 0.9694 | 38 (time) |
| 10% | m10 + cap 30 | 24.0 | 114 | 0.8297 | 0.8254 | 0.9284 | 0.9729 | 0 |
| 10% | m10 + cap 30 + gate 400 kb | 20.2 | 108 | 0.8300 | 0.8264 | 0.9278 | 0.9727 | 19 (6 / 7 / 6) |

Non-TRF F1 is identical within each set. Q100 container representation on the 15% set vs the
defaults: +0.00008 (27 better, 19 worse) for both m10 variants.

Conclusion:

- -m 10,50 + read cap 30 makes the wall-clock container limit unnecessary. Against the defaults:
  CPU −35% (10%) and −45% (15%), run wall −33% / −32%. F1 is equal, MC equal or better, Q100
  identity better.
- The 400 kb gate on top saves another 14–16% CPU (−45% / −53% in all). It is deterministic, the
  trio members lose similar numbers of containers, and its F1 and MC differences to m10 + cap 30
  are within noise (MC −0.0005 on both sets).
- The remaining wall-clock dependence is the escalation's per-tool timeouts (1:20,4:60,12:120).
  They are now rarer (297 → 121 escalations on the 15% set) but still decide when a container is
  redone with more threads and, at the last level, when an assembly is given up.
