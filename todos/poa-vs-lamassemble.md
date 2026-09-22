# POA + in-RAM polishing vs. lamassemble — evaluation, 2026-09-22

Written after implementing the reviewer's suggestion (replace `lamassemble`
with POA, abPOA specifically) as `--consensus-method poa` and benchmarking it
against the `lamassemble` and `racon` backends. Branch: `worktree-poa-consensus`.

**Verdict: at parity with lamassemble genome-wide, and better on exactly the
variant class svirlpool exists to call -- at 1.57x the runtime.** On a 10%
subset of GRCh38 (chr1 + chr20), POA edges lamassemble on GIAB F1 (0.7106 vs
0.7084), genotype concordance (0.8438 vs 0.8351) and trio Mendelian consistency
(92.81% vs 92.38%, with fewer violations in absolute terms despite 10% more
calls). Its clearest win, in both evaluations independently, is 350-1000 bp
indels. It costs 1.57x the wall clock and has a higher, much less predictable
memory ceiling.

**The default stays `lamassemble`** -- the quality differences are real but
small, and they do not yet justify changing a default on one tenth of one
genome. What the numbers do justify is treating POA as a serious candidate
rather than an experiment, and running the remaining 90%.

> **Correction, written after the genome-scale run.** An earlier version of this
> document concluded "NOT ready to replace lamassemble", on the strength of the
> muc1 example alone, where POA lost variant-supporting reads and flipped
> genotypes to `0/0`. **That conclusion does not generalise.** muc1 is a single
> pathological VNTR whose clusters are extreme tiling cases; genome-wide, most
> clusters are stacks, where POA matches or beats lamassemble. The muc1 section
> below is kept because the *mechanism* it identifies is real and still bounds
> where POA is weak -- but it is not representative, and sizing a decision on it
> was a mistake. Benchmark 3 is the result to quote.

---

## What was built

`src/svirlpool/localassembly/poa_consensus.py`, selected via
`--consensus-method poa`. It touches no temporary files and spawns no
subprocesses — everything is Python strings through `pyabpoa` and `mappy`
(minimap2's bindings). Only the final consensus FASTA is written, because
`final_consensus()` downstream re-aligns the cut reads against it.

Per cluster:

1. **Orient** every read onto a common strand (POA has no notion of reverse
   complements; a flipped read would be threaded as unrelated sequence).
2. **Two backbone candidates**, scored by `layout_score` — the number of read
   bases the candidate explains — and the better one wins:
   - the abPOA consensus over the cluster (capped at `max_poa_reads`);
   - the read that hosts the cluster best, extended along reads overhanging
     its ends.
3. **Polish** the winner with windowed POA, racon's algorithm: map every read
   with `mappy`, cut the backbone into 500 bp windows, re-run abPOA per window,
   concatenate.

## Two defects found and fixed during development

Both were found by benchmarking, not by reading code, and both are recorded as
regression tests in `tests/test_poa_consensus.py`.

- **A plain abPOA consensus collapses a tiled cluster.** svirlpool's clusters
  are *not* stacks of reads that all span the locus; they are cut reads that
  *tile* it. muc1 cluster 09 is two 20 kb reads plus 13 short reads piled on
  one flank. abPOA's consensus is the heaviest path through the graph, so it
  follows the pile: 508 bp instead of 21 161 bp, i.e. the whole insertion lost.
  Neither the backbone choice (median vs. longest read) nor the alignment mode
  (`g`/`l`) changes this — it is structural. This is why the layout step exists.
- **Windowed polishing silently did nothing** when reads were shorter than the
  window, because fragments were required to span a whole window. In a cut-read
  cluster that is most reads. Fragments are now padded with the backbone's own
  sequence out to the window edges.

Backbone extension initially welded reads on at spurious repeat-phase overlaps,
producing consensuses their own reads only aligned to in pieces (and on the
reverse strand). It is now guarded by a long, high-identity overlap requirement
*and* by `layout_score` having to improve — the first round that does not is
discarded.

## Benchmark 1 — 18 real muc1 clusters, identical read sets

Clustering happens before the backend is chosen, so every backend saw exactly
the same clusters (preserved via `--tmp-dir-path`). These are *cut reads that
tile* a locus.

| backend | clusters | time | mean identity | mean read cov | mapped | total bp |
| --- | --- | --- | --- | --- | --- | --- |
| lamassemble | 18 | 6.98 s | 0.9290 | **0.946** | 0.996 | 52 777 |
| racon | 17 | **0.52 s** | 0.9352 | 0.796 | 0.995 | 51 390 |
| poa + polish | 18 | 3.62 s | **0.9377** | 0.906 | 0.996 | 49 474 |

`identity` is mean gap-compressed identity of the cluster's reads against the
consensus; `read cov` the mean fraction of each read's bases the consensus
accounts for; `mapped` the fraction of reads that align at all. racon fails
outright on the 2-read cluster.

Whole consensus stage on the muc1 example, 8 threads:

| backend | wall clock | peak RSS |
| --- | --- | --- |
| lamassemble | 15.5 s | 230 MB |
| racon | 7.5 s | 230 MB |
| poa + polish | 10.0 s | 593 MB |

## Benchmark 2 — 12 real HG002 20x loci

`HG002.minimap2.20x.softclipped.bam`, GRCh38. Loci were found by scanning
chr1/chr2/chr15 for bins where >=4 reads carry a >=300 bp insertion, keeping
the 12 strongest well-separated ones. Reads spanning the locus +/-3 kb were cut
to that window the way svirlpool cuts them. This skips the spectral clustering,
so a pool may mix haplotypes -- the harder case, and identical for every
backend. Crucially, these pools are **stacks**: every read spans the window.

| backend | loci | time | mean identity | mean read cov | mapped | total bp |
| --- | --- | --- | --- | --- | --- | --- |
| lamassemble | 12 | 17.5 s | 0.9565 | 0.996 | 1.000 | 108 498 |
| racon | 12 | **0.7 s** | 0.9172 | 0.996 | 1.000 | 116 449 |
| poa + polish | 12 | 15.2 s | **0.9604** | **0.997** | 1.000 | 102 868 |

**On stacked pools POA matches or beats lamassemble on every quality axis and
is slightly faster.** The muc1 deficit is specific to *tiled* cut-read
clusters, where the layout -- not the polishing -- is the limiting step. That
is the most useful result here: the reviewer's suggestion is sound wherever
reads span the locus, and the whole difficulty lives in the case where they
do not.

## Benchmark 3 — HG002 over 10% of GRCh38, against GIAB (the headline result)

`--regions` restricted to chr1 + chr20 = 313.4 Mb, **10.15%** of the GRCh38
primary assembly. HG002 20x, GRCh38, full svirlpool pipeline per backend, then
`truvari bench 4.3.1` against `HG002_GRCh38_v5.0q_stvar`, `>=50 bp`, inside the
GIAB benchmark BED intersected with the run regions. 24 threads.

| | lamassemble | poa |
| --- | --- | --- |
| TP / FN / FP | 1674 / 1134 / 244 | 1697 / 1111 / 271 |
| precision | **0.8728** | 0.8623 |
| recall | 0.5962 | **0.6043** |
| **F1** | 0.7084 | **0.7106** |
| **GT concordance** | 0.8351 | **0.8438** |
| wall clock | **13:26** | 18:45 |
| peak RSS | **13.6 GB** | 17.4 GB |

Essentially a tie, with POA marginally ahead on F1 and genotype concordance and
marginally behind on precision. The size breakdown shows the trade is
structured rather than noise:

| svtype | size | lam recall / prec | poa recall / prec | dTP |
| --- | --- | --- | --- | --- |
| INS | 350-1000 | 0.560 / 0.774 | **0.598 / 0.808** | **+14** |
| INS | 1000-10000 | 0.611 / 0.864 | **0.632** / 0.858 | +5 |
| DEL | 1000-10000 | 0.817 / 0.906 | **0.859** / 0.897 | +3 |
| DEL | 50-100 | 0.565 / 0.899 | **0.580 / 0.905** | +6 |
| INS | 100-350 | 0.582 / **0.890** | 0.579 / 0.840 | -2 (+30 FP) |
| INS | 50-100 | **0.539** / 0.895 | 0.527 / 0.881 | -6 |

POA wins on mid-to-large insertions -- better recall *and* precision in the
350-1000 bp bin, which is the class svirlpool exists to call -- and loses on
small insertions, with its extra false positives concentrated almost entirely
in the 100-350 bp INS bin. Both directions match the mechanism: windowed-POA
polishing reconstructs large events better, while the layout occasionally
invents a short insertion.

The timing is the reverse of muc1, where POA was ~2x faster. At genome scale
most clusters are deep stacks, and building both backbone candidates and
ranking them costs more than lamassemble's single pass. If POA is adopted,
skipping the layout candidate when the reads demonstrably co-span is the
obvious optimisation -- the information needed to decide is already computed.

## Benchmark 4 — Mendelian consistency, HG002/HG003/HG004 trio

Same 10% subset, all three samples run through both backends, joint-called per
backend, then `svirlpool.analysis.mendelian_consistency` with `--min-size 50`.
HG002 is the child. Genotype consistency is measured over *informative* sites
only (consistent + inconsistent); `non_informative` and `no_call` are reported
separately rather than folded in.

| | lamassemble | poa |
| --- | --- | --- |
| trio records | 3 600 | 3 973 |
| informative (denominator) | 3 032 | 3 116 |
| **consistent** | 2 801 — 92.38% | 2 892 — **92.81%** |
| inconsistent | 231 — 7.62% | **224** — **7.19%** |
| non-informative | **329** | 549 |
| no-call rate | **6.64%** | 7.75% |

POA calls 10% more variants yet produces *fewer* Mendelian violations in
absolute terms (224 vs 231), so the extra calls are not noise. It pays for that
with more no-calls and markedly more non-informative sites.

By type and size (percent consistent, informative sites):

| stratum | lamassemble | poa |
| --- | --- | --- |
| DEL (all) | 91.10 | **91.82** |
| INS (all) | 93.08 | **93.39** |
| <100 bp | 94.23 | **94.34** |
| 100-350 bp | **94.32** | 94.20 |
| **350-1000 bp** | 86.99 | **89.43** |
| 1k-10k bp | 92.06 | **92.14** |
| >10k bp | **72.73** (16/22) | 65.38 (17/26) |

**The 350-1000 bp bin is the result worth keeping.** Two independent
evaluations -- a truth-set comparison against GIAB and a pedigree-based
consistency check that uses no truth set at all -- both single out the same
stratum as POA's clearest gain: +14 true positives with better precision *and*
better recall in Benchmark 3, and +2.4 points of Mendelian consistency here.
Agreement between a truth-set metric and a truth-free one is much harder to
explain away as a benchmark artefact than either result alone.

The >10k bin moves the other way, but on 22 vs 26 informative sites it settles
nothing; a wider run is needed before reading anything into it.

## Runtime over the 10% subset (24 threads)

| sample | lamassemble | poa |
| --- | --- | --- |
| HG002 | 13:26 / 13.6 GB | 18:45 / 17.4 GB |
| HG003 | 12:27 / 12.8 GB | 20:57 / 13.3 GB |
| HG004 | 11:33 / 13.2 GB | 19:03 / 22.4 GB |
| **total** | **37:26** | 58:45 |

POA costs **1.57x the wall clock**, and its peak memory is both higher and far
more variable (13.3-22.4 GB against a steady 12.8-13.6 GB). The variance is the
more awkward property for cluster scheduling: lamassemble's ceiling is
predictable, POA's is not.

## The muc1 failure mode: read coverage becomes lost variant support

On muc1, `read cov` is the metric POA loses, and it shows up undiluted in the
VCF. Same pipeline, same `TC`, only the consensus backend differs. Read this as
a description of POA's *worst case*, not as its expected behaviour -- Benchmark
3 shows the genome-wide picture is a wash:

| locus | lamassemble | poa |
| --- | --- | --- |
| 93931 / 94010 | DV=17/21 → **1/1** | DV=4/21 → 0/0 |
| 94431 / 94425 | DV=19/22 → **1/1** | DV=4/22 → 0/0 |
| 109709 / 109710 | DV=14/16 → **1/1** | DV=2/16 → 0/0 |
| 110272 / 110222 | DV=15/18 → **1/1** | DV=12/18 → 0/1 |
| 248259 | DV=23/26 → **0/1** | split into 347 bp + 82 bp, DV=0 → 0/0 |
| 278587 | DV=10/11 → **1/1** | absent |
| 117941 / 117930 | DV=19/23 → 1/1 | DV=19/23 → 1/1 (agrees) |

(Re-verified against the final implementation; the picture is unchanged.)

`TC` is unchanged throughout, so the clusters and the coverage are the same:
it is `DV` that collapses. Reads that align to only part of the consensus do
not contribute their SV signal, so they never count as supporting. lamassemble
also calls one insertion per locus where POA splits it in two.

## Where the remaining gap actually is

Not in the polishing — in the **layout**. The clusters where POA still trails
(muc1 08, 13, 14, 15) are short, noisy, repeat-rich fragments that minimap2's
`map-ont` cannot place on the backbone but LAST, with lamassemble's trained
`promethion.mat`, can. Making the minimap2 index more sensitive was tried
(k=13/w=5, k=11/w=5, k=11/w=3) and moves mean coverage by 0.004 while costing
`mapped` — it is not the answer.

lamassemble solves these by a genuine multiple alignment over LAST overlaps: on
cluster 15 it merges the short-read pile and the long-read region into one
1018 bp consensus that explains more read bases than any single-backbone
solution the current search can reach. Closing this needs a real
overlap-layout-consensus step, not a better POA.

## If this is picked up again

1. The objective (`layout_score`) is right and already in place — the *search*
   over layouts is what is too weak. A proper OLC over `mappy` all-vs-all
   overlaps, scored by `layout_score`, is the natural next step.
2. Alternatively keep lamassemble for the layout and use the in-RAM windowed
   POA purely as its polisher; identity is the axis where POA already wins on
   both datasets (0.9377 vs 0.9290 on muc1, 0.9604 vs 0.9565 on HG002).
3. **Close the runtime gap before the wider run.** POA builds both backbone
   candidates unconditionally. At genome scale most clusters are deep stacks,
   where the layout candidate cannot win -- and `layout_score` already computes
   what is needed to know that. Skipping it when the reads demonstrably co-span
   should recover most of the 1.57x without touching quality, and would also
   cut the memory variance.
4. **Run the remaining 90%.** One tenth of one genome, on differences this
   small, is not enough to move a default. The trio pipeline is scripted end to
   end in `/home/mayv_c/svirlpool-poa-bench/drive_all.sh`; pointing
   `subset10.bed` at the rest of the assembly is the only change needed.
5. The >10k stratum disagrees with every other bin (POA worse: 65.4% vs 72.7%
   Mendelian consistency) on 22 vs 26 informative sites. Either it is noise or
   it is the tiled-cluster weakness from muc1 showing through at the largest
   events -- worth resolving, since that is where the layout argument predicts
   POA *should* struggle.

## Reproducing this

The genome-scale pipeline (Benchmarks 3 and 4) is scripted in
`/home/mayv_c/svirlpool-poa-bench/`, outside the repo because it carries
multi-gigabyte working directories:

| file | what it does |
| --- | --- |
| `drive_all.sh` | the whole comparison: both backends x 3 samples, truvari, joint calling, Mendelian |
| `run_sample.sh` | one sample, one backend, over a regions bed |
| `eval_sv.sh` | truvari 4.3.1 against GIAB, confident regions = GIAB bed ∩ run regions |
| `stratify.py` | the size/type breakdown in Benchmark 3 |
| `subset10.bed` | chr1 + chr20; point this elsewhere to widen the run |
| `trio.ped` | HG002 child, HG003 father, HG004 mother |
| `benchenv/` | pixi env holding truvari (deliberately not a project dependency) |

truvari is **not** added to `pyproject.toml` -- it is an evaluation tool, not a
runtime dependency, and it drags in a conflicting htslib stack.

The cluster-level harnesses for Benchmarks 1 and 2 were not kept; they are
straightforward to recreate from the tables (assemble every preserved
`reads.*.fasta` with each backend, score with `mappy`).
