# POA + in-RAM polishing vs. lamassemble — evaluation, 2026-09-22

Written after implementing the reviewer's suggestion (replace `lamassemble`
with POA, abPOA specifically) as `--consensus-method poa` and benchmarking it
against the `lamassemble` and `racon` backends on the shipped `examples/muc1`
dataset. Branch: `worktree-poa-consensus`.

**Verdict: implemented and available, but NOT ready to replace lamassemble.**
It is faster and produces a higher mean per-base identity, yet it loses
variant-supporting reads and flips genotypes on this example. The default
stays `lamassemble`.

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

## Why it is not ready: read coverage becomes lost variant support

`read cov` is the metric that matters, and it is the one POA loses. It shows up
undiluted in the VCF. Same pipeline, same `TC`, only the consensus backend
differs:

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
3. Benchmark 2 skips the spectral clustering. Once genome-wide TRF and
   mononucleotide beds are available locally (they live in the separate
   `bihealth/svirlpool-data` repo, not on this machine), run the full pipeline
   on an HG002 subset so the tiled-cluster case is measured at scale too --
   that is the case that decides this.

Benchmark scripts used are not committed; they live in the session scratchpad
and are straightforward to recreate from the tables above (`bench_clusters.py`
assembles every preserved `reads.*.fasta` with each backend and scores the
result with `mappy`).
