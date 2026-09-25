# Read phasing in the all-vs-all alignments (experimental arm)

Branch `exp/ava-read-phasing`, on top of `fix/all-fixes-integration` 93c0f79
(= flank200 + haplotype-aware merge by default). HG002 ONT 2025.01 20x, the
svp_improvements 5% genome (2004 candidate-region containers).

**Question.** The consensus module splits a container's reads into alleles
from their SV-scale signals only, with k = local copy number. Reads of one
haplotype also share every small variant (SNVs) over their whole length, and
the reads are far longer than the candidate region. Can phasing the reads in
the all-vs-all (AVA) alignments separate alleles better, and give the number
of alleles without a copy-number input?

**Answer: yes.** Phasing on reads cut to CR +-10 kb phases the assigned reads of
96% of containers without a single error (trio truth), the allele count comes out of the
phasing, and end to end HG002 F1 goes up everywhere, most outside tandem
repeats (+2.2 / +3.3 points V5 / T2TQ100).

## Method (`src/svirlpool/localassembly/read_phasing.py`)

1. Reads of the container, cut to the CRs +- `--phasing-flank` (10 kb), are
   aligned all-vs-all with `minimap2 --eqx` *without* `-X` (ava presets imply
   it), so every read is the target of all others; best alignment per pair.
2. **SNV sites** per target read t: a column where >= 3 reads share one
   alternative base and >= 2 carry t's base (minor side >= 20%), outside
   homopolymers >= 4 and >= 2 bp from indels in that read.
3. **Recurrence**: a site is kept only if its read split recurs (>= 90%
   concordant) at another position of t. All het SNVs of a haplotype split the
   reads the same way; coinciding sequencing errors form random subsets.
4. **SV sites**: windows of t with a net indel >= 20 bp or an interior clip in
   some read; reads spanning them disagree (|net| >= 20) or agree (< 10).
5. **Low-quality reads** (median divergence to their partners > 3x the pool's
   and > +5 points) are left out of phasing: a few failed reads (10-23%
   divergence against GRCh38) share wrong bases often enough to pose as a
   recurrent haplotype.
6. Pair weights = log((1-e)/e) * (agree - disagree), where an agreement is
   weighted by the site's balance (at a 2-of-20 site every read "agrees" with
   t regardless of haplotype). Greedy **correlation clustering** (merge while
   the summed inter-cluster weight is positive) -- no k needed.
7. **Allele count**: clusters >= 3 reads not separated by >= 2
   *discriminating* positions (both near-unanimous, different alleles) are
   merged. The remaining clusters are the alleles (capped at 4).

`--consensus-clustering-mode phased` (consensus CLI: `--clustering-mode`)
builds one consensus per allele from the reads cut as before (CR +- 200 bp);
unassigned reads are re-added by alignment as before. With < 2 alleles:
`--phasing-fallback single` (default) builds one consensus from all good reads,
`legacy` runs the old clustering. Each consensus carries the phasing in
`clustering_meta_data` (status, n_alleles, haplotype, site counts). The CN
track is only used for the existing skip-if-CN>4 gate.

## Truth used

* `data/trio_read_labels.tsv.gz` -- HG002 read haplotype (pat/mat) from
  trio-informative SNVs in HG003/HG004 (`trio_read_labels.py`), independent of
  any SV signal. 54,534 sites, 77% of reads labelled, 0.77% minority votes,
  99.97% agreement with phased GIAB v5 deletions.
* `data/container_truth.tsv` -- per container: T2TQ100/V5 per-haplotype
  alleles (het / hom / cpx_het / ref), legacy consensus counts, calls, truvari
  status (`container_truth.py`); `data/test_set.tsv` a stratified set of 58.

## Results

### Read level (446 containers, `batch_eval.py` + `summarize.py`/`allele_eval.py`)

| phasing flank | pair acc. (assigned reads) | containers perfect | unphasable | s / container |
|---|---|---|---|---|
| 2 kb | 0.955 | 0.80 | 46 | 0.8 |
| 5 kb | 0.977 | 0.89 | 21 | 0.9 |
| **10 kb** | **0.987** | **0.96** | **2** | **1.2** |
| 20 kb | 0.991 | 0.97 | 0 | 2.6 |

Allele level (all 2004 containers, 5 kb, early version; T2TQ100 categories,
trio-labelled reads): allele recovered (every truth allele is the majority of a
cluster) new vs legacy -- cpx_het non-TRF 0.84 vs 0.41, cpx_het TRF 0.96 vs
0.78, het TRF 0.98 vs 0.87, het non-TRF 0.98 vs 0.92; cluster purity 0.996 vs
0.944. These are the legacy failure modes found in `container_truth.tsv`:
small (20-49 bp) hets in repeats not split, compound hets of similar size
merged, hom SVs split off into a tiny second cluster.

Signal sources (100 containers, 5 kb): SNV only 0.900 mean pair accuracy,
SV sites only 0.776, both 0.906, legacy 0.761 (haplotype metric).

### End to end (svp_improvements, truvari refined, HG002)

| variant | V5 all | V5 non-TRF | T2TQ100 all | T2TQ100 non-TRF | GT conc. V5 / T2T | wall |
|---|---|---|---|---|---|---|
| flank200 (default) | 0.8565 | 0.9287 | 0.8538 | 0.9169 | 0.879 / 0.860 | 5.5 min |
| flank200 + override subset | 0.8552 | 0.9287 | 0.8585 | 0.9249 | 0.913 / 0.849 | |
| ava_phase (5 kb, v1) | 0.8570 | 0.9489 | 0.8629 | 0.9486 | 0.934 / 0.857 | |
| **ava_phase2 (10 kb)** | **0.8619** | **0.9503** | **0.8685** | **0.9498** | 0.931 / 0.863 | 11.5 min |
| ava_phase, override all | 0.8581 | 0.9489 | 0.8369 | 0.9171 | 0.752 / 0.834 | |

Phased mode needs `sv-calling --multi-assembly-override subset`: a hom SV is
now carried by both haplotype consensuses, and the default override turns it
into 0/1 (T2TQ100 FN 357 -> 424). Outside TRF the FPs halve (33 -> 16). Inside
TRF the clustering improves as much, but calls barely move: there the
bottleneck is the representation of the consensus alignment (see the
svp_improvements README, "allele representation in repeats").

Variants: `svp_variants.yaml` (`bash run.sh stage_benchmark --configfile
<it>`), tables `results/e2e_summary.tsv`.

### Repeat gate -- tried and dropped

A switch that took the phasing arm only for containers in which <= 50% of the
SV signals lie in a tandem repeat (`repeatID` from the TRF annotation; commit
a17bae7, removed again). Benchmarked as `ava_phase3`:

| variant | V5 all | V5 non-TRF | T2TQ100 all | T2TQ100 non-TRF | consensus CPU |
|---|---|---|---|---|---|
| legacy + override subset | 0.8552 | 0.9287 | 0.8585 | 0.9249 | 2134 s |
| phased everywhere (ava_phase2) | **0.8619** | **0.9503** | **0.8685** | **0.9498** | 6838 s |
| phased where repeat fraction <= 0.5 | 0.8561 | 0.9312 | 0.8611 | 0.9297 | 3330 s |

It halves the cost but gives back most of the gain, also in truvari's
non-TRF stratum: 20 of the 22 unrefined V5 non-TRF FPs that phasing removes
come from containers whose signals lie in repeats (e.g. 13 over-split DELs
around chr4:7.86 Mb, container 1303). Truvari's stratum says where the call
lands, not whether the container is repetitive.

*Correction (`../consensus_perf/README.md`):* in ava_phase2 the all-vs-all of
container 1303 hit the 20 s timeout, so the container got one consensus from
the fallback. Its 13 FPs went away through the timeout, not through the
phasing. With the faster phasing it completes (phased correctly), and they
are back.

### Copy number

Phasing counts alleles from the reads: 1711/2004 containers give 2 alleles,
149 one, 92 no information (mostly hom loci without het SNVs, where one
consensus is right), 45 >= 3. Where the CN track says 4 but depth is diploid
(crIDs 791-794, 24-33 reads) phasing finds two balanced haplotypes; where
depth is 4x (790, 80 reads) it finds 6 groups.

## Open ends

* **Runtime**: consensus stage 3.2x CPU (2134 -> 6838 s), longest batch 141 ->
  436 s. Phasing is ~1.2 s per container; the tail is repeat containers where
  AVA hits the 20 s timeout. Addressed on `perf/consensus-profiling`
  (`../consensus_perf`): 6744 -> 4617 s.
* **Genotype from the phasing**: with haplotype consensuses the GT is the
  number of haplotypes carrying a call; `clustering_meta_data` has what
  sv-calling needs to do that instead of read counting + override.
* **Phase blocks**: reads shared between neighbouring containers link their
  haplotypes -> phased VCF (PS).
* Over-splitting in segdup-like / collapsed regions (>= 3 groups in ~2%).
* Near-identical alleles in repeats with no het SNV within 10 kb stay merged.
* No trio / Mendelian run of the phased mode yet.

## Files

`proto.py` prototype (same algorithm, instrumented), `batch_eval.py` runner
(`batch_eval.py <crIDs> <out.tsv> <flank> key=value ...`), `summarize.py`
read-level summary, `allele_eval.py` allele-level evaluation, `results/`
tables (`sub_s1_*`, `sub_s2_*`: sweeps; `all_f5000_both`, `b0_*`: early
versions without the low-quality filter / refinement). `data/flank200/`
(ignored) is a snapshot of the flank200 container DB the scripts read.
