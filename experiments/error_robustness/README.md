# Robustness of the phased consensus to read errors (exploratory)

How much additional per-base read error the consensus stage (phased mode,
`perf/consensus-profiling` 0dcbc04) tolerates on the svp_improvements HG002
ONT 20x data. Exploratory; a formal test follows later.

## Setup

* `consensus.add_read_errors` / `added_error_rate` (threaded
  `crs_containers_to_consensus` -> `process_consensus_container` ->
  `trim_reads`): every base of the cut reads is an error with probability
  `rate`, one third each substitution, insertion of a random base, deletion.
  Errors are drawn per (read name, position on the full read), so the reads cut
  for phasing (CR +-10 kb) and for the assembly (CR +-200 bp) carry the same
  errors where they overlap, and the errors of a lower rate are a subset of
  those of a higher one. Checked: realised divergence 0.98 / 4.8 / 9.2 % at
  1 / 5 / 10 %. Not mutated: the consensus padding (raw read flanks) and the
  read-vs-reference alignments (KMeans features; unused in phased mode).
* Native read divergence of the data: median 1.1 % (minimap2 `de`,
  gap-compressed; p10 0.6 %, p90 4.5 %), so the added rates come on top of ~1 %.
* 200 containers (10 % of the 2004 of `esc2_phased`, random, seed 0,
  `results/crIDs.txt`); 190 evaluable against T2TQ100 (in the benchmark
  regions; 77 het, 39 cpx_het, 63 hom, 11 ref).
* `run_rates.py`: consensus stage with the run's config (phased, escalation
  1:20,4:60,12:120) at added rates 0, 0.01, ..., 0.10; `eval_rates.py`:
  each consensus aligned to its locus +-10 kb (minimap2 map-ont), net indel
  over the locus span matched to the T2TQ100 haplotypes (as
  `ava_phasing/container_truth.py`); trio pair accuracy of the read groups.

## Result: the phasing breaks at ~5-6 % added error, the assembly does not

| added error | phased | pair acc | alleles recovered | het recovered | cpx_het | hom | non-TRF cont. | TRF cont. | consensus diff % | CPU s | escalated |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 180 | 0.941 | 0.882 | 0.909 | 0.590 | 0.921 | 0.982 | 0.791 | 1.02 | 280 | 0 |
| 0.01 | 181 | 0.930 | 0.886 | 0.909 | 0.590 | 0.921 | 0.982 | 0.791 | 1.20 | 311 | 0 |
| 0.02 | 176 | 0.910 | 0.869 | 0.857 | 0.590 | 0.921 | 0.964 | 0.769 | 1.25 | 312 | 0 |
| 0.03 | 173 | 0.879 | 0.856 | 0.857 | 0.513 | 0.921 | 0.929 | 0.761 | 1.45 | 310 | 0 |
| 0.04 | 154 | 0.815 | 0.837 | 0.818 | 0.538 | 0.921 | 0.911 | 0.754 | 1.86 | 316 | 0 |
| 0.05 | 133 | 0.738 | 0.801 | 0.766 | 0.436 | 0.921 | 0.857 | 0.716 | 2.63 | 313 | 0 |
| 0.06 | 97 | 0.665 | 0.696 | 0.584 | 0.359 | 0.905 | 0.714 | 0.634 | 3.15 | 324 | 0 |
| 0.07 | 61 | 0.578 | 0.569 | 0.364 | 0.103 | 0.921 | 0.607 | 0.485 | 3.68 | 350 | 0 |
| 0.08 | 28 | 0.536 | 0.477 | 0.169 | 0.026 | 0.905 | 0.536 | 0.366 | 4.24 | 366 | 0 |
| 0.09 | 9 | 0.512 | 0.418 | 0.065 | 0.000 | 0.905 | 0.482 | 0.313 | 4.58 | 406 | 0 |
| 0.10 | 1 | 0.506 | 0.376 | 0.000 | 0.000 | 0.873 | 0.446 | 0.269 | 5.07 | 438 | 0 |

phased: containers (of 200) with >= 2 alleles; pair acc: trio pair accuracy
of the reads' consensus membership; recovered columns: fraction of evaluable
containers with every truth allele matched by a consensus (alleles
recovered: of all truth alleles); consensus diff: mismatches + indels < 5 bp
per aligned reference bp of the consensuses. Full table
`results/eval.summary.tsv`, per container `results/eval.tsv`.

* Up to +3 % nothing much changes (het 0.91 -> 0.86, 95 % of containers
  recover the same alleles as without added errors). +4-5 % is the knee;
  above +6 % (~7 % total) het alleles are lost fast, at +10 % none.
* hom and ref loci keep their allele (hom 0.92 -> 0.87): lamassemble still
  assembles the SV from noisy reads; its consensus gets less accurate (1.0 ->
  5.1 % small differences) but the SV size holds. So what breaks is the
  allele separation, not the assembly.
* Cost stays flat: +56 % CPU at +10 %, no timeouts, no escalations.
* When the phasing still phases, it is right: phasing alone (`phase_diag.py`),
  trio accuracy of phased het containers 0.997 / 0.980 / 0.961 / 0.971 at
  0 / 4 / 6 / 8 %. It fails by giving up (status `no_information`, all reads
  unassigned: 9 -> 193 containers), the fallback is one consensus.

## Mechanism

`phase_steps.py` on container 1558 (het, 22 reads):

| added | SNV sites | split along trio haplotypes | mixed | W within hap (mean, <0) | clusters |
|---|---|---|---|---|---|
| 0 | 892 | 762 | 130 | +110, 0 % | 12, 8 |
| 0.04 | 1857 | 581 | 1276 | +65, 8 % | 12, 8, 2 |
| 0.06 | 2912 | 556 | 2356 | +18, 19 % | 11, 8, 2, 1 |
| 0.08 | 4371 | 508 | 3863 | -35, 74 % | 6 + singletons -> one allele |
| 0.10 | 6487 | 447 | 6040 | -104, 91 % | all singletons -> no_information |

Sequencing errors become SNV sites and the recurrence filter does not remove
them (4402 -> 4371 at 8 %): with thousands of candidate sites per target
read, some other site almost always matches a noise split at >= 90 %
concordance. Their disagreements pull the pair weights of reads of the
*same* haplotype negative, until the correlation clustering merges nothing.
Between haplotypes the weights stay strongly negative throughout: the
haplotype signal is still there, it is drowned.

Relaxing single filters does not help (`phase_diag.py`, fraction of het
containers phased at 0 / 4 / 6 / 8 %): base 0.98 / 0.86 / 0.60 / 0.17;
recurrence concordance 0.8: 0.98 / 0.87 / 0.57 / 0.16; `min_alt` 5:
0.96 / 0.83 / 0.54 / 0.17; recurrence 2: 0.98 / 0.88 / 0.60 / 0.21; no indel
masking: 0.97 / 0.46 / 0.11 / 0.01 (the masking protects).

Candidate directions (not tried): a recurrence test that compares the minor
read sets of two sites (not the concordance over all reads, which the
majority dominates) or corrects for the number of candidate sites; weighting
disagreements by site balance / support like the agreements; estimating
`error_rate` of the pair weights from the observed read divergence instead of
the fixed 0.1.

## Legacy clustering on the same containers

Same 200 containers, same work dir (`esc2_phased`), `run_rates.py --mode
legacy` (KMeans on summed indels, spectral fallback, -U 25,35). Side by side
(`compare.py`, `results/compare.tsv`; recovery = containers with every truth
allele matched):

| added | het phased | het legacy | cpx_het phased | cpx_het legacy | hom phased | hom legacy | alleles phased | alleles legacy | pair acc phased | pair acc legacy | CPU phased | CPU legacy |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 0.909 | 0.753 | 0.590 | 0.359 | 0.921 | 0.905 | 0.882 | 0.778 | 0.941 | 0.728 | 280 | 144 |
| 0.02 | 0.857 | 0.766 | 0.590 | 0.359 | 0.921 | 0.905 | 0.869 | 0.791 | 0.910 | 0.730 | 312 | 153 |
| 0.04 | 0.818 | 0.766 | 0.538 | 0.308 | 0.921 | 0.905 | 0.837 | 0.784 | 0.815 | 0.722 | 316 | 156 |
| 0.05 | 0.766 | 0.753 | 0.436 | 0.359 | 0.921 | 0.905 | 0.801 | 0.781 | 0.738 | 0.716 | 313 | 158 |
| 0.06 | 0.584 | 0.753 | 0.359 | 0.333 | 0.905 | 0.905 | 0.696 | 0.788 | 0.665 | 0.714 | 324 | 159 |
| 0.07 | 0.364 | 0.779 | 0.103 | 0.333 | 0.921 | 0.905 | 0.569 | 0.797 | 0.578 | 0.716 | 350 | 162 |
| 0.08 | 0.169 | 0.766 | 0.026 | 0.333 | 0.905 | 0.889 | 0.477 | 0.784 | 0.536 | 0.713 | 366 | 168 |
| 0.10 | 0.000 | 0.753 | 0.000 | 0.359 | 0.873 | 0.905 | 0.376 | 0.788 | 0.506 | 0.702 | 438 | 181 |

Legacy splits 105 of 200 containers into >= 2 consensuses at every rate
(phased: 179 -> 1); 92 % of its containers recover the same alleles at +10 %
as at 0. Phased is ahead up to +5 % (the crossover in alleles recovered is
between 5 and 6 %), legacy from +6 %.

**Caveat, the comparison favours legacy:** legacy's allele decision (KMeans
on the summed indels >= 12 bp per read) is computed from the read-vs-reference
alignments of the BAM, which the simulation does not touch; only its
consensus reads and its spectral fallback see the errors. Phased decides on
the mutated reads. So legacy's flat curve shows that its decision never saw
the errors, not that it would survive real noisier reads. A fair test
realigns the mutated reads to the reference (or mutates the BAM) before the
indel features are taken.

## Files

`run_rates.py` consensus runs per rate (`--mode` overrides the clustering
mode), `eval_rates.py` evaluation, `compare.py` modes side by side,
`breakdown.py` per category + phasing metadata, `phase_diag.py` phasing
alone with filter variants on dumped read sets
(`consensus_perf/dump_phasing_reads.py --orient`), `phase_steps.py` one
container step by step, `py.sh` runner, `results/`.

```
bash py.sh run_rates.py <svirlpool workdir> <out> --procs 12   # ~7 min
bash py.sh eval_rates.py <out> <out>/eval.tsv
bash py.sh breakdown.py <out>
bash py.sh compare.py phased=<out>/eval.tsv legacy=<out_legacy>/eval.tsv --out cmp.tsv
```
