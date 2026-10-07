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
