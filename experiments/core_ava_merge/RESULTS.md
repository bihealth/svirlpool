# Core-vs-core alignment merge: results (prototype)

GIAB HG002 (son) / HG003 / HG004, 30x ONT, hg38. Local copies of the evaluation
DBs (built by 8d5ebd7), `bench.sh` (strict Mendel scorer pinned to c8804b1,
truvari 5.4.0 on HG002's column, no refine), `--threads 8`. All arms ran on a
machine shared with two other agents, so wall times are noisy; the merge step's
own time is logged separately. Scratch and per-arm analyses:
`/home/mayv_c/development/svp_mendel/proto/`, arms in `.../svp_mendel/arms/`.
Scripts: `scripts/` next to this file (hard-coded local paths).

## Arms

| arm | commit | flags | PASS | cons | incons | MC | P (T2T) | R (T2T) | F1 (T2T) | H (T2T) | F1 (V5) | H (V5) | wall |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| base_c8804b1 | c8804b1 | - | 40,696 | 35,218 | 3,204 | 0.9166 | 0.8755 | 0.7649 | 0.8165 | 0.8473 | 0.8163 | 0.8472 | 16:31 |
| proto_align | f89e59b | align (v1: window around x only) | 39,637 | 34,831 | 2,505 | 0.9329 | 0.8854 | 0.7617 | 0.8189 | 0.8537 | 0.8180 | 0.8530 | 20:54 |
| proto_align_v2 | ca3fab9 | align (hull of both SVs) | 40,608 | 35,574 | 2,582 | 0.9323 | 0.8863 | 0.7678 | 0.8228 | 0.8563 | 0.8223 | 0.8560 | 21:10 |
| proto_or_v2 | ca3fab9 | or | 38,788 | 34,248 | 2,386 | 0.9349 | 0.8800 | 0.7494 | 0.8094 | 0.8473 | 0.8087 | 0.8468 | 20:59 |
| proto_align_v3 | f4768cd | align (+ twins) | 39,849 | 35,071 | 2,408 | 0.9358 | 0.8884 | 0.7623 | 0.8206 | 0.8557 | 0.8199 | 0.8552 | 21:50 |
| proto_align_v4 | b626531 | align (+ transitive gate, smaller-size tol.) | 41,327 | 36,098 | 2,584 | 0.9332 | 0.8860 | 0.7722 | 0.8252 | 0.8583 | 0.8249 | 0.8581 | 21:54 |
| proto_align_v4_tol02 | b626531 | align, `--core-alignment-frac-tol 0.2` | 39,902 | 34,956 | 2,467 | 0.9341 | 0.8855 | 0.7719 | 0.8248 | 0.8583 | 0.8242 | 0.8578 | 22:00 |
| **proto_align_v5** | **0f3908f** | **align (+ net indels < 50 bp apart)** | 41,019 | 35,896 | 2,554 | **0.9336** | 0.8861 | 0.7719 | **0.8251** | **0.8583** | 0.8248 | 0.8581 | 22:04 |
| proto_align_v6 | 14881fe | align (+ same-consensus guard, guarded union for INV) | 41,449 | 36,044 | 2,673 | 0.9310 | 0.8893 | 0.7755 | 0.8285 | 0.8600 | 0.8283 | 0.8599 | 18:52 |

v5 vs base: MC +1.70 pt (3,204 -> 2,554 inconsistent, consistent +678), precision
+1.06 pt, recall +0.70 pt, F1 +0.86 pt, H +1.10 pt (T2T; V5 the same). truvari
T2T: TP-base 21,561 -> 21,759, FP 2,816 -> 2,593.
For comparison, the D3 size-tolerance arms of the parallel agent (same bench):
tol 0.15 MC 0.9204 / F1 0.8118; tol 0.3 MC 0.9296 / F1 0.7917 (they lose recall).
`or` mode keeps the old rule's over-merges and loses recall (H = base): the
alignment has to *replace* the pairwise rule where it can judge, not complement it.

## Cost (v5; base in brackets)

- Cores: 136,136 cores (261 Mb) of consensuses with INS/DEL composites, streamed
  from the three DBs in 6-8 s.
- Vertical merge step: 278 s wall at 8 workers [45 s]. minimap2 time summed over
  workers 1,417 s. 31,719 candidate groups aligned, 95.8k minimap2 processes (one
  per target core), 223,904 consensus pairs aligned, 955,945 composite pairs
  judged (554k evaluable, 58%). 14-18 groups per run hit the 10 s per-group cap and
  34 groups (> 300 kb of cores) are skipped (fallback to the pairwise rule).
- Whole call: user CPU 2,581 s [1,161 s], peak RSS 20.31 GB [20.07 GB], wall
  22:04 [16:31; the merge accounts for ~3:50, the rest is machine load].
- Nondeterminism: the wall-clock cap makes the set of aligned groups depend on load
  (evaluable pairs varied 539k-555k between runs).

## Pair links vs the existing rule (v5, genome-wide, `CORE_ALN_MERGE|SUMMARY`)

| | pairs |
|---|---:|
| alignment link (corresponding places, haplotypes agree) | 191,385 |
| of these also linked by the old rule | 170,762 |
| alignment link the old rule rejected (added) | 20,623 |
| old-rule link the alignment rejected (removed, evaluable pairs) | 33,616 |
| twins (same haplotype difference placed differently) | 6,040 |
| not evaluable -> old rule linked (fallback) | 25,130 |
| links skipped by the transitive haplotype gate | 7,756 |

Records (harness on the dumped composites, v5 code, records >= 50 bp):
47,211 -> 47,601; cross-sample pattern pairs in one record 144,793 -> 121,411
(12,625 added, 36,007 removed).

## Mendel transitions (record level, v5 vs base)

- base inconsistent -> consistent: 904; still inconsistent: 2,219.
- new violations (arm inconsistent, base counterpart consistent): 298
  - 217 split: the old record held more consensuses (the alignment or the gate
    took one out). Typical: chr3:194,527,040, child consensus .0 (+74) differs from
    the child's .1 and both parents by a 10-12 bp indel (+4 bp) in the window;
    the old rule merged all four (1/1 everywhere); now .0 stands alone. Also VNTRs
    cut into 4-5 pieces per haplotype whose haplotypes differ (chr8:19,461,780):
    the old rule paired equal-size pieces, consistent by accident.
  - 61 joined, 16 regrouped, 4 same consensuses.

## tr_substring labels (base seed records; matched by CHROM/POS/SVTYPE/SVLEN, else
## by the child's consensus within 2 kb)

Size-discordant child-ALT records (class a). "joined" = the arm record holds a
consensus of the transmitting parent (Q100 PAT|MAT); "new" = not in the base
record; "wrong only" = newly joined to the non-transmitting parent and not to
the transmitting one.

| class | n | consistent v5 | joined v5 (new) | wrong only | consistent v4_tol02 | offline containment rule (README) |
|---|---:|---:|---:|---:|---:|---|
| A1 (same Q100 haplotype in child and transmitting parent) | 592 | **410 (69%)** | 435 (386) | 2 | 421 | 400 fixed (68%) |
| A2 | 50 | 14 | 15 (12) | 0 | 20 | 26 merged |
| B (candidate de novo) | 21 | **0** | 4 (2) | 0 | 3 | 15 merged (71%) |
| C (noise ~ difference) | 44 | 15 | 16 (14) | 1 | 20 | 34 merged |
| D_cons (child assembly wrong) | 96 | 35 | 22 (21) | 1 | 37 | 35 merged (36%) |
| D_gt | 52 | 22 | 5 (5) | 0 | 28 | 28 merged |
| controls A1 (class c) | 254 | 250 (base 254) | | | | |

Wrong-parent attachments of the offline containment rule: 24% of its merges.
Here 2 A1 records are newly attached only to the non-transmitting parent (79 are
newly joined to both parents, i.e. both carry an allele the alignment judges the
same).

Pair-level reachability (v4 code, `scripts/reachable.py`): for 579 of the 592
A1 records a transmitting-parent composite is in the child's candidate group;
the best child-parent verdict is link 495, twin 44, distinct 13, unlinked 7,
no alignment 16.

What keeps A1 records inconsistent although the pair links (v2 analysis,
`scripts/a1_gap.py`): the parent is in the record but genotyped 0/0 (51; DV
3-7 against DR 20-30, the D2 problem), another piece of the child's consensus
is the one that got joined (59; decomposition into several records), the
record moved > 2 kb or was lost (16).

## Case-level evidence (v4/v5 code, `scripts/loci.py`)

Fixed:
- chr1:108,023,172 (A1, 2-bp STR): child INS 76 joins the mother's 76+13 (MAT);
  the old rule kept 76 vs 89 apart (15%).
- chr13:73,209,280 (A1, plan D1 repeatID example): child .0 pairs with the
  mother's .0 piece by piece (256, 46, 66+13), child .1 with the father's .0
  (89+13+100; 74 with 46). The old rule linked the father's 74 to the child's 79
  through the shared repeatID.
- chr10:103,251,035 (A1, 18-bp VNTR, PAT +2,198 cut into 1,326+701+75+96): the
  child's PAT pieces join the father's 2,048; the child's MAT 1,020 joins the
  mother's 1,020; the father's other consensus 2,136 joins the mother's 2,172.
  v1 (window around x alone) had put both child alleles into one record.
- chr10:20,546,130 (README example of a wrong-parent containment merge): the
  child's MAT 1,113 joins only the mother's 1,113; the father's 1,033 (old rule:
  merged into the child's MAT record) stands alone. The child's PAT 646 stays
  with the father's 646, but the father is still genotyped 0/0 there (D2).
- Synthetic (tests): one haplotype decomposed 149+199+596 vs 944, both strands;
  the same TR allele with breakpoints 208 bp apart; the same DEL placed at two
  copies of an imperfect repeat.

Kept apart (as they should be):
- B: chr1:54,510,617 (child +280 vs mother +330/+375), chr14:23,280,741.
- C: chr2:194,106,856 (child 95 with the father; mother 115 apart).
- chr11:62,386,440 (plan D3's 121 vs 104, label C): 17 bp between the
  haplotypes > 10.4 bp -> distinct. Q100 says the child is 104/121 and the
  child's single consensus collapsed to 121, so a merge would only hide that.

New errors:
- 298 new violations, mostly splits of records the old rule had made consistent
  by size (see above). Some are probably consensus errors of 10-20 bp in a TR
  (chr3:194,527,040), others real haplotype differences in VNTRs.
- Two consensuses of the child linked directly although their haplotypes differ
  within the tolerance (the transitive gate needs a direct link, so collapses can
  only come from the tolerance): records holding two child consensuses whose
  summed sizes differ > 10% went 405 -> 508; the ones inspected were one allele
  decomposed differently in the two consensuses (e.g. chr1:2,440,512: 1,372 vs
  865+166+..., 36 bp between the haplotypes), not two alleles.

## D5 probe: genotyping a parent from the same alignments (`scripts/d5_probe.py`, v4 code)

Child-carrying records (>= 50 bp) in which a parent has no pattern (today a
forced 0/0 unless its reads say otherwise), the parent's consensuses that carry
INS/DEL composites in the same candidate group, child's allele judged on them:

| parent | present | absent | not evaluable | no consensus in the group |
|---|---:|---:|---:|---:|
| HG003 (9,765 records) | 1,016 (10%) | 3,453 (35%) | 1,671 (17%) | 3,625 (37%) |
| HG004 (9,279 records) | 964 (10%) | 3,347 (36%) | 1,699 (18%) | 3,269 (35%) |

"present" records are allele carriers the record layer missed (could be genotyped
alt from that consensus's reads); "absent" gives consensus-level evidence for 0/0.
The 35-37% without a consensus in the group need the consensuses that carry no
SV at the locus (reference-like haplotypes); they are in `consensus_sequences`
but are not placed on the reference in the DB except through `region=`.

## v6: same-consensus guard (follow-up to the D3 study)

D3 (exp/merge-repeat-aware 8189a30) found that most of its `--merge-cannot-link`
gain comes from pairs of the SAME consensus bridged into one record through
another sample's composite. v1-v5 skipped same-consensus pairs (so pieces of a
decomposed allele could join through a partner), which also let partners bridge
two distinct events of one consensus.

v6 (14881fe): two composites of one (sample, consensusID) get a cannot-link
unless the alignment links both of them (corresponding places, agreeing
haplotypes) to the same partner composite: that is the signature of one allele
carried in pieces on this consensus and in one piece on the other (the partner's
range projects over both pieces and the haplotypes agree across the hull).
Links of the pairwise rule, twins, or links to two different partners do not
qualify. Rule links are applied after the alignment's, closest in size first (as
D3's `_cannot_link_components`). Inversions (not judged by the alignment) go
through the same guarded union. Genome-wide: 201,932 same-consensus pairs
guarded, 2,352 allowed as pieces, 17,571 unions blocked (v5: 7,756).

| arm | MC | P | R | F1 | H | strict F1 (pctsize/pctseq 0.9: P / R / F1) |
|---|---:|---:|---:|---:|---:|---|
| base_c8804b1 | 0.9166 | 0.8755 | 0.7649 | 0.8165 | 0.8473 | 0.8145 / 0.7095 / 0.7584 |
| d3_guard (8189a30, `--merge-cannot-link`) | 0.9135 | 0.8918 | 0.7834 | 0.8341 | 0.8590 | 0.8330 / 0.7309 / 0.7786 |
| proto_align_v5 | 0.9336 | 0.8861 | 0.7719 | 0.8251 | 0.8583 | 0.8256 / 0.7165 / 0.7672 |
| proto_align_v6 | 0.9310 | 0.8893 | 0.7755 | 0.8285 | 0.8600 | 0.8298 / 0.7209 / 0.7715 |

(T2T-Q100; strict numbers of base and d3_guard from d3/final_table_strict.tsv,
base re-run here with the same result.) Labels (class a, consistent / newly
joined to the transmitting parent): A1 395 / 375 (v5 410 / 386; d3_guard 6 / 0),
B 1 / 1 (v5 0 / 2; d3_guard 1 / 0), C 15, D_cons 34, D_gt 22, wrong-parent-only
attachments 2. Record transitions vs base: 883 inconsistent -> consistent, 388
new violations (v5: 904 / 298). Merge step 238 s, wall 18:52, peak RSS 20.3 GB.
