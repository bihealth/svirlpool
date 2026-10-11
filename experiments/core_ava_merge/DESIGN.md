# Core-vs-core alignment as the cross-sample merge criterion

Prototype on branch `exp/core-ava-merge` (v1-v6), productionized on
`feat/core-alignment-merge`, on by default since it went into
`fix/all-fixes-integration` (`--no-core-alignment-merge`: the pairwise rule only).
Switch: `svirlpool sv-calling --core-alignment-merge` (v6 design: the alignment
decides every pair it can judge, the pairwise rule the others; the prototype's
`or` mode is gone).
Code: `src/svirlpool/svcalling/core_alignment_merge.py`, hooked into
`svcomposite_merging.merge_svComposites` (`merge_svComposites_with_core_alignment`
schedules the candidate groups, `merge_group_with_core_alignment` merges one).
Tests: `tests/test_core_alignment_merge.py`. Results: `RESULTS.md` next to this
file.

## Question

Can "same allele" across samples be decided by comparing the consensus haplotypes
themselves (where the SVs sit on them), instead of SV size (10%), a 9-mer test and
reference proximity / shared repeatIDs? Everything must come from the per-sample
DBs (no BAMs).

## Inputs (all from the DBs)

- `consensus_sequences`: padded consensus + `core=S-E`. Only the core is used
  (the padding comes from single raw reads). Streaming the table and slicing the
  cores of the consensuses that carry an INS/DEL composite: 136k cores, 262 Mb
  for the trio, 4-6 s per sample from page cache.
- `svPatterns`: `read_start/read_end` of the primitives are forward positions on
  the padded consensus (verified: for all 122k INS patterns the stored inserted
  sequence equals the core slice, reverse-complemented when `aln_is_reverse`;
  DELs are points, `read_end == read_start` in 99.8%).
- `aln_is_reverse` gives the strand of the consensus-to-reference alignment;
  `consensus_aln_interval` gives the core's reference placement (its start
  matches core position 0 to within a few bp for the median composite).

## Which pairs are aligned

The existing candidate groups are kept (composites whose core reference
intervals overlap, short ones padded by `--near`). Inside a group, every
consensus carrying an INS or DEL composite is aligned against every other one:
different samples, the two haplotype consensuses of one sample, and overlapping
CRs of one sample. Nothing outside the group (no paralog search).

Each core is first turned to the reference strand (majority `aln_is_reverse` of
its composites), so all pairs align on `+`.

## How

minimap2 in-process (mappy, preset asm20): each core of the group is indexed once
and the later cores are mapped onto it. Primary and supplementary alignments
are kept (the flanks of a large indel come back as separate blocks), secondary
ones are not. mappy does not report the DP score, so each block is scored with
asm20's scoring (-A1 -B4 -O6,26 -E2,1) from its cs string. Groups with more than
300 kb of cores are not aligned (their pairs fall back to the pairwise rule).
Nothing depends on time, so the result is deterministic.

History: the prototype ran the minimap2 binary, one process per target core
(about 96k per trio), with `-U c,c`, `c = max(10, 3 * cores in the group)`, and a
10 s wall-clock cap per group, which made the set of aligned groups depend on
machine load. mappy cannot set the seed-occurrence cap (it fixes `mid_occ` at
1000 for an index built from a sequence); on chr21 this changed nothing that
matters: 2,181 of 2,195 records identical, the satellite groups no slower,
alignment time 32.3 s -> 23.8 s. The alignment time is concentrated in a few
large groups (chr21: the top 10 of 677 groups take 69%), so the groups of all
chromosomes are merged in parallel group by group, largest first, and
reassembled in chromosome and group order.

Why not one all-vs-all call: `-X` all-vs-all reports every secondary chain; in
tandem repeats that is thousands of chains per pair (one 6-core, 44 kb group:
19,143 PAF lines, 22 s).

Per pair, the blocks inside a band around the diagonal the two reference
placements imply (`|(t - q) - (ref_q - ref_t)| <= 1000 + the summed SV sizes of
both cores`) are chained (best co-linear, non-overlapping chain). The band drops
blocks to another copy of a duplicated sequence inside a core (seen: an exact
same-position DEL projected 11 kb away). Gaps between blocks count as
unaligned on both sides.

## Decision for a pair of composites x (on core A) and y (on core B)

Same SV type only (INS with INS, DEL with DEL), different consensuses.

1. Placement range of x on A: its interval (the inserted bases; a point for a
   DEL), extended over the bases across which the same event can slide (exact
   shift of the inserted bases / of the deleted sequence along the core, at most
   1 kb each side). Inside a perfect tandem repeat this covers the repeat.
2. Projection: x's range through the A-B chain onto B (P_x), y's range through
   the inverted chain onto A (P_y). Correspondence: P_x overlaps y's range or
   P_y overlaps x's range (+- 50 bp).
3. Haplotype agreement: the indel bases of the A-B alignment inside
   hull(x's range, P_y) +- 150 bp flank (clipped to the aligned part) are
   <= max(10 bp, 0.1 * smaller SV size), and x's inserted bases carry <= 20% (+2)
   mismatches. A stretch unaligned on both sides (between blocks) counts
   max(gap on A, gap on B); insertions/deletions less than 50 bp apart count
   their net length (v5: inside a repeat minimap2 writes one difference as
   alternating 2-bp insertions and deletions, chr7:87626349: 36 indel bases,
   net 11). At least 30 aligned bp are needed on each side of the hull.
   The hull matters: a first version judged a window around x alone ("is x
   present on B?"); in long, imperfect VNTRs that window looked identical while
   the other haplotype carried 1.1 kb of extra units next to it, and the
   child's two alleles at chr10:103251035 (MAT +1062, PAT +2198) ended up in
   one record.
4. **link** = correspondence and agreement; **distinct** = correspondence and
   no agreement.
5. Both ranges must project. When one does not although the reference spans of
   both cores cover both SVs, the alignment breaks inside the locus (a large
   difference or an unalignable repeat): judged, not linked. When the cores do
   not both cover the SVs, or there is no alignment, or the SV has no aligned
   flanks, or the composite spans two consensuses/strands, the pair is not
   evaluable and the existing pairwise rule decides. (The prototype's `or`
   mode, which kept every link of the existing rule and added the alignment's,
   lost F1 and is gone.)
6. Twins: no correspondence, but the haplotypes agree across both SVs within
   10 bp (no relative tolerance), and neither SV is linked to another composite
   of the other consensus: one haplotype difference placed differently by the
   two consensus-to-reference alignments (an imperfect repeat; the 40-bp DEL
   placed 169 bp apart at chr22:30631699/30631908). Linked. (v2 handed these to
   the pairwise rule, which rejects them on size; 35 A1 records stayed unjoined
   for that reason.)

### Tolerance ("near-identical")

Indel bases between the two *haplotypes* around the SVs, not the difference of
the reference-alignment indel sizes. `max(10, 0.1 * size of the smaller SV)`
(v1-v3: the larger; a 973-bp DEL then lent 97 bp of tolerance to a 47-bp
piece 57 bp different):
- 10 bp absorb consensus homopolymer/STR slippage (ONT R10; the tr_substring
  study finds median consensus identity 1.000 to Q100, read sigma 1.5-7 bp);
- 0.1 * size mirrors the existing 10% size rule, but on sequence. A tandem-repeat
  length difference of 17 bp in a 104/121 bp allele (14%) stays a different
  allele; `--core-alignment-frac-tol 0.2` merges it. Different decompositions of
  one haplotype have ~0 indel bases between the haplotypes and link at any
  tolerance; that is the case the user's proposal targets (tr_substring A1: the
  same Q100 haplotype in child and transmitting-parent consensus in 68% of the
  size-discordant records). B (parent reads lack the child's allele, typically
  several STR units, 10-55 units) should stay apart at 0.1.
- Measured (RESULTS.md): 0.2 instead of 0.1 gives MC +0.09 pt, the same F1/H,
  11 more A1 records consistent, and 3 of the 21 B records consistent (0 at 0.1).
  0.1 stays the default.

### Different decompositions of one haplotype

Composites of the same consensus are never linked directly (they are distinct
events on one haplotype), but each piece links to the partner's single event
(INS 149, 199, 596 on the child each land inside the mother's INS 944), so
union-find puts all of them into one record. Representative unchanged: the
best-supported (sample, consensus) group; its pieces are summed (size) and
concatenated (ALT), its strongest pattern gives the position.

### Constraint (transitive haplotype gate)

Links are applied best first (alignment links by indel/size, twins after them,
then the rule's). Two composites of one sample from different consensuses (two
haplotypes, or two overlapping CRs) may share a record only if they are linked
directly; a link that would join them through partners is skipped. This keeps
a compound-heterozygous child's two alleles apart. (v1-v3 blocked only pairs
judged `distinct`; the harness then showed 405 -> 657 records holding two child
consensuses whose sizes differ by more than 10%; most were one allele
decomposed differently, some were chains.)

Two composites of one consensus (v6) share a record only as pieces of one allele
that another consensus carries in one piece: both linked by the alignment to the
same partner. Any other path (fallback links, twins, links to different
partners) must not bridge them. Inversions, which the alignment does not judge,
get the pairwise rule with the same transitive gate (two composites of one
sample whose pair the rule rejects never share a record; accepted pairs are
applied closest in size first).

### Samples without a consensus at the locus

Unchanged (their genotype is the D5 problem). The same alignments could genotype
them: a consensus of S on which x is *present* (step 3) is evidence for x in S,
one on which x is evaluable and absent is evidence against; see RESULTS.md.
