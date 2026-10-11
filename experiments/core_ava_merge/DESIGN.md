# Core-vs-core alignment as the cross-sample merge criterion (prototype)

Branch `exp/core-ava-merge` (prototype; default behaviour unchanged).
Switch: `svirlpool sv-calling --core-alignment-merge {off,align,or}`.
Code: `src/svirlpool/svcalling/core_alignment_merge.py`, hooked into
`svcomposite_merging.merge_svComposites` (per-chromosome worker
`merge_svComposites_for_chromosome_core_aln`). Tests:
`tests/test_core_alignment_merge.py`. Results: `RESULTS.md` next to this file.

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

minimap2 2.31 (the env's binary), `-x asm20 -c --eqx --secondary=no`, one call
per target core with the remaining cores of the group as queries, `-U c,c` with
`c = max(10, 3 * cores in the group)` (minimizers more frequent than that are
not anchors), at most 10 s per group (then: no alignment, fallback).

Why this and not one all-vs-all call: `-X` all-vs-all reports every secondary
chain; in tandem repeats that is thousands of chains per pair (one 6-core, 44 kb
group: 19,143 PAF lines, 22 s). Per-target `--secondary=no` returns the primary
alignment plus supplementary blocks (the two flanks of a large indel) per pair:
6x faster on the slow tail. No edlib/mappy/parasail in the env.

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
   <= max(10 bp, 0.1 * larger SV size), and x's inserted bases carry <= 20% (+2)
   mismatches. Gaps between blocks count max(gap on A, gap on B). At least
   30 aligned bp are needed on each side of the hull.
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
   evaluable and the existing pairwise rule decides (`align` mode). In `or`
   mode the existing rule's links are always kept and alignment links added.
6. One case is handed back to the existing rule: no correspondence, but the
   haplotypes agree across both SVs, and neither SV is linked to another
   composite of the other consensus. That is one event placed differently by
   the two consensus-to-reference alignments (an imperfect repeat; the 40-bp
   DEL placed 169 bp apart at chr22:30631699/30631908).

### Tolerance ("near-identical")

Indel bases between the two *haplotypes* around the SV, not the difference of
the reference-alignment indel sizes. `max(10, 0.1 * size)`:
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

### Different decompositions of one haplotype

Composites of the same consensus are never linked directly (they are distinct
events on one haplotype), but each piece links to the partner's single event
(INS 149, 199, 596 on the child each land inside the mother's INS 944), so
union-find puts all of them into one record. Representative unchanged: the
best-supported (sample, consensus) group; its pieces are summed (size) and
concatenated (ALT), its strongest pattern gives the position.

### Constraint

Links are applied best first (alignment links by indel/size, then the rule's).
A link is skipped when it would put two composites of the same sample from
different consensuses into one record although the alignment found them to be
different alleles (`distinct`). This keeps a compound-heterozygous child's two
alleles apart even when partners would chain them.

### Samples without a consensus at the locus

Unchanged (their genotype is the D5 problem). The same alignments could genotype
them: a consensus of S on which x is *present* (step 3) is evidence for x in S,
one on which x is evaluable and absent is evidence against; see RESULTS.md.
