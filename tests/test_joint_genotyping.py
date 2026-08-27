"""Unit tests for joint genotyping in ``svirlpool.svcalling.multisample_sv_calling``.

These tests pin the behaviour targeted by defect **F9** ("absence of evidence is
emitted as a confident homozygous-reference call"):

* a locus with **no read coverage at all** must not produce a genotype — it must
  produce a no-call (``./.``) with ``TC = 0`` and ``GQ = 0``;
* a locus with **N reference reads and no alternate reads** may produce ``0/0``,
  but its ``GQ`` must be derived from the depth ``N`` (growing with ``N`` and
  reaching the ceiling of 60 only at high depth), not hard-coded to 60;
* ``genotype_likelihood`` with ``n_total_reads == 0`` must signal *no
  information* rather than certainty of ``0/0``;
* the pre-fix force-call must remain reproducible **exactly** behind the
  ``--legacy-force-wildtype-genotypes`` control flag.

Everything is built from genuine objects: real ``IntervalTree`` coverage trees
whose interval payloads are the same ``xxh64`` read-name hashes the pipeline
uses, and real ``SVprimitive`` / ``SVpattern`` / ``SVcomposite`` instances
carrying ``GenotypeMeasurement``s.
"""

import math

import pytest
from intervaltree import IntervalTree
from xxhash import xxh64

from svirlpool.localassembly import SVpatterns, SVprimitives
from svirlpool.svcalling import genotyping
from svirlpool.svcalling.multisample_sv_calling import (
    Genotype,
    SVcalls_from_SVcomposite,
    create_dummy_covtrees_from_reference,
    create_wild_type_genotype,
    generate_header,
    genotype_likelihood,
    genotype_of_sample,
)
from svirlpool.svcalling.SVcomposite import SVcomposite

SAMPLE = "sample1"
CHR = "chr1"
GQ_CEILING = 60

# ---------------------------------------------------------------------------
# Helper factories — real objects, no mocks
# ---------------------------------------------------------------------------


def _hash(readname: str) -> int:
    """The exact read-name hash the pipeline stores in covtrees / alt-read sets."""
    return int(xxh64(readname).intdigest())


def _covtree(
    reads: list[str],
    start: int = 900,
    end: int = 1600,
    chrname: str = CHR,
    samplename: str = SAMPLE,
) -> dict[str, dict[str, IntervalTree]]:
    """Build a real per-sample covtree in which every named read spans [start, end)."""
    tree = IntervalTree()
    for readname in reads:
        tree.addi(start, end, _hash(readname))
    return {samplename: {chrname: tree}}


def _covtree_empty(
    chrname: str = CHR, samplename: str = SAMPLE
) -> dict[str, dict[str, IntervalTree]]:
    """A covtree that knows the chromosome but holds no read intervals at the locus."""
    tree = IntervalTree()
    # A read far away from the loci used below, so the chromosome key exists but
    # a query at 1000-1500 returns nothing.
    tree.addi(10_000_000, 10_000_500, _hash("far_away_read"))
    return {samplename: {chrname: tree}}


def _genotype_of(
    reads_in_covtree: list[str],
    alt_reads: list[str],
    *,
    covtrees: dict[str, dict[str, IntervalTree]] | None = None,
    start: int = 1000,
    end: int = 1500,
    cn_tracks: dict | None = None,
    **kwargs,
) -> Genotype:
    if covtrees is None:
        covtrees = _covtree(reads_in_covtree)
    return genotype_of_sample(
        samplename=SAMPLE,
        chrname=CHR,
        start=start,
        end=end,
        raw_alt_reads={_hash(r) for r in alt_reads},
        covtrees=covtrees,
        cn_tracks=cn_tracks if cn_tracks is not None else {},
        **kwargs,
    )


def _reads(n: int, prefix: str = "read") -> list[str]:
    return [f"{prefix}{i}" for i in range(n)]


def _make_svprimitive_ins(
    *,
    chrname: str = CHR,
    ref_start: int = 1000,
    size: int = 500,
    samplename: str = SAMPLE,
    consensusID: str = "1.0",
    reads: list[str],
) -> SVprimitives.SVprimitive:
    return SVprimitives.SVprimitive(
        ref_start=ref_start,
        ref_end=ref_start + 1,
        read_start=0,
        read_end=size,
        size=size,
        sv_type=0,  # INS
        chr=chrname,
        repeatIDs=[],
        original_alt_sequences=["A" * size],
        original_ref_sequences=[],
        samplename=samplename,
        consensusID=consensusID,
        alignmentID=0,
        svID=0,
        aln_is_reverse=False,
        consensus_aln_interval=(chrname, ref_start - 500, ref_start + 500),
        genotypeMeasurement=genotyping.GenotypeMeasurement(
            start_on_consensus=0,
            supporting_reads_start=list(reads),
        ),
    )


def _make_insertion_composite(
    *,
    size: int = 500,
    ref_start: int = 1000,
    samplename: str = SAMPLE,
    consensusID: str = "1.0",
    reads: list[str],
) -> SVcomposite:
    svp = _make_svprimitive_ins(
        ref_start=ref_start,
        size=size,
        samplename=samplename,
        consensusID=consensusID,
        reads=reads,
    )
    pattern = SVpatterns.SVpatternInsertion(SVprimitives=[svp], size_distortions=None)
    pattern.set_sequence("A" * size)
    return SVcomposite.from_SVpattern(pattern)


# ===========================================================================
# 1. The catalogue's "Verification" bar
# ===========================================================================


class TestNoCoverageIsNotAReferenceCall:
    """A locus outside any read interval must yield ``./.`` with ``GQ = 0``."""

    def test_locus_outside_any_read_interval_is_a_no_call(self):
        gt = _genotype_of([], ["altread"], covtrees=_covtree_empty())
        assert gt.genotype == "./."
        assert gt.genotype_quality == 0
        assert gt.total_coverage == 0
        assert gt.ref_reads == 0
        assert gt.var_reads == 0

    def test_no_call_does_not_claim_a_genotype_probability(self):
        """``GP`` must not assert probability 1.0 for a genotype we did not call."""
        gt = _genotype_of([], ["altread"], covtrees=_covtree_empty())
        assert gt.gt_likelihood is None
        assert gt.to_format_field("GT:GQ:TC:DR:DV:GP") == "./.:0:0:0:0:."

    def test_chromosome_absent_from_covtree_is_a_no_call(self):
        """Call site 1: the chromosome is missing from the sample's coverage tree."""
        gt = _genotype_of(
            [],
            ["altread"],
            covtrees={SAMPLE: {"chrOTHER": IntervalTree()}},
        )
        assert gt.genotype == "./."
        assert gt.genotype_quality == 0
        assert gt.total_coverage == 0

    def test_sample_absent_from_covtrees_is_a_no_call(self):
        gt = _genotype_of([], ["altread"], covtrees={})
        assert gt.genotype == "./."
        assert gt.genotype_quality == 0


class TestReferenceCallQualityFollowsDepth:
    """N reference reads and no alt reads -> ``0/0`` with a depth-derived ``GQ``."""

    def test_reference_reads_with_no_alt_reads_give_wild_type(self):
        """Call site 2: coverage present, no alternate support."""
        gt = _genotype_of(_reads(10), ["unrelated_alt_read"])
        assert gt.genotype == "0/0"
        assert gt.total_coverage == 10
        assert gt.ref_reads == 10
        assert gt.var_reads == 0

    @pytest.mark.parametrize("depth", [1, 2, 3, 5, 10, 15, 20, 30, 60])
    def test_gq_never_exceeds_the_ceiling(self, depth):
        gt = _genotype_of(_reads(depth), ["unrelated_alt_read"])
        assert 0 <= gt.genotype_quality <= GQ_CEILING

    def test_gq_increases_with_depth(self):
        """Monotonicity: more reference support must never mean lower quality."""
        qualities = [
            _genotype_of(_reads(n), ["unrelated_alt_read"]).genotype_quality
            for n in range(1, 41)
        ]
        assert qualities == sorted(qualities), qualities
        # and it must actually move, not be a flat line
        assert qualities[0] < qualities[-1]

    def test_gq_is_low_at_low_depth(self):
        """A single reference read is not a confident homozygous-reference call."""
        assert _genotype_of(_reads(1), ["unrelated_alt_read"]).genotype_quality < 10

    def test_gq_reaches_the_ceiling_only_at_high_depth(self):
        assert _genotype_of(_reads(5), ["unrelated_alt_read"]).genotype_quality < (
            GQ_CEILING
        )
        assert (
            _genotype_of(_reads(40), ["unrelated_alt_read"]).genotype_quality
            == GQ_CEILING
        )

    def test_gp_is_a_real_posterior_not_a_hard_coded_one(self):
        low = _genotype_of(_reads(2), ["unrelated_alt_read"])
        high = _genotype_of(_reads(30), ["unrelated_alt_read"])
        assert low.gt_likelihood is not None and high.gt_likelihood is not None
        assert low.gt_likelihood < 1.0
        assert low.gt_likelihood < high.gt_likelihood


# ===========================================================================
# 2. The three documented call sites, exercised distinctly
# ===========================================================================


class TestCallSitesAreDistinguishable:
    def test_no_alt_reads_branch_is_the_one_reached_when_coverage_is_empty(
        self, caplog
    ):
        """The third branch (``len(all_reads) == 0``) is unreachable.

        ``alt_reads = all_reads & raw_alt_reads``, so an empty ``all_reads``
        forces an empty ``alt_reads`` and the ``len(alt_reads) == 0`` branch
        returns first.  The two branches log different messages, which is what
        this test discriminates on.
        """
        with caplog.at_level("WARNING"):
            _genotype_of([], ["altread"], covtrees=_covtree_empty())
        messages = " ".join(r.getMessage() for r in caplog.records)
        assert "No alt reads found" in messages
        assert "No ref/alt reads found" not in messages

    def test_alt_reads_are_always_a_subset_of_all_reads(self):
        """The invariant that makes the third branch unreachable."""
        # An alt read that is *not* in the covtree cannot become an alt read.
        gt = _genotype_of(_reads(5), ["not_in_covtree"])
        assert gt.var_reads == 0
        assert gt.ref_reads == 5

    def test_alt_supported_call_is_unaffected(self):
        """Negative control: a well-supported call must not change at all."""
        alt = _reads(5, prefix="alt")
        ref = _reads(5, prefix="ref")
        gt = _genotype_of(alt + ref, alt)
        assert gt.genotype == "0/1"
        assert gt.var_reads == 5
        assert gt.ref_reads == 5
        assert gt.total_coverage == 10
        assert gt.genotype_quality > 0

    def test_all_reads_alt_gives_homozygous(self):
        alt = _reads(20, prefix="alt")
        gt = _genotype_of(alt, alt)
        assert gt.genotype == "1/1"
        assert gt.var_reads == 20
        assert gt.ref_reads == 0


# ===========================================================================
# 3. genotype_likelihood with no observations
# ===========================================================================


class TestGenotypeLikelihoodWithoutObservations:
    def test_zero_total_reads_is_not_certainty_of_reference(self):
        probs = genotype_likelihood(n_alt_reads=0, n_total_reads=0, cn=2)
        assert probs["0/0"] < 1.0
        assert max(probs.values()) < 1.0

    def test_zero_total_reads_is_uninformative(self):
        """No observations -> no genotype is preferred over any other."""
        probs = genotype_likelihood(n_alt_reads=0, n_total_reads=0, cn=2)
        assert set(probs) == {"0/0", "0/1", "1/1"}
        assert all(math.isclose(p, 1.0 / 3.0, abs_tol=1e-9) for p in probs.values())

    def test_zero_total_reads_is_uninformative_for_haploid(self):
        probs = genotype_likelihood(n_alt_reads=0, n_total_reads=0, cn=1)
        assert set(probs) == {"0", "1"}
        assert all(math.isclose(p, 0.5, abs_tol=1e-9) for p in probs.values())

    def test_probabilities_still_sum_to_one(self):
        probs = genotype_likelihood(n_alt_reads=0, n_total_reads=0, cn=2)
        assert math.isclose(sum(probs.values()), 1.0, abs_tol=1e-9)

    def test_homozygous_deletion_locus_is_still_certain(self):
        """CN=0 genuinely implies no copies; that certainty is not the defect."""
        assert genotype_likelihood(n_alt_reads=0, n_total_reads=0, cn=0) == {"0": 1.0}

    @pytest.mark.parametrize("depth", list(range(1, 41)))
    def test_reference_posterior_grows_with_depth(self, depth):
        probs = genotype_likelihood(n_alt_reads=0, n_total_reads=depth, cn=2)
        assert probs["0/0"] == max(probs.values())

    def test_reference_posterior_is_monotone_in_depth(self):
        posteriors = [
            genotype_likelihood(n_alt_reads=0, n_total_reads=n, cn=2)["0/0"]
            for n in range(0, 41)
        ]
        assert posteriors == sorted(posteriors), posteriors


# ===========================================================================
# 4. create_wild_type_genotype directly
# ===========================================================================


class TestCreateWildTypeGenotype:
    def test_zero_coverage_is_refused(self):
        gt = create_wild_type_genotype(samplename=SAMPLE, total_coverage=0)
        assert gt.genotype == "./."
        assert gt.genotype_quality == 0
        assert gt.total_coverage == 0

    def test_positive_coverage_gives_full_reference_support(self):
        gt = create_wild_type_genotype(samplename=SAMPLE, total_coverage=12)
        assert gt.genotype == "0/0"
        assert gt.ref_reads == 12
        assert gt.var_reads == 0
        assert gt.total_coverage == 12

    def test_haploid_locus_gets_a_haploid_wild_type(self):
        gt = create_wild_type_genotype(
            samplename=SAMPLE, total_coverage=12, copy_number=1
        )
        assert gt.genotype == "0"

    def test_quality_is_monotone_in_coverage(self):
        qualities = [
            create_wild_type_genotype(
                samplename=SAMPLE, total_coverage=n
            ).genotype_quality
            for n in range(1, 41)
        ]
        assert qualities == sorted(qualities)
        assert qualities[0] < qualities[-1]


# ===========================================================================
# 5. The --legacy-force-wildtype-genotypes control
# ===========================================================================


class TestLegacyControlReproducesThePreFixOutput:
    """The control must reproduce the pre-fix record field for field."""

    LEGACY_FIELDS = ("genotype", "gt_likelihood", "genotype_quality")

    def test_legacy_wild_type_helper_is_the_old_force_call(self):
        for coverage in (0, 1, 7, 30):
            gt = create_wild_type_genotype(
                samplename=SAMPLE, total_coverage=coverage, legacy_force_call=True
            )
            assert gt == Genotype(
                samplename=SAMPLE,
                genotype="0/0",
                gt_likelihood=1.0,
                genotype_quality=60,
                total_coverage=coverage,
                ref_reads=coverage,
                var_reads=0,
            )

    def test_legacy_flag_restores_force_call_for_uncovered_locus(self):
        gt = _genotype_of(
            [], ["altread"], covtrees=_covtree_empty(), legacy_force_wildtype=True
        )
        assert gt.genotype == "0/0"
        assert gt.genotype_quality == 60
        assert gt.gt_likelihood == 1.0
        assert gt.total_coverage == 0
        assert gt.ref_reads == 0
        assert gt.var_reads == 0

    def test_legacy_flag_restores_force_call_for_missing_chromosome(self):
        gt = _genotype_of(
            [],
            ["altread"],
            covtrees={SAMPLE: {"chrOTHER": IntervalTree()}},
            legacy_force_wildtype=True,
        )
        assert gt.genotype == "0/0"
        assert gt.genotype_quality == 60
        assert gt.gt_likelihood == 1.0
        assert gt.total_coverage == 0

    @pytest.mark.parametrize("depth", [1, 3, 10, 40])
    def test_legacy_flag_restores_gq_60_at_every_depth(self, depth):
        gt = _genotype_of(
            _reads(depth), ["unrelated_alt_read"], legacy_force_wildtype=True
        )
        assert gt.genotype == "0/0"
        assert gt.genotype_quality == 60
        assert gt.gt_likelihood == 1.0
        assert gt.total_coverage == depth
        assert gt.ref_reads == depth
        assert gt.var_reads == 0

    def test_legacy_flag_does_not_change_this_alt_supported_call(self):
        """DV=5/TC=10 scores 33, below the ceiling, so capping is a no-op here.

        Not a general property of alternate-supported calls: N16 gave the
        control the pre-fix *uncapped* quality on that path too, so the two
        agree only where the pre-fix value was already <= GQ_CEILING.  See
        ``TestLegacyControlReproducesTheUncappedQuality`` for the cases where
        they diverge.
        """
        alt = _reads(5, prefix="alt")
        ref = _reads(5, prefix="ref")
        fixed = _genotype_of(alt + ref, alt)
        legacy = _genotype_of(alt + ref, alt, legacy_force_wildtype=True)
        assert fixed.genotype_quality == 33
        assert fixed == legacy

    def test_legacy_genotype_likelihood_returns_the_old_certainty(self):
        assert genotype_likelihood(
            n_alt_reads=0, n_total_reads=0, cn=2, legacy_zero_coverage=True
        ) == {"0/0": 1.0, "0/1": 0.0, "1/1": 0.0}


# ===========================================================================
# 6. The --skip-covtrees / dummy coverage path
# ===========================================================================


@pytest.fixture
def dummy_covtrees(tmp_path):
    """Real ``create_dummy_covtrees_from_reference`` output, from a real .fai."""
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text(f"{CHR}\t248956422\t6\t60\t61\n")
    return create_dummy_covtrees_from_reference(
        reference=ref, samplenames=[SAMPLE], default_coverage=30
    )


class TestDummyCoverageCannotProduceAConfidentCall:
    """``--skip-covtrees`` is an investigation tool and gets no special-casing.

    The catalogue worried that the fabricated tracks would silently reinstate the
    defect, because the coverage query always reports "covered".  Deriving GQ
    from depth makes the flag self-limiting instead: the dummy tree holds a
    single interval per chromosome whose payload is the coverage *value*, so the
    read-hash query yields exactly one fabricated read genome-wide and the
    resulting quality is the lowest the model can produce, not the ceiling.
    """

    def test_dummy_covtrees_report_coverage_everywhere(self, dummy_covtrees):
        """Sanity check: the dummy track claims to cover the locus."""
        assert len(dummy_covtrees[SAMPLE][CHR][1000:1500]) == 1

    def test_dummy_coverage_yields_a_single_fabricated_read(self, dummy_covtrees):
        """One interval per chromosome, payload 30 -> one "read", at every locus."""
        for start, end in [(1000, 1500), (5_000_000, 5_000_100), (99, 100)]:
            gt = _genotype_of(
                [], ["altread"], covtrees=dummy_covtrees, start=start, end=end
            )
            assert gt.total_coverage == 1

    def test_dummy_coverage_does_not_produce_a_confident_wild_type(
        self, dummy_covtrees
    ):
        """The depth-derived quality collapses to the floor instead of the ceiling."""
        gt = _genotype_of([], ["altread"], covtrees=dummy_covtrees)
        assert gt.genotype == "0/0"
        assert gt.genotype_quality < 10
        assert gt.genotype_quality != GQ_CEILING

    def test_legacy_flag_restores_the_ceiling_even_on_dummy_coverage(
        self, dummy_covtrees
    ):
        """The pre-fix behaviour under the flag: GQ=60 from a fabricated read."""
        gt = _genotype_of(
            [], ["altread"], covtrees=dummy_covtrees, legacy_force_wildtype=True
        )
        assert gt.genotype == "0/0"
        assert gt.genotype_quality == GQ_CEILING


# ===========================================================================
# 7. The VCF-level path: samples with no Genotype object at all
# ===========================================================================


def _svcalls(composite, covtrees, **kwargs):
    return SVcalls_from_SVcomposite(
        svComposite=composite,
        covtrees=covtrees,
        cn_tracks={},
        find_leftmost_reference_position=False,
        symbolic_threshold=100_000,
        **kwargs,
    )


class TestJointCallAcrossSamples:
    """A locus seen in one sample and absent from another — the joint-call case."""

    def test_uncovered_second_sample_is_a_no_call_in_the_vcf(self):
        alt = _reads(6, prefix="alt")
        composite = _make_insertion_composite(reads=alt)
        covtrees = _covtree(alt)
        # sample2 knows chr1 but has no reads at the locus
        covtrees["sample2"] = _covtree_empty(samplename="sample2")["sample2"]

        calls = _svcalls(composite, covtrees)
        assert len(calls) == 1
        line = calls[0].to_vcf_line(
            vcfIDnumber=0,
            samplenames=[SAMPLE, "sample2"],
            covtrees=covtrees,
            refdict={f"{CHR}:1000": "A"},
            symbolic_threshold=100_000,
        )
        assert line is not None
        fields = line.split("\t")
        gt_sample1, gt_sample2 = fields[-2], fields[-1]
        assert gt_sample1.startswith("0/1:") or gt_sample1.startswith("1/1:")
        assert gt_sample2 == "./.:0:0:0:0:."

    def test_covered_second_sample_is_a_justified_wild_type_in_the_vcf(self):
        alt = _reads(6, prefix="alt")
        composite = _make_insertion_composite(reads=alt)
        covtrees = _covtree(alt)
        covtrees["sample2"] = _covtree(
            _reads(11, prefix="s2read"), samplename="sample2"
        )["sample2"]

        calls = _svcalls(composite, covtrees)
        line = calls[0].to_vcf_line(
            vcfIDnumber=0,
            samplenames=[SAMPLE, "sample2"],
            covtrees=covtrees,
            refdict={f"{CHR}:1000": "A"},
            symbolic_threshold=100_000,
        )
        assert line is not None
        gt_sample2 = line.split("\t")[-1]
        gt, gq, tc, dr, dv, _gp = gt_sample2.split(":")
        assert gt == "0/0"
        assert int(tc) == 11
        assert int(dr) == 11
        assert int(dv) == 0
        assert 0 < int(gq) < GQ_CEILING

    def test_depth_counts_reads_not_alignment_fragments(self):
        """A read contributing two RAF intervals at the locus is one read.

        ``get_ref_reads_from_covtrees`` deduplicates by read hash everywhere
        else; the VCF-level fallback used to count intervals, which inflates the
        depth of long reads spanning large events (78 reads reported as TC=210
        at the 18 kb chr6 insertion).
        """
        alt = _reads(6, prefix="alt")
        composite = _make_insertion_composite(reads=alt)
        covtrees = _covtree(alt)
        tree = IntervalTree()
        tree.addi(900, 1600, _hash("s2read0"))
        tree.addi(950, 1200, _hash("s2read0"))  # second fragment of the same read
        tree.addi(900, 1600, _hash("s2read1"))
        covtrees["sample2"] = {CHR: tree}

        calls = _svcalls(composite, covtrees)
        line = calls[0].to_vcf_line(
            vcfIDnumber=0,
            samplenames=[SAMPLE, "sample2"],
            covtrees=covtrees,
            refdict={f"{CHR}:1000": "A"},
            symbolic_threshold=100_000,
        )
        assert line is not None
        gt, _gq, tc, dr, dv, _gp = line.split("\t")[-1].split(":")
        assert (gt, int(tc), int(dr), int(dv)) == ("0/0", 2, 2, 0)

        # ... and the control still reproduces the pre-fix (interval-counted) depth
        legacy_line = calls[0].to_vcf_line(
            vcfIDnumber=0,
            samplenames=[SAMPLE, "sample2"],
            covtrees=covtrees,
            refdict={f"{CHR}:1000": "A"},
            symbolic_threshold=100_000,
            legacy_force_wildtype=True,
        )
        assert legacy_line is not None
        assert legacy_line.split("\t")[-1] == "0/0:60:3:3:0:1.0"

    def test_legacy_flag_restores_the_force_call_in_the_vcf(self):
        alt = _reads(6, prefix="alt")
        composite = _make_insertion_composite(reads=alt)
        covtrees = _covtree(alt)
        covtrees["sample2"] = _covtree_empty(samplename="sample2")["sample2"]

        calls = _svcalls(composite, covtrees, legacy_force_wildtype=True)
        line = calls[0].to_vcf_line(
            vcfIDnumber=0,
            samplenames=[SAMPLE, "sample2"],
            covtrees=covtrees,
            refdict={f"{CHR}:1000": "A"},
            symbolic_threshold=100_000,
            legacy_force_wildtype=True,
        )
        assert line is not None
        assert line.split("\t")[-1] == "0/0:60:0:0:0:1.0"


# ===========================================================================
# 8. N15 — the breakpoint-placement margin is deletion-only
# ===========================================================================
#
# A read that carries an insertion at position ``p`` has, by construction,
# sequence at ``p`` that does not align to the reference at ``p``.  When the
# insertion is large or its flanks are repetitive, the aligner emits two
# fragments — one ending before ``p``, one resuming after it — and *no* fragment
# covering ``p``.  The very property that makes the read alt-supporting is what
# removes it from a point query at ``p``.
#
# The two breakpoints are also placed independently: once by the read-to-
# reference alignment (which produced the coverage track) and once by the
# consensus-to-reference alignment (which produced the call position).  In
# repetitive sequence those two placements disagree by a repeat-period-scale
# amount.  Measured on the muc1 example and eight HG002 tiles, that
# disagreement is 26–158 bp for 28 of 42 lost reads.
#
# ``breakpoint_mode`` recognises exactly this geometry — but is wired to
# deletions only (``is_deletion`` in ``svcall_object_from_svcomposite``).  These
# tests pin the asymmetry and its repair.

# The margin defaults to the same value as ``--near`` (150), the distance within
# which the merger already declares two independently placed breakpoints to be
# the same event.
DEFAULT_BREAKPOINT_MARGIN = 150


def _covtree_from_intervals(
    intervals: dict[str, list[tuple[int, int]]],
    chrname: str = CHR,
    samplename: str = SAMPLE,
) -> dict[str, dict[str, IntervalTree]]:
    """Build a covtree with explicit per-read alignment-fragment intervals."""
    tree = IntervalTree()
    for readname, spans in intervals.items():
        for begin, end in spans:
            tree.addi(begin, end, _hash(readname))
    return {samplename: {chrname: tree}}


def _fragmented_locus(
    *,
    gap: int,
    n_supporting: int = 3,
    n_spanning: int = 14,
    locus: int = 1000,
) -> tuple[dict[str, dict[str, IntervalTree]], list[str]]:
    """The muc1 1:69668 configuration, in miniature.

    ``n_supporting`` alt reads whose alignment fragment stops ``gap`` bases short
    of ``locus`` and resumes far downstream, plus ``n_spanning`` unrelated reads
    that cover ``locus`` contiguously.  The locus therefore looks well covered
    and quiet, which is precisely why the defect is invisible.
    """
    supporting = _reads(n_supporting, prefix="alt")
    intervals: dict[str, list[tuple[int, int]]] = {
        # ends at locus - gap + 1 (exclusive end), i.e. a gap of exactly `gap`
        name: [(locus - 400, locus - gap + 1), (locus + 10_000, locus + 10_300)]
        for name in supporting
    }
    for name in _reads(n_spanning, prefix="ref"):
        intervals[name] = [(locus - 400, locus + 400)]
    return _covtree_from_intervals(intervals), supporting


def _genotype_insertion(supporting, covtrees, **kwargs):
    """``genotype_of_sample`` with the geometry an insertion locus dispatches.

    ``svcall_object_from_svcomposite`` selects the query geometry from the SV
    type: ``breakpoint_mode`` for deletions, ``apply_breakpoint_margin`` for
    insertions.  These tests exercise the insertion branch directly; the
    composite-level tests below check the dispatch itself.
    """
    return _genotype_of(
        [],
        supporting,
        covtrees=covtrees,
        start=1000,
        end=1001,
        apply_breakpoint_margin=True,
        **kwargs,
    )


class TestInsertionBreakpointMargin:
    """Supporting reads whose alignment stops short of the insertion breakpoint."""

    def test_a_point_query_loses_every_supporting_read(self):
        """Characterisation: with no margin the sample is 0/0 for its own variant."""
        covtrees, supporting = _fragmented_locus(gap=118)
        gt = _genotype_insertion(supporting, covtrees, breakpoint_margin=0)
        assert gt.genotype == "0/0"
        assert gt.var_reads == 0
        assert gt.total_coverage == 14

    def test_the_margin_recovers_them_by_default(self):
        """The repair: the default margin counts the reads that produced the call."""
        covtrees, supporting = _fragmented_locus(gap=118)
        gt = _genotype_insertion(supporting, covtrees)
        assert gt.var_reads == 3
        assert gt.total_coverage == 17
        assert gt.ref_reads == 14

    def test_a_deletion_in_the_same_configuration_was_already_handled(self):
        """The asymmetry itself: ``breakpoint_mode`` recovers what the default loses.

        Same covtree, same reads, same locus — only the query geometry differs.
        """
        covtrees, supporting = _fragmented_locus(gap=118)
        as_insertion = _genotype_of(
            [],
            supporting,
            covtrees=covtrees,
            start=1000,
            end=1001,
        )
        as_deletion = _genotype_of(
            [],
            supporting,
            covtrees=covtrees,
            start=1000,
            end=1001,
            breakpoint_mode=True,
            breakpoint_margin=DEFAULT_BREAKPOINT_MARGIN,
        )
        assert as_insertion.var_reads == 0
        assert as_deletion.var_reads == 3

    @pytest.mark.parametrize(
        ("gap", "expected_alt"),
        [
            (1, 3),
            (50, 3),
            (DEFAULT_BREAKPOINT_MARGIN - 1, 3),
            (DEFAULT_BREAKPOINT_MARGIN, 3),  # just inside
            (DEFAULT_BREAKPOINT_MARGIN + 1, 0),  # just outside
            (DEFAULT_BREAKPOINT_MARGIN + 100, 0),
        ],
    )
    def test_margin_boundary_is_exact(self, gap, expected_alt):
        covtrees, supporting = _fragmented_locus(gap=gap)
        gt = _genotype_insertion(supporting, covtrees)
        assert gt.var_reads == expected_alt

    def test_margin_is_configurable(self):
        covtrees, supporting = _fragmented_locus(gap=400)
        assert _genotype_insertion(supporting, covtrees).var_reads == 0
        assert (
            _genotype_insertion(supporting, covtrees, breakpoint_margin=400).var_reads
            == 3
        )

    def test_composite_level_insertion_gets_the_margin(self):
        """End to end: an INS composite must not be 0/0 for its own supporting reads.

        This is the muc1 ``1:69668`` record: a PRECISE insertion assembled from
        six reads, genotyped ``0/0`` with ``DV=0``.
        """
        covtrees, supporting = _fragmented_locus(gap=118)
        composite = _make_insertion_composite(reads=supporting, ref_start=1000)
        calls = _svcalls(composite, covtrees)
        assert len(calls) == 1
        gt = calls[0].genotypes[SAMPLE]
        assert gt.var_reads == 3, "the insertion's own supporting reads are invisible"
        assert calls[0].pass_altreads is True

    def test_composite_level_deletion_dispatch_is_unchanged(self):
        """Deletions keep the two-breakpoint geometry; only the margin value moved."""
        covtrees, supporting = _fragmented_locus(gap=118)
        gt_bp = _genotype_of(
            [],
            supporting,
            covtrees=covtrees,
            start=1000,
            end=1001,
            breakpoint_mode=True,
        )
        assert gt_bp.var_reads == 3


class TestMarginDoesNotManufactureSupport:
    """The false-positive containment — the real danger of widening a window."""

    def test_a_nearby_unrelated_read_is_never_counted_as_alt(self):
        """A read near the breakpoint that the assembly did not name is not support.

        ``alt_reads = all_reads & raw_alt_reads``, so widening the *coverage*
        query can only add reads to ``TC``/``DR``.  A read can enter ``DV`` only
        if the assembly already listed it as supporting this composite.
        """
        covtrees, supporting = _fragmented_locus(gap=118)
        # A read whose alignment sits 120 bp before the breakpoint but which is
        # not among the composite's supporting reads.
        covtrees[SAMPLE][CHR].addi(600, 881, _hash("bystander"))
        gt = _genotype_insertion(supporting, covtrees)
        assert gt.var_reads == 3
        assert _hash("bystander") not in {_hash(r) for r in supporting}
        # the bystander is counted as reference depth, never as alt support
        assert gt.total_coverage == 18
        assert gt.ref_reads == 15

    @pytest.mark.parametrize("margin", [0, 1, 50, 150, 500, 5000])
    def test_dv_never_exceeds_tc(self, margin):
        """``DV > TC`` would mean the window double-counts; it cannot, by construction."""
        covtrees, supporting = _fragmented_locus(gap=118)
        gt = _genotype_insertion(supporting, covtrees, breakpoint_margin=margin)
        assert gt.var_reads <= gt.total_coverage
        assert gt.var_reads + gt.ref_reads == gt.total_coverage

    @pytest.mark.parametrize("margin", [0, 1, 50, 150, 500, 5000])
    def test_widening_never_loses_a_read_a_point_query_found(self, margin):
        """Monotonicity: the widened query is a strict superset of the point query."""
        covtrees, supporting = _fragmented_locus(gap=118)
        # give one supporting read an extra fragment that *does* span the locus
        covtrees[SAMPLE][CHR].addi(990, 1010, _hash("alt0"))
        gt = _genotype_insertion(supporting, covtrees, breakpoint_margin=margin)
        assert gt.var_reads >= 1

    def test_a_read_with_no_interval_anywhere_stays_uncounted(self):
        """Mechanism 2 of N15: reads absent from the coverage track altogether.

        ``alignments_to_rafs`` drops whole reads (the non-separated-read filter,
        ``min_mapq``, ``min_segment_size``) that the local assembler happily used
        as supporting reads.  Those reads have *no* interval on the chromosome,
        so no widening of the genotyping window can ever recover them.  This is a
        separate defect from the short-fragment one and is documented, not fixed,
        here.
        """
        covtrees, supporting = _fragmented_locus(gap=118)
        raw_alt = {_hash(r) for r in supporting} | {_hash("never_in_the_covtree")}
        gt = genotype_of_sample(
            samplename=SAMPLE,
            chrname=CHR,
            start=1000,
            end=1001,
            raw_alt_reads=raw_alt,
            covtrees=covtrees,
            cn_tracks={},
            apply_breakpoint_margin=True,
            breakpoint_margin=5000,
        )
        assert gt.var_reads == 3  # not 4 — the absent read is unrecoverable


class TestMarginLeavesTheF9PathsIntact:
    """The no-call and copy-number behaviour introduced by F9 must be unchanged."""

    def test_no_coverage_at_all_is_still_a_no_call(self):
        gt = _genotype_of(
            [], ["altread"], covtrees=_covtree_empty(), start=1000, end=1001
        )
        assert gt.genotype == "./."
        assert gt.genotype_quality == 0
        assert gt.total_coverage == 0

    def test_covered_locus_without_support_is_still_a_wild_type(self):
        gt = _genotype_of(_reads(10), ["unrelated_alt_read"])
        assert gt.genotype == "0/0"
        assert gt.var_reads == 0
        assert 0 < gt.genotype_quality <= GQ_CEILING

    @pytest.mark.parametrize(
        ("copy_number", "expected"),
        [(1, "1"), (2, "1/1"), (3, "1/1/1"), (4, "1/1/1/1")],
    )
    def test_non_diploid_copy_number_still_works_with_the_margin(
        self, copy_number, expected
    ):
        covtrees, supporting = _fragmented_locus(gap=118, n_spanning=0)
        cn_tree = IntervalTree()
        cn_tree.addi(0, 2_000_000, copy_number)
        gt = _genotype_insertion(
            supporting, covtrees, cn_tracks={SAMPLE: {CHR: cn_tree}}
        )
        assert gt.var_reads == 3
        assert gt.genotype == expected

    def test_haploid_locus_with_no_support_is_still_reference(self):
        cn_tree = IntervalTree()
        cn_tree.addi(0, 2_000_000, 1)
        gt = _genotype_of(
            _reads(10),
            ["unrelated_alt_read"],
            cn_tracks={SAMPLE: {CHR: cn_tree}},
        )
        assert gt.genotype == "0"
        assert gt.var_reads == 0

    def test_legacy_control_is_unaffected_by_the_margin(self):
        """``--legacy-force-wildtype-genotypes`` must still reproduce the old call."""
        gt = _genotype_of(
            [], ["altread"], covtrees=_covtree_empty(), legacy_force_wildtype=True
        )
        assert gt.genotype == "0/0"
        assert gt.genotype_quality == GQ_CEILING
        assert gt.total_coverage == 0


class TestCompositeLevelDispatch:
    """The behaviour a real run sees: geometry chosen from the SV type."""

    @pytest.mark.parametrize(
        ("gap", "expected_alt", "expected_filter"),
        [
            (118, 3, True),  # muc1 1:69668
            (DEFAULT_BREAKPOINT_MARGIN, 3, True),
            (DEFAULT_BREAKPOINT_MARGIN + 1, 0, False),
        ],
    )
    def test_insertion_composite_counts_short_reads_within_the_margin(
        self, gap, expected_alt, expected_filter
    ):
        covtrees, supporting = _fragmented_locus(gap=gap)
        composite = _make_insertion_composite(reads=supporting, ref_start=1000)
        call = _svcalls(composite, covtrees)[0]
        assert call.genotypes[SAMPLE].var_reads == expected_alt
        assert call.pass_altreads is expected_filter

    def test_insertion_composite_never_invents_support(self):
        """A bystander read inside the widened window stays reference depth."""
        covtrees, supporting = _fragmented_locus(gap=118)
        for i in range(5):
            covtrees[SAMPLE][CHR].addi(600, 881, _hash(f"bystander{i}"))
        composite = _make_insertion_composite(reads=supporting, ref_start=1000)
        gt = _svcalls(composite, covtrees)[0].genotypes[SAMPLE]
        assert gt.var_reads == 3
        assert gt.ref_reads == 19
        assert gt.total_coverage == 22

    def test_insertion_composite_margin_is_configurable_end_to_end(self):
        covtrees, supporting = _fragmented_locus(gap=400)
        composite = _make_insertion_composite(reads=supporting, ref_start=1000)
        assert _svcalls(composite, covtrees)[0].genotypes[SAMPLE].var_reads == 0
        assert (
            _svcalls(composite, covtrees, breakpoint_margin=500)[0]
            .genotypes[SAMPLE]
            .var_reads
            == 3
        )

    def test_a_well_supported_insertion_is_not_degraded(self):
        """Negative control, the shape of the chr6 18 kb INS (DR=2, DV=76).

        Widening the window may only add reference depth; it must not move a
        call that a point query already resolved.
        """
        alt = _reads(30, prefix="alt")
        covtrees = _covtree(alt + _reads(2, prefix="ref"))
        composite = _make_insertion_composite(reads=alt, ref_start=1000)
        gt = _svcalls(composite, covtrees)[0].genotypes[SAMPLE]
        assert gt.var_reads == 30
        assert gt.ref_reads == 2
        assert gt.genotype == "1/1"


class TestLegacyControlRestoresThePreFixWindow:
    """The campaign requires a setting that reproduces the pre-fix output exactly.

    The pre-fix geometry was *asymmetric* — 100 bases at a deletion's two
    breakpoints, none at all for an insertion — so no single value of
    ``--genotype-breakpoint-margin`` can express it.  The legacy control flag,
    which exists to reproduce pre-fix genotype fields, restores it.
    """

    def test_legacy_control_loses_the_short_reads_again(self):
        covtrees, supporting = _fragmented_locus(gap=118)
        composite = _make_insertion_composite(reads=supporting, ref_start=1000)
        assert _svcalls(composite, covtrees)[0].genotypes[SAMPLE].var_reads == 3
        legacy = _svcalls(composite, covtrees, legacy_force_wildtype=True)[0]
        assert legacy.genotypes[SAMPLE].var_reads == 0
        assert legacy.genotypes[SAMPLE].total_coverage == 14

    def test_legacy_control_overrides_an_explicit_margin(self):
        covtrees, supporting = _fragmented_locus(gap=118)
        composite = _make_insertion_composite(reads=supporting, ref_start=1000)
        legacy = _svcalls(
            composite,
            covtrees,
            breakpoint_margin=5000,
            legacy_force_wildtype=True,
        )[0]
        assert legacy.genotypes[SAMPLE].var_reads == 0


# ===========================================================================
# 9. N16 — the alternate-supported genotype quality
# ===========================================================================
#
# ``genotype_of_sample`` has two exits that carry a quality.  The
# homozygous-reference exit routes through ``phred_from_probability`` and is
# therefore bounded by ``GQ_CEILING`` and monotone in the evidence.  The
# alternate-supported exit phred-scaled the posterior inline, uncapped, with a
# literal ``60`` for the case where the posterior saturates to exactly ``1.0``
# in float64.  That produced two defects on one scale:
#
#   * GQ ran past the ceiling the VCF header declares (``GQ=154`` at
#     chr11:11,246,978 with DV=23/TC=50; ``GQ=63`` at muc1 1:248938 with
#     DV=2/TC=32);
#   * GQ was *non-monotone*: a moderately supported call scored 141 (DV=20,
#     TC=40) while a maximally supported one scored the literal 60 (DV=30,
#     TC=60).  Stronger evidence produced a lower quality, so any consumer
#     filtering or ranking on GQ got the ordering backwards at the top of the
#     range.
#
# These tests pin one scale — [0, GQ_CEILING], non-decreasing in support — for
# both exits, and pin the pre-fix expression behind the legacy control.


def _alt_supported(
    *, n_total: int, n_alt: int, prefix_alt: str = "alt", **kwargs
) -> Genotype:
    """A genotype from the alternate-supported exit: ``n_alt`` of ``n_total``.

    ``n_alt >= 1``, so ``alt_reads`` is non-empty and the wild-type exit is not
    taken — even where the argmax genotype is still ``0/0`` (the muc1 1:248938
    shape, DV=2 out of TC=32).
    """
    assert 1 <= n_alt <= n_total
    alt = _reads(n_alt, prefix=prefix_alt)
    ref = _reads(n_total - n_alt, prefix="ref")
    return _genotype_of(alt + ref, alt, **kwargs)


# (TC, DV, pre-fix GQ) — the uncapped expression's actual output, including the
# two field observations and the saturating cases that collapsed to the literal.
PREFIX_QUALITIES = [
    (10, 5, 33),
    (20, 10, 69),
    (32, 2, 63),  # muc1 1:248938 — argmax is 0/0, but DV=2 takes the alt exit
    (40, 20, 141),
    (50, 23, 154),  # chr11:11,246,978
    (50, 25, 60),  # posterior saturates -> the literal, below the 141 above
    (60, 30, 60),  # likewise
    (100, 50, 60),  # likewise
    (10, 10, 27),
    (20, 20, 55),
    (30, 30, 83),
    (50, 50, 139),
    (100, 100, 60),  # likewise
]


class TestAlternateCallQualityIsCapped:
    """The alternate-supported exit is bounded by the ceiling the header declares."""

    @pytest.mark.parametrize(
        ("n_total", "n_alt"), [(t, a) for t, a, _ in PREFIX_QUALITIES]
    )
    def test_gq_never_exceeds_the_ceiling(self, n_total, n_alt):
        gt = _alt_supported(n_total=n_total, n_alt=n_alt)
        assert gt.var_reads == n_alt
        assert 0 <= gt.genotype_quality <= GQ_CEILING

    @pytest.mark.parametrize("n_total", list(range(2, 61, 2)))
    def test_heterozygous_sweep_stays_within_the_ceiling(self, n_total):
        gt = _alt_supported(n_total=n_total, n_alt=n_total // 2)
        assert 0 <= gt.genotype_quality <= GQ_CEILING

    @pytest.mark.parametrize("n_total", list(range(1, 61)))
    def test_homozygous_sweep_stays_within_the_ceiling(self, n_total):
        gt = _alt_supported(n_total=n_total, n_alt=n_total)
        assert 0 <= gt.genotype_quality <= GQ_CEILING

    def test_the_two_observed_regressions_are_capped(self):
        """The two GQs seen on real data, both above the declared ceiling."""
        assert _alt_supported(n_total=50, n_alt=23).genotype_quality == GQ_CEILING
        muc1 = _alt_supported(n_total=32, n_alt=2)
        assert muc1.genotype == "0/0"  # the alt exit, despite the argmax
        assert muc1.var_reads == 2
        assert muc1.genotype_quality == GQ_CEILING

    def test_both_exits_share_one_scale(self):
        """Reference and alternate calls must be comparable, not two scales."""
        reference = _genotype_of(_reads(40), ["unrelated_alt_read"])
        alternate = _alt_supported(n_total=40, n_alt=20)
        assert reference.var_reads == 0 and alternate.var_reads == 20
        assert 0 <= reference.genotype_quality <= GQ_CEILING
        assert 0 <= alternate.genotype_quality <= GQ_CEILING


class TestAlternateCallQualityIsMonotone:
    """More alternate support must never mean a lower quality."""

    def test_saturating_call_does_not_score_below_a_moderate_one(self):
        """The inversion the literal ``60`` branch produced, in its sharpest form.

        DV=20/TC=40 and DV=30/TC=60 are the same allele fraction; the second has
        half again as much evidence.  Pre-fix the first scored 141 and the
        second 60.
        """
        moderate = _alt_supported(n_total=40, n_alt=20)
        strong = _alt_supported(n_total=60, n_alt=30)
        assert strong.genotype == moderate.genotype == "0/1"
        assert strong.genotype_quality >= moderate.genotype_quality

    def test_homozygous_saturating_call_does_not_score_below_a_moderate_one(self):
        """The same inversion on the 1/1 branch: DV=TC=50 scored 139, DV=TC=100 scored 60."""
        moderate = _alt_supported(n_total=50, n_alt=50)
        strong = _alt_supported(n_total=100, n_alt=100)
        assert strong.genotype == moderate.genotype == "1/1"
        assert strong.genotype_quality >= moderate.genotype_quality

    def test_gq_is_non_decreasing_across_a_heterozygous_depth_sweep(self):
        qualities = [
            _alt_supported(n_total=n, n_alt=n // 2).genotype_quality
            for n in range(2, 121, 2)
        ]
        assert qualities == sorted(qualities), qualities
        assert qualities[0] < qualities[-1]  # and it actually moves

    def test_gq_is_non_decreasing_across_a_homozygous_depth_sweep(self):
        qualities = [
            _alt_supported(n_total=n, n_alt=n).genotype_quality for n in range(1, 121)
        ]
        assert qualities == sorted(qualities), qualities
        assert qualities[0] < qualities[-1]

    def test_gq_reaches_the_ceiling_only_with_real_support(self):
        """The scale still discriminates: it is not a constant 60."""
        assert _alt_supported(n_total=4, n_alt=2).genotype_quality < GQ_CEILING
        assert _alt_supported(n_total=60, n_alt=30).genotype_quality == GQ_CEILING


class TestSingleEvidenceGenotypeQuality:
    """``--single-evidence-gt`` returns a fixed likelihood of exactly 1.0.

    Both the pre-fix expression (``else 60``) and ``phred_from_probability``
    (``probability >= 1.0 -> cap``) map that to 60, so this path is unchanged by
    the repair and needs no legacy special-casing.
    """

    @pytest.mark.parametrize(
        ("n_total", "n_alt"), [(10, 5), (10, 10), (50, 23), (60, 30)]
    )
    def test_single_evidence_quality_is_the_ceiling(self, n_total, n_alt):
        gt = _alt_supported(n_total=n_total, n_alt=n_alt, single_evidence_gt=True)
        assert gt.gt_likelihood == 1.0
        assert gt.genotype_quality == GQ_CEILING

    @pytest.mark.parametrize(("n_total", "n_alt"), [(10, 5), (10, 10), (50, 23)])
    def test_legacy_control_does_not_move_the_single_evidence_path(
        self, n_total, n_alt
    ):
        fixed = _alt_supported(n_total=n_total, n_alt=n_alt, single_evidence_gt=True)
        legacy = _alt_supported(
            n_total=n_total,
            n_alt=n_alt,
            single_evidence_gt=True,
            legacy_force_wildtype=True,
        )
        assert fixed == legacy


class TestLegacyControlReproducesTheUncappedQuality:
    """``--legacy-force-wildtype-genotypes`` must reproduce the pre-fix GQ exactly.

    The flag's contract is a controlled comparison against a pre-fix run, so it
    has to restore the uncapped expression on the alternate-supported exit too,
    not only the force-called wild types and the pre-fix query geometry.
    """

    @pytest.mark.parametrize(("n_total", "n_alt", "prefix_gq"), PREFIX_QUALITIES)
    def test_legacy_control_restores_the_pre_fix_quality(
        self, n_total, n_alt, prefix_gq
    ):
        legacy = _alt_supported(
            n_total=n_total, n_alt=n_alt, legacy_force_wildtype=True
        )
        assert legacy.genotype_quality == prefix_gq

    @pytest.mark.parametrize(("n_total", "n_alt", "prefix_gq"), PREFIX_QUALITIES)
    def test_only_the_quality_differs_between_control_and_repair(
        self, n_total, n_alt, prefix_gq
    ):
        """Field for field: GQ is the only thing the repair moves, and only down."""
        fixed = _alt_supported(n_total=n_total, n_alt=n_alt)
        legacy = _alt_supported(
            n_total=n_total, n_alt=n_alt, legacy_force_wildtype=True
        )
        for field in (
            "genotype",
            "gt_likelihood",
            "total_coverage",
            "ref_reads",
            "var_reads",
        ):
            assert getattr(fixed, field) == getattr(legacy, field), field
        assert fixed.genotype_quality == min(prefix_gq, GQ_CEILING)
        assert fixed.genotype_quality <= legacy.genotype_quality


class TestGenotypeQualityHeaderDescribesOneScale:
    """The header text has to describe the scale that is actually emitted."""

    def _gq_header(self, tmp_path) -> str:
        (ref := tmp_path / "ref.fa").write_text(">chr1\nACGT\n")
        (tmp_path / "ref.fa.fai").write_text(f"{CHR}\t248956422\t6\t60\t61\n")
        lines = [
            line
            for line in generate_header(reference=ref, samplenames=[SAMPLE])
            if line.startswith("##FORMAT=<ID=GQ,")
        ]
        assert len(lines) == 1
        return lines[0]

    def test_header_does_not_restrict_the_cap_to_reference_calls(self, tmp_path):
        header = self._gq_header(tmp_path)
        assert "homozygous-reference calls are capped" not in header

    def test_header_states_the_ceiling_and_the_no_call_value(self, tmp_path):
        header = self._gq_header(tmp_path)
        assert str(GQ_CEILING) in header
        assert "no-call" in header
