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

    def test_legacy_flag_does_not_change_alt_supported_calls(self):
        alt = _reads(5, prefix="alt")
        ref = _reads(5, prefix="ref")
        fixed = _genotype_of(alt + ref, alt)
        legacy = _genotype_of(alt + ref, alt, legacy_force_wildtype=True)
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


class TestSyntheticCoverageDoesNotReintroduceTheDefect:
    def test_dummy_covtrees_report_coverage_everywhere(self, dummy_covtrees):
        """Sanity check: the dummy track claims to cover the locus."""
        assert len(dummy_covtrees[SAMPLE][CHR][1000:1500]) == 1

    def test_dummy_coverage_does_not_produce_a_confident_wild_type(
        self, dummy_covtrees
    ):
        """Under ``--skip-covtrees`` the depth is fabricated: refuse to call 0/0."""
        gt = _genotype_of(
            [],
            ["altread"],
            covtrees=dummy_covtrees,
            synthetic_coverage=True,
        )
        assert gt.genotype == "./."
        assert gt.genotype_quality == 0
        assert gt.total_coverage == 0

    def test_without_the_synthetic_marker_the_dummy_track_would_look_covered(
        self, dummy_covtrees
    ):
        """This is exactly the silent degradation the caveat warns about."""
        gt = _genotype_of([], ["altread"], covtrees=dummy_covtrees)
        assert gt.genotype == "0/0"
        assert gt.total_coverage == 1  # the single fabricated interval

    def test_legacy_flag_still_wins_over_the_synthetic_marker(self, dummy_covtrees):
        gt = _genotype_of(
            [],
            ["altread"],
            covtrees=dummy_covtrees,
            synthetic_coverage=True,
            legacy_force_wildtype=True,
        )
        assert gt.genotype == "0/0"
        assert gt.genotype_quality == 60


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
