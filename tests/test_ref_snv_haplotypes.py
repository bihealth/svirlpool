"""ref_snv_haplotypes: the cheap haplotype split from het SNVs in the reads'
reference alignments."""

import random

import pysam
import pytest

from svirlpool.localassembly.ref_snv_haplotypes import ref_snv_haplotypes

REF_LEN = 3000
SNV_POS = (400, 900, 1300, 1800, 2200, 2600)


def _random_seq(rng, n):
    # no homopolymers >= 4 so that every SNV column is eligible
    seq = []
    while len(seq) < n:
        b = rng.choice("ACGT")
        if len(seq) >= 3 and seq[-1] == seq[-2] == seq[-3] == b:
            continue
        seq.append(b)
    return "".join(seq)


@pytest.fixture
def reference(tmp_path):
    rng = random.Random(3)
    seq = _random_seq(rng, REF_LEN)
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chrT\n{seq}\n")
    pysam.faidx(str(fa))
    with pysam.FastaFile(str(fa)) as f:
        yield f, seq


def _hap_b(seq):
    s = list(seq)
    for p in SNV_POS:
        s[p] = next(b for b in "ACGT" if b != s[p] and b not in (s[p - 1], s[p + 1]))
    return "".join(s)


def _reads(seq, prefix, n, rng, error_rate=0.01):
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chrT", "LN": REF_LEN}]})
    out = []
    for i in range(n):
        read = list(seq)
        for p in range(REF_LEN):  # private errors: never at the same column twice
            if rng.random() < error_rate:
                read[p] = rng.choice([b for b in "ACGT" if b != read[p]])
        a = pysam.AlignedSegment(header)
        a.query_name = f"{prefix}{i}"
        a.query_sequence = "".join(read)
        a.reference_id = 0
        a.reference_start = 0
        a.cigartuples = [(0, REF_LEN)]
        a.mapping_quality = 60
        out.append(a)
    return out


def test_two_haplotypes_split_by_their_snvs(reference):
    ref, seq = reference
    rng = random.Random(1)
    alns = _reads(seq, "a", 7, rng) + _reads(_hap_b(seq), "b", 7, rng)
    split = ref_snv_haplotypes(
        alns=alns, reads={a.query_name for a in alns}, ref=ref,
        chrom="chrT", start=1400, end=1600,
    )  # fmt: skip
    assert split.is_split
    assert split.n_sites == len(SNV_POS)
    haps = split.haplotypes
    assert len({haps[f"a{i}"] for i in range(7)}) == 1
    assert len({haps[f"b{i}"] for i in range(7)}) == 1
    assert haps["a0"] != haps["b0"]


def test_one_haplotype_does_not_split(reference):
    ref, seq = reference
    rng = random.Random(2)
    alns = _reads(seq, "a", 14, rng)
    split = ref_snv_haplotypes(
        alns=alns, reads={a.query_name for a in alns}, ref=ref,
        chrom="chrT", start=1400, end=1600,
    )  # fmt: skip
    assert not split.is_split
    assert split.n_sites == 0


def test_sequencing_errors_alone_are_no_sites(reference):
    ref, seq = reference
    rng = random.Random(4)
    alns = _reads(seq, "a", 14, rng, error_rate=0.05)
    split = ref_snv_haplotypes(
        alns=alns, reads={a.query_name for a in alns}, ref=ref,
        chrom="chrT", start=1400, end=1600,
    )  # fmt: skip
    assert not split.is_split
