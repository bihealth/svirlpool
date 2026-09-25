"""Tests for the experimental read phasing (localassembly/read_phasing.py)."""

import random
import shutil

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from svirlpool.localassembly import read_phasing

pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None, reason="minimap2 not available"
)


def _mutate(seq: str, rate: float, rng: random.Random) -> str:
    out = []
    for b in seq:
        r = rng.random()
        if r < rate / 2:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
        elif r < rate * 3 / 4:
            continue  # deletion
        elif r < rate:
            out.append(b + rng.choice("ACGT"))  # insertion
        else:
            out.append(b)
    return "".join(out)


def _haplotypes(rng: random.Random, n_snvs_every: int | None) -> tuple[str, str]:
    hap1 = "".join(rng.choice("ACGT") for _ in range(12_000))
    hap2 = list(hap1)
    if n_snvs_every:
        for pos in range(250, len(hap1), n_snvs_every):
            hap2[pos] = {"A": "C", "C": "G", "G": "T", "T": "A"}[hap2[pos]]
    return hap1, "".join(hap2)


def _reads(
    haps: dict[str, str], n_per_hap: int, rng: random.Random
) -> dict[str, SeqRecord]:
    reads = {}
    for name, hap in haps.items():
        for i in range(n_per_hap):
            rn = f"{name}_{i}"
            reads[rn] = SeqRecord(
                Seq(_mutate(hap, 0.005, rng)), id=rn, name=rn, description=""
            )
    return reads


def test_two_haplotypes_are_separated(tmp_path):
    rng = random.Random(1)
    hap1, hap2 = _haplotypes(rng, n_snvs_every=500)
    reads = _reads({"h1": hap1, "h2": hap2}, n_per_hap=8, rng=rng)
    res = read_phasing.phase_reads(reads, tmp_dir_path=tmp_path)
    assert res.status == "phased"
    assert res.n_alleles == 2
    clusters = [set(c) for c in res.clusters().values()]
    for c in clusters:
        # every cluster is haplotype-pure
        assert len({rn.split("_")[0] for rn in c}) == 1
    assert sum(len(c) for c in clusters) >= 14


def test_identical_haplotypes_are_not_split(tmp_path):
    rng = random.Random(2)
    hap1, _ = _haplotypes(rng, n_snvs_every=None)
    reads = _reads({"h1": hap1, "h2": hap1}, n_per_hap=8, rng=rng)
    res = read_phasing.phase_reads(reads, tmp_dir_path=tmp_path)
    assert res.n_alleles <= 1
    assert res.status in ("single", "no_information")


def test_too_few_reads_is_no_information(tmp_path):
    rng = random.Random(3)
    hap1, hap2 = _haplotypes(rng, n_snvs_every=500)
    reads = _reads({"h1": hap1, "h2": hap2}, n_per_hap=2, rng=rng)
    res = read_phasing.phase_reads(reads, tmp_dir_path=tmp_path)
    assert res.status == "no_information"
    assert res.groups == {}
