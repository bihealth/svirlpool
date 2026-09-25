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


def test_max_reads_phases_the_longest_and_leaves_the_rest(tmp_path):
    rng = random.Random(5)
    hap1, hap2 = _haplotypes(rng, n_snvs_every=500)
    reads = _reads({"h1": hap1, "h2": hap2}, n_per_hap=10, rng=rng)
    params = read_phasing.PhasingParams(max_reads=12)
    res = read_phasing.phase_reads(reads, params=params, tmp_dir_path=tmp_path)
    longest = sorted(reads, key=lambda n: (-len(reads[n].seq), n))[:12]
    assert set(res.groups) <= set(longest)
    assert set(reads) - set(longest) <= set(res.unassigned)
    assert res.status == "phased"
    for c in res.clusters().values():
        assert len({rn.split("_")[0] for rn in c}) == 1


@pytest.mark.parametrize("reverse", [False, True])
def test_inverted_alignment_matches_the_read_bases(tmp_path, reverse):
    """parse_pair_both: the inverted view (t on q) puts every difference at q's
    forward coordinate with t's base, also for reverse-strand pairs."""
    rng = random.Random(4)
    hap = "".join(rng.choice("ACGT") for _ in range(8_000))
    a = list(hap)
    a[3000] = {"A": "C", "C": "G", "G": "T", "T": "A"}[a[3000]]
    a = "".join(a[:5000] + a[5030:])  # SNV at 3000, 30 bp deletion at 5000
    b = hap
    if reverse:
        b = b.translate(read_phasing._RC)[::-1]
    reads = {
        n: SeqRecord(Seq(s), id=n, name=n, description="")
        for n, s in (("a", a), ("b", b))
    }
    seqs = {"a": a, "b": b}
    alns = read_phasing.run_ava(
        reads, tmp_path, read_phasing.PhasingParams(), threads=1, timeout=60
    )
    assert len(alns) == 1
    fwd, inv = read_phasing.parse_pair_both(alns[0], seqs)
    assert {fwd.t, inv.t} == {"a", "b"}
    for view in (fwd, inv):
        # "a" lacks 30 bp of "b"
        net = sum(s for _, s in view.indels)
        assert net == (30 if view.t == "a" else -30)
        assert view.t_start == 0 and view.t_end == len(seqs[view.t])
    # "a" as the target: the SNV sits at a's position 3000 with hap's base
    view_a = fwd if fwd.t == "a" else inv
    assert view_a.mism == {3000: hap[3000]}
    # "b" as the target: at b's forward coordinate of hap position 3000
    view_b = fwd if fwd.t == "b" else inv
    pos_b = len(hap) - 1 - 3000 if reverse else 3000
    base_b = a[3000].translate(read_phasing._RC) if reverse else a[3000]
    assert view_b.mism == {pos_b: base_b}


def test_too_few_reads_is_no_information(tmp_path):
    rng = random.Random(3)
    hap1, hap2 = _haplotypes(rng, n_snvs_every=500)
    reads = _reads({"h1": hap1, "h2": hap2}, n_per_hap=2, rng=rng)
    res = read_phasing.phase_reads(reads, tmp_dir_path=tmp_path)
    assert res.status == "no_information"
    assert res.groups == {}
