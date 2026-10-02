"""Tests for the reference-based ablation arm of the read phasing
(localassembly/ref_read_phasing.py)."""

import random
import shutil
import subprocess

import pysam
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from svirlpool.localassembly import read_phasing, ref_read_phasing

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
            continue
        elif r < rate:
            out.append(b + rng.choice("ACGT"))
        else:
            out.append(b)
    return "".join(out)


def _setup(tmp_path, hap2_snvs_every, hap2_deletion, seed=1, n_per_hap=8):
    rng = random.Random(seed)
    ref = "".join(rng.choice("ACGT") for _ in range(14_000))
    hap1 = ref
    h2 = list(ref)
    if hap2_snvs_every:
        for pos in range(300, len(h2), hap2_snvs_every):
            h2[pos] = {"A": "C", "C": "G", "G": "T", "T": "A"}[h2[pos]]
    hap2 = "".join(h2)
    if hap2_deletion:
        hap2 = hap2[:7000] + hap2[7300:]
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chrT\n{ref}\n")
    pysam.faidx(str(fa))
    reads = {}
    for name, hap in (("h1", hap1), ("h2", hap2)):
        for i in range(n_per_hap):
            rn = f"{name}_{i}"
            s = _mutate(hap[500:-500], 0.005, rng)
            reads[rn] = SeqRecord(Seq(s), id=rn, name=rn, description="")
    rfa = tmp_path / "reads.fa"
    rfa.write_text("".join(f">{n}\n{r.seq}\n" for n, r in reads.items()))
    sam = tmp_path / "aln.sam"
    with open(sam, "w") as f:
        subprocess.run(["minimap2", "-a", "-x", "map-ont", str(fa), str(rfa)],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)
    with pysam.AlignmentFile(str(sam)) as f:
        alns = [a for a in f if not a.is_unmapped]
    return reads, alns, pysam.FastaFile(str(fa))


def _phase(reads, alns, ref, params=None):
    return ref_read_phasing.phase_reads_reference(
        reads, alns, ref, "chrT", 6900, 7400, flank=10_000, params=params
    )


def test_snv_haplotypes_are_separated(tmp_path):
    res = _phase(*_setup(tmp_path, hap2_snvs_every=500, hap2_deletion=False))
    assert res.status == "phased"
    assert res.n_alleles == 2
    for members in res.clusters().values():
        assert len({m.split("_")[0] for m in members}) == 1


def test_identical_haplotypes_are_not_split(tmp_path):
    res = _phase(*_setup(tmp_path, hap2_snvs_every=None, hap2_deletion=False))
    assert res.status in ("single", "no_information")


def test_a_deletion_alone_separates_the_haplotypes(tmp_path):
    # one SV is one discriminating position: the default refinement
    # (min_discriminating = 2) merges the two groups, in read_phasing too
    reads, alns, ref = _setup(tmp_path, hap2_snvs_every=None, hap2_deletion=True)
    assert _phase(reads, alns, ref).status == "single"
    res = _phase(reads, alns, ref, read_phasing.PhasingParams(min_discriminating=1))
    assert res.n_snv_sites == 0
    assert res.n_sv_sites > 0
    assert res.status == "phased"
    for members in res.clusters().values():
        assert len({m.split("_")[0] for m in members}) == 1
