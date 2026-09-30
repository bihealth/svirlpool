"""The KMeans gate (consensus.kmeans_partition) and --kmeans-fast-path-min-k."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from svirlpool.localassembly import consensus


def _pool(names, length=5000):
    return {rn: SeqRecord(Seq("A" * length), id=rn, name=rn) for rn in names}


def _gate(indels, max_k=2):
    return consensus.kmeans_partition(
        dict_summed_indels=indels,
        pool=_pool(indels),
        max_k=max_k,
        variance_threshold=29.0,
        distance_threshold=29.0,
    )


def test_a_homogeneous_pool_is_one_cluster():
    indels = {f"r{i}": [300 + i % 3, 0] for i in range(12)}
    names, labels, k = _gate(indels)
    assert k == 1
    assert sorted(names) == sorted(indels)
    assert set(labels) == {0}


def test_two_separated_alleles_are_two_clusters():
    indels = {f"a{i}": [300 + i % 3, 0] for i in range(8)}
    indels |= {f"b{i}": [0 + i % 3, 0] for i in range(8)}
    names, labels, k = _gate(indels)
    assert k == 2
    by_read = dict(zip(names, labels, strict=True))
    assert len({by_read[f"a{i}"] for i in range(8)}) == 1
    assert len({by_read[f"b{i}"] for i in range(8)}) == 1
    assert by_read["a0"] != by_read["b0"]


def test_a_diffuse_pool_is_rejected():
    indels = {f"r{i}": [i * 40, (i * 97) % 400] for i in range(12)}
    assert _gate(indels) is None

