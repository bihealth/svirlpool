"""The KMeans gate (consensus.kmeans_partition) of --clustering-strategy balanced / fast."""

import pytest
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



_REQUIRED = ["-s", "S", "-i", "i.db", "-a", "a.bam", "-cn", "c.bed.gz", "-o", "o.jsonl",
             "-r", "r.fa"]  # fmt: skip


def test_the_clustering_strategy_defaults_to_accurate():
    parser = consensus.get_consensus_parser()
    assert parser.parse_args(_REQUIRED).clustering_strategy == "accurate"
    for strategy in ("balanced", "fast"):
        args = parser.parse_args(_REQUIRED + ["--clustering-strategy", strategy])
        assert args.clustering_strategy == strategy


@pytest.mark.parametrize("value", ["off", "kmeans", "kmeans+snv"])
def test_the_former_fast_clustering_values_are_refused(value):
    with pytest.raises(SystemExit):
        consensus.get_consensus_parser().parse_args(
            _REQUIRED + ["--clustering-strategy", value]
        )
