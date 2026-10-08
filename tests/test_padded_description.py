"""The padded consensus FASTA records say where their core lies, and that
survives svirltile.db and `svirlpool get-consensus`."""

import json
from pathlib import Path

from Bio import SeqIO

from svirlpool.localassembly import consensus_align, svirltile
from svirlpool.localassembly.consensus_class import (
    Consensus,
    ConsensusPadding,
    CrsContainerResult,
    parse_padded_description,
)
from svirlpool.scripts.get_consensus_sequences import write_consensus_sequences


def make_consensus(cid="7.0", phasing="phased") -> Consensus:
    core = "acgtacgtac"
    padded = "ttttt" + core.upper() + "gggggggg"
    return Consensus(
        ID=cid,
        crIDs=[7],
        original_regions=[("chr4", 1000, 1500), ("chr4", 2000, 2100)],
        consensus_sequence=core,
        consensus_padding=ConsensusPadding(
            sequence=padded,
            readname_left="r1",
            readname_right="r2",
            padding_size_left=5,
            padding_size_right=8,
            consensus_interval_on_sequence_with_padding=(5, 15),
        ),
        clustering_meta_data={"phasing_status": phasing} if phasing else {},
    )


def test_description_round_trip():
    c = make_consensus()
    description = c.padded_description()
    assert description == "core=5-15 region=chr4:1000-1500 phasing=phased"
    fields = parse_padded_description(f"{c.ID} {description}")
    assert fields == {
        "core_start": 5,
        "core_end": 15,
        "chrom": "chr4",
        "region_start": 1000,
        "region_end": 1500,
        "phasing": "phased",
    }
    padded = c.consensus_padding.sequence
    assert (
        padded[fields["core_start"] : fields["core_end"]].lower()
        == c.consensus_sequence
    )


def test_description_optional_fields():
    c = make_consensus(phasing=None)
    c.original_regions = []
    assert c.padded_description() == "core=5-15"
    assert parse_padded_description("7.0 core=5-15") == {
        "core_start": 5,
        "core_end": 15,
    }
    # a DB built before the field existed stores the bare ID
    assert parse_padded_description("7.0") == {}


def test_contig_names_with_colons():
    fields = parse_padded_description("1.0 core=0-3 region=HLA-A*01:01:01:01:10-20")
    assert fields["chrom"] == "HLA-A*01:01:01:01"
    assert (fields["region_start"], fields["region_end"]) == (10, 20)


def test_padded_fasta_records_carry_the_core(tmp_path: Path):
    """write_consensus_files_for_parallel_processing writes the description."""
    a, b = make_consensus("7.0"), make_consensus("7.1", phasing="single")
    results = tmp_path / "results.jsonl"
    result = CrsContainerResult(consensus_dicts={a.ID: a, b.ID: b}, unused_reads={})
    results.write_text(json.dumps(result.unstructure()) + "\n")
    consensus_paths = {0: tmp_path / "consensus.0.tsv"}
    padded_paths = {0: tmp_path / "padded.0.fasta"}
    consensus_align.write_consensus_files_for_parallel_processing(
        input_consensus_container_results=results,
        consensus_paths=consensus_paths,
        padded_sequences_paths=padded_paths,
    )
    records = {r.id: r for r in SeqIO.parse(padded_paths[0], "fasta")}
    assert set(records) == {"7.0", "7.1"}
    fields = parse_padded_description(records["7.1"].description)
    assert (fields["core_start"], fields["core_end"]) == (5, 15)
    assert fields["phasing"] == "single"
    assert str(records["7.1"].seq)[5:15] == "ACGTACGTAC"


def test_description_survives_svirltile_and_export(tmp_path: Path):
    c = make_consensus()
    fasta = tmp_path / "consensus.fasta"
    fasta.write_text(
        f">{c.ID} {c.padded_description()}\n{c.consensus_padding.sequence}\n"
    )
    db = tmp_path / "svirltile.db"
    svirltile._add_consensus_sequences_to_db(db_path=db, fasta_path=fasta)
    out = tmp_path / "exported.fa"
    write_consensus_sequences(db, out)
    (record,) = SeqIO.parse(out, "fasta")
    assert record.id == c.ID
    fields = parse_padded_description(record.description)
    assert (
        str(record.seq)[fields["core_start"] : fields["core_end"]].lower()
        == c.consensus_sequence
    )
    assert fields["chrom"] == "chr4"
    # the header is written once, not with the ID repeated
    assert out.read_text().splitlines()[0] == f">{c.ID} {c.padded_description()}"
