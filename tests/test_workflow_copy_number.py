"""main.smk: the consensus batches do not depend on the copy-number tracks, and
the copy-number plot is not an output of the tracks rule.

Both made every consensus batch rerun when a run resumed (evaluation hg38,
2026-10-09): a failed run's outputs were removed, among them the QC directory;
the missing QC/copy_number_tracks.png reran the tracks rule, the rewritten
copy_number_tracks.bed.gz then made every consensus batch out of date."""

import re
from pathlib import Path

import svirlpool

SMK = (Path(svirlpool.__file__).parent / "workflows" / "main.smk").read_text()


def rule_block(name: str) -> str:
    m = re.search(rf"^rule {name}:\n(.*?)(?=^rule |\Z)", SMK, re.S | re.M)
    assert m, f"rule {name} not found"
    return m.group(1)


def section(block: str, key: str) -> str:
    m = re.search(rf"^    {key}:\n(.*?)(?=^    [a-z_]+:|\Z)", block, re.S | re.M)
    return m.group(1) if m else ""


def test_consensus_does_not_read_the_copy_number_tracks():
    block = rule_block("consensus_consensus")
    assert "copy_number" not in section(block, "input")
    assert "-cn " not in section(block, "shell")


def test_copy_number_plot_is_its_own_rule():
    tracks = section(rule_block("signalprocessing_generate_copynumber_tracks"), "output")
    assert "copy_number_tracks.bed.gz" in tracks
    assert ".png" not in tracks
    plot = rule_block("QC_copy_number_tracks")
    assert "QC/copy_number_tracks.png" in section(plot, "output")
    assert "copy_number_tracks.db" in section(plot, "input")
