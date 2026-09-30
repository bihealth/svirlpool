"""Stored data written before the size distortions were removed must still load.

``SVpattern.size_distortions`` and ``Consensus.cut_read_alignment_signals`` fed
the removed per-read noise model. Call-only benchmark arms re-run ``sv-calling``
on existing svPatterns / svirltile databases, so records that still carry those
attributes have to load with the current classes; the attributes are dropped.
(The consensus side is covered in ``test_consensus_align.py`` with a real
record.)
"""

import pickle
import sqlite3

from test_can_merge_svComposites import _make_svprimitive_ins

from svirlpool.localassembly import SVpatterns

OLD_DISTORTIONS = {"read1": 3.5, "read2": 0.0, "read3": 12.25}


def _insertion() -> SVpatterns.SVpatternInsertion:
    pattern = SVpatterns.SVpatternInsertion(
        SVprimitives=[_make_svprimitive_ins(read_start=0, read_end=500)]
    )
    pattern.set_sequence("ACGT" * 125)
    return pattern


def _old_record(pattern: SVpatterns.SVpatternType) -> dict:
    """The unstructured record as the previous code wrote it."""
    record = SVpatterns.converter.unstructure(pattern)
    assert "size_distortions" not in record["data"]
    record["data"]["size_distortions"] = dict(OLD_DISTORTIONS)
    return record


def _assert_same_pattern(loaded, pattern) -> None:
    assert type(loaded) is type(pattern)
    assert not hasattr(loaded, "size_distortions")
    assert loaded.get_sequence() == pattern.get_sequence()
    assert loaded.get_size() == pattern.get_size()
    assert loaded.get_supporting_reads() == pattern.get_supporting_reads()
    assert SVpatterns.converter.unstructure(loaded) == SVpatterns.converter.unstructure(
        pattern
    )


def test_an_old_unstructured_record_structures():
    pattern = _insertion()
    loaded = SVpatterns.converter.structure(
        _old_record(pattern), SVpatterns.SVpatternType
    )
    _assert_same_pattern(loaded, pattern)


def test_an_old_json_line_structures():
    """The per-partition JSONL that ``write_svPatterns_to_db`` reads."""
    import json

    pattern = _insertion()
    loaded = SVpatterns.SVpattern.from_json(json.dumps(_old_record(pattern)))
    _assert_same_pattern(loaded, pattern)


def test_an_old_svpatterns_db_row_loads(tmp_path):
    """A row as stored in svpatterns.db / svirltile.db by the previous code."""
    pattern = _insertion()
    db = tmp_path / "svpatterns.db"
    SVpatterns.create_svPatterns_db(db)
    with sqlite3.connect(db) as conn:
        conn.execute(
            "INSERT INTO svPatterns (svPatternID, consensusID, crID, svPattern) "
            "VALUES (?, ?, ?, ?)",
            ("0-INS-1.0", "1.0", 1, pickle.dumps(_old_record(pattern))),
        )
        conn.commit()

    loaded = SVpatterns.read_svPatterns_from_db(db)
    assert 1 == len(loaded)
    _assert_same_pattern(loaded[0], pattern)


def _new(cls):
    return cls.__new__(cls)


class _OldPickle:
    """Pickles as an SVpatternInsertion whose state still has ``size_distortions``."""

    def __init__(self, state: dict):
        self.state = state

    def __reduce__(self):
        return (_new, (SVpatterns.SVpatternInsertion,), self.state)


def test_an_old_directly_pickled_object_loads():
    """attrs' slotted ``__setstate__`` skips attributes the class no longer has."""
    pattern = _insertion()
    state = pattern.__getstate__()
    state["size_distortions"] = dict(OLD_DISTORTIONS)

    loaded = pickle.loads(pickle.dumps(_OldPickle(state)))
    _assert_same_pattern(loaded, pattern)
