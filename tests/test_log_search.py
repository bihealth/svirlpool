"""
Tests for log_search module.

Tests the extraction of crIDs and regions from log lines, and the matching logic.
"""

import gzip
from pathlib import Path

import pytest

from svirlpool.debugtools.log_search import (
    Region,
    extract_crIDs_from_line,
    extract_regions_from_line,
    is_relevant_line,
    line_matches,
    parse_region,
    search_log,
)

# ---- Region parsing ----


class TestParseRegion:
    def test_simple_region(self):
        r = parse_region("chr1:100-200")
        assert r.chrom == "chr1"
        assert r.start == 100
        assert r.end == 200

    def test_region_with_commas(self):
        r = parse_region("chr1:1,000-2,000")
        assert r.chrom == "chr1"
        assert r.start == 1000
        assert r.end == 2000

    def test_complex_chrom_name(self):
        r = parse_region("R8_14_91775678_92309099:247964-248188")
        assert r.chrom == "R8_14_91775678_92309099"
        assert r.start == 247964
        assert r.end == 248188

    def test_invalid_no_colon(self):
        with pytest.raises(ValueError, match="no colon"):
            parse_region("chr1_100_200")

    def test_invalid_no_dash(self):
        with pytest.raises(ValueError, match="no dash"):
            parse_region("chr1:100200")


# ---- Region overlaps ----


class TestRegionOverlaps:
    def test_overlapping_regions(self):
        a = Region("chr1", 100, 300)
        b = Region("chr1", 200, 400)
        assert a.overlaps(b)
        assert b.overlaps(a)

    def test_non_overlapping_same_chrom(self):
        a = Region("chr1", 100, 200)
        b = Region("chr1", 300, 400)
        assert not a.overlaps(b)
        assert not b.overlaps(a)

    def test_adjacent_non_overlapping(self):
        a = Region("chr1", 100, 200)
        b = Region("chr1", 200, 300)
        assert not a.overlaps(b)

    def test_contained_region(self):
        a = Region("chr1", 100, 500)
        b = Region("chr1", 200, 300)
        assert a.overlaps(b)
        assert b.overlaps(a)

    def test_different_chromosomes(self):
        a = Region("chr1", 100, 300)
        b = Region("chr2", 100, 300)
        assert not a.overlaps(b)


# ---- crID extraction ----


class TestExtractCRIDs:
    def test_single_crID(self):
        line = "DROPPED::some_function crID=42 something"
        assert extract_crIDs_from_line(line) == {42}

    def test_crIDs_curly_braces(self):
        line = "TRANSFORMED:: crIDs={1,2,3} more stuff"
        assert extract_crIDs_from_line(line) == {1, 2, 3}

    def test_crIDs_square_brackets(self):
        line = "TRANSFORMED:: crIDs=[5, 19, 42] more stuff"
        assert extract_crIDs_from_line(line) == {5, 19, 42}

    def test_crIDs_curly_with_spaces(self):
        line = "crIDs={10, 20, 30}"
        assert extract_crIDs_from_line(line) == {10, 20, 30}

    def test_multiple_crID_fields(self):
        line = "crID=5 something crID=19 and crIDs={42,99}"
        assert extract_crIDs_from_line(line) == {5, 19, 42, 99}

    def test_svpattern_log_id_format(self):
        line = "TRANSFORMED::func: sample=S1|consensusID=19.0|crID=19|type=DEL|size=100|region=chr1:1000-2000"
        assert extract_crIDs_from_line(line) == {19}

    def test_svcomposite_log_id_format(self):
        line = (
            "DEL|crIDs={5,19}|regions=chr1:1000-2000|consensusIDs=['S1:5.0', 'S1:19.0']"
        )
        assert extract_crIDs_from_line(line) == {5, 19}

    def test_no_crIDs(self):
        line = "DROPPED::some function with no cr IDs"
        assert extract_crIDs_from_line(line) == set()

    def test_consensusID_with_crID(self):
        # crID should be extracted from crID= field, not from consensusID=19.0
        line = "consensusID=19.0|crID=19|type=INS"
        crIDs = extract_crIDs_from_line(line)
        assert 19 in crIDs


# ---- Region extraction ----


class TestExtractRegions:
    def test_single_region(self):
        line = "TRANSFORMED:: region=chr1:1000-2000 other stuff"
        regions = extract_regions_from_line(line)
        assert len(regions) == 1
        assert regions[0] == Region("chr1", 1000, 2000)

    def test_regions_semicolon_separated(self):
        line = "regions=chr1:100-200;chr2:300-400"
        regions = extract_regions_from_line(line)
        assert len(regions) == 2
        assert Region("chr1", 100, 200) in regions
        assert Region("chr2", 300, 400) in regions

    def test_svpattern_log_id_region(self):
        line = (
            "sample=S1|consensusID=19.0|crID=19|type=DEL|size=100|region=chr1:1000-2000"
        )
        regions = extract_regions_from_line(line)
        assert len(regions) == 1
        assert regions[0] == Region("chr1", 1000, 2000)

    def test_svcomposite_regions(self):
        line = "DEL|crIDs={5}|regions=chr1:1000-2000;chr1:3000-4000|consensusIDs=['S1:5.0']"
        regions = extract_regions_from_line(line)
        assert len(regions) == 2

    def test_no_regions(self):
        line = "DROPPED::some function with no regions"
        regions = extract_regions_from_line(line)
        assert len(regions) == 0

    def test_multiple_region_fields(self):
        line = "region=chr1:100-200 region=chr2:300-400"
        regions = extract_regions_from_line(line)
        assert len(regions) == 2

    def test_genotype_result_format(self):
        line = "GENOTYPE|RESULT    sample=S1   region=chr1:1000-2000    GT=0/1"
        regions = extract_regions_from_line(line)
        assert len(regions) == 1
        assert regions[0] == Region("chr1", 1000, 2000)


# ---- is_relevant_line ----


class TestIsRelevantLine:
    def test_dropped_line(self):
        assert is_relevant_line("2025-01-01 - DEBUG - DROPPED::some_function::reason")

    def test_transformed_line(self):
        assert is_relevant_line(
            "2025-01-01 - DEBUG - TRANSFORMED::merge_insertions::vertical_merge"
        )

    def test_irrelevant_line(self):
        assert not is_relevant_line("2025-01-01 - INFO - Processing sample S1")

    def test_empty_line(self):
        assert not is_relevant_line("")


# ---- line_matches ----


class TestLineMatches:
    def test_crID_match(self):
        line = "DROPPED:: crID=42 region=chr1:100-200"
        assert line_matches(line, {42}, [])
        assert not line_matches(line, {99}, [])

    def test_region_overlap_match(self):
        line = "TRANSFORMED:: region=chr1:100-200"
        query_region = Region("chr1", 150, 250)
        assert line_matches(line, set(), [query_region])

    def test_region_no_overlap(self):
        line = "TRANSFORMED:: region=chr1:100-200"
        query_region = Region("chr1", 300, 400)
        assert not line_matches(line, set(), [query_region])

    def test_crID_or_region_match(self):
        line = "TRANSFORMED:: crID=5 region=chr1:100-200"
        # Match by crID even if region doesn't match
        assert line_matches(line, {5}, [Region("chr2", 100, 200)])
        # Match by region even if crID doesn't match
        assert line_matches(line, {99}, [Region("chr1", 150, 250)])

    def test_no_queries(self):
        line = "DROPPED:: crID=42"
        assert not line_matches(line, set(), [])

    def test_svcomposite_log_id_match(self):
        line = (
            "TRANSFORMED::merge_deletions::vertical_merge:(to merged SVcomposite)"
            "svComposites=[DEL|crIDs={5,19}|regions=chr1:1000-2000|consensusIDs=['S1:5.0', 'S1:19.0']] "
            "-->   DEL|crIDs={5,19}|regions=chr1:1000-2000|consensusIDs=['S1:5.0', 'S1:19.0']"
        )
        assert line_matches(line, {5}, [])
        assert line_matches(line, {19}, [])
        assert not line_matches(line, {99}, [])
        assert line_matches(line, set(), [Region("chr1", 1500, 2500)])

    def test_real_dropped_line(self):
        line = (
            "2025-01-01 10:00:00 - svirlpool.svcalling.multisample_sv_calling - DEBUG - "
            "DROPPED::multisample_sv_calling::MIN SV SIZE NOT REACHED: 25 < min_sv_size 30, "
            "svComposite=INS|crIDs={42}|regions=chr3:5000-6000|consensusIDs=['S1:42.0']"
        )
        assert line_matches(line, {42}, [])
        assert line_matches(line, set(), [Region("chr3", 5500, 7000)])
        assert not line_matches(line, {1}, [])


# ---- search_log ----


class TestSearchLog:
    def _write_log(self, lines: list[str], path: Path) -> None:
        with open(path, "w") as f:
            for line in lines:
                f.write(line + "\n")

    def _write_log_gz(self, lines: list[str], path: Path) -> None:
        with gzip.open(path, "wt") as f:
            for line in lines:
                f.write(line + "\n")

    def test_search_by_crID(self, tmp_path):
        log_lines = [
            "INFO - Starting processing",
            "DEBUG - DROPPED::func crID=5 region=chr1:100-200",
            "DEBUG - TRANSFORMED::func crID=10 region=chr2:300-400",
            "DEBUG - DROPPED::func crID=20 region=chr3:500-600",
            "INFO - Done",
        ]
        log_file = tmp_path / "test.log"
        self._write_log(log_lines, log_file)

        matches = search_log(log_file, {5}, [])
        assert len(matches) == 1
        assert "crID=5" in matches[0]

    def test_search_by_region(self, tmp_path):
        log_lines = [
            "DEBUG - DROPPED::func crID=5 region=chr1:100-200",
            "DEBUG - TRANSFORMED::func crID=10 region=chr2:300-400",
        ]
        log_file = tmp_path / "test.log"
        self._write_log(log_lines, log_file)

        matches = search_log(log_file, set(), [Region("chr1", 150, 250)])
        assert len(matches) == 1
        assert "crID=5" in matches[0]

    def test_search_gzipped(self, tmp_path):
        log_lines = [
            "DEBUG - DROPPED::func crID=5 region=chr1:100-200",
            "DEBUG - TRANSFORMED::func crID=10 region=chr2:300-400",
        ]
        log_file = tmp_path / "test.log.gz"
        self._write_log_gz(log_lines, log_file)

        matches = search_log(log_file, {10}, [])
        assert len(matches) == 1
        assert "crID=10" in matches[0]

    def test_search_show_all(self, tmp_path):
        log_lines = [
            "INFO - Starting",
            "DEBUG - DROPPED::func1 crID=5",
            "DEBUG - TRANSFORMED::func2 crID=10",
            "INFO - Done",
        ]
        log_file = tmp_path / "test.log"
        self._write_log(log_lines, log_file)

        matches = search_log(log_file, set(), [], show_all_relevant=True)
        assert len(matches) == 2

    def test_search_no_matches(self, tmp_path):
        log_lines = [
            "DEBUG - DROPPED::func crID=5 region=chr1:100-200",
            "INFO - Processing",
        ]
        log_file = tmp_path / "test.log"
        self._write_log(log_lines, log_file)

        matches = search_log(log_file, {999}, [Region("chrX", 1, 2)])
        assert len(matches) == 0

    def test_search_multiple_crIDs(self, tmp_path):
        log_lines = [
            "DEBUG - DROPPED::func crID=5 region=chr1:100-200",
            "DEBUG - TRANSFORMED::func crID=10 region=chr2:300-400",
            "DEBUG - DROPPED::func crID=20 region=chr3:500-600",
        ]
        log_file = tmp_path / "test.log"
        self._write_log(log_lines, log_file)

        matches = search_log(log_file, {5, 20}, [])
        assert len(matches) == 2

    def test_search_combined_crID_and_region(self, tmp_path):
        log_lines = [
            "DEBUG - DROPPED::func crID=5 region=chr1:100-200",
            "DEBUG - TRANSFORMED::func crID=10 region=chr2:300-400",
            "DEBUG - DROPPED::func crID=20 region=chr3:500-600",
        ]
        log_file = tmp_path / "test.log"
        self._write_log(log_lines, log_file)

        # Match crID=5 OR region overlapping chr2:350-450
        matches = search_log(log_file, {5}, [Region("chr2", 350, 450)])
        assert len(matches) == 2
        assert any("crID=5" in m for m in matches)
        assert any("crID=10" in m for m in matches)

    def test_realistic_log_lines(self, tmp_path):
        log_lines = [
            "2025-01-01 10:00:00 - svirlpool - DEBUG - TRANSFORMED::svPatterns_to_horizontally_merged_svComposites::from_SVpattern: sample=S1|consensusID=19.0|crID=19|type=DEL|size=100|region=chr1:1000-2000 -> DEL|crIDs={19}|regions=chr1:1000-2000|consensusIDs=['S1:19.0']",
            "2025-01-01 10:00:01 - svirlpool - DEBUG - DROPPED::multisample_sv_calling::MIN SV SIZE NOT REACHED: 25 < min_sv_size 30, svComposite=INS|crIDs={42}|regions=chr3:5000-6000|consensusIDs=['S1:42.0']",
            "2025-01-01 10:00:02 - svirlpool - DEBUG - TRANSFORMED::merge_adjacencies::vertical_merge:(to merged SVcomposite) crIDs=[7, 8, 9] consensusIDs=['S1:7.0', 'S1:8.0', 'S1:9.0'] n_merged=3 regions=chr2:100-500;chr2:600-900",
            "2025-01-01 10:00:03 - svirlpool - INFO - Processing complete",
        ]
        log_file = tmp_path / "test.log"
        self._write_log(log_lines, log_file)

        # Search for crID=19
        matches = search_log(log_file, {19}, [])
        assert len(matches) == 1
        assert "crID=19" in matches[0]

        # Search for crID=42
        matches = search_log(log_file, {42}, [])
        assert len(matches) == 1
        assert "MIN SV SIZE" in matches[0]

        # Search for crID=8 (from crIDs=[7, 8, 9])
        matches = search_log(log_file, {8}, [])
        assert len(matches) == 1
        assert "merge_adjacencies" in matches[0]

        # Search by region overlapping chr2:150-250
        matches = search_log(log_file, set(), [Region("chr2", 150, 250)])
        assert len(matches) == 1
        assert "merge_adjacencies" in matches[0]

        # Search by region overlapping chr3:5500-7000
        matches = search_log(log_file, set(), [Region("chr3", 5500, 7000)])
        assert len(matches) == 1
        assert "crIDs={42}" in matches[0]
