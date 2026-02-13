#!/usr/bin/env python
"""
Build a reduced reference genome and lifted annotation beds for end-to-end testing.

Steps:
1) Use source regions to create a new reference FASTA (each region becomes a contig).
2) Intersect TRF and mononucleotide beds with the source regions and lift coordinates.
3) Extract reads aligned to the source regions and realign them to the reduced reference.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..localassembly import consensus
from ..util import util


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    name: str


def parse_bed_regions(path: Path) -> list[Region]:
    regions: list[Region] = []
    with open(path, "r") as f:
        for idx, line in enumerate(f):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(f"Invalid BED line: {line.rstrip()}")
            chrom, start_s, end_s = fields[0], fields[1], fields[2]
            start, end = int(start_s), int(end_s)
            if end <= start:
                raise ValueError(f"Invalid BED interval: {line.rstrip()}")
            name = fields[3] if len(fields) > 3 else f"R{idx+1}_{chrom}_{start}_{end}"
            regions.append(Region(chrom=chrom, start=start, end=end, name=name))
    if not regions:
        raise ValueError(f"No regions found in {path}")
    return regions


def write_regions_bed(path: Path, regions: Iterable[Region]) -> None:
    with open(path, "w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        for region in regions:
            writer.writerow([region.chrom, region.start, region.end, region.name])


def create_reduced_reference(
    reference_fasta: Path, regions: list[Region], output_fasta: Path
) -> dict[str, Region]:
    mapping: dict[str, Region] = {}
    with pysam.FastaFile(str(reference_fasta)) as ref, open(output_fasta, "w") as out:
        for region in regions:
            seq = ref.fetch(region.chrom, region.start, region.end)
            record = SeqRecord(
                Seq(seq),
                id=region.name,
                name=region.name,
                description=f"{region.chrom}:{region.start}-{region.end}",
            )
            SeqIO.write(record, out, "fasta")
            mapping[region.name] = region
    return mapping


def run_bedtools_intersect(input_bed: Path, regions_bed: Path) -> list[list[str]]:
    cmd = [
        "bedtools",
        "intersect",
        "-wa",
        "-wb",
        "-a",
        str(input_bed),
        "-b",
        str(regions_bed),
    ]
    try:
        result = subprocess.run(
            cmd, check=True, capture_output=True, text=True
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "bedtools not found. Please install bedtools and ensure it is in PATH."
        ) from exc
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    return [line.split("\t") for line in lines]


def lift_bed_annotations(
    input_bed: Path,
    regions: list[Region],
    output_bed: Path,
) -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        regions_bed = Path(tmp_dir) / "regions.bed"
        write_regions_bed(regions_bed, regions)
        intersect_rows = run_bedtools_intersect(input_bed, regions_bed)

    region_lookup = {region.name: region for region in regions}

    lifted_rows: list[list[str]] = []
    for row in intersect_rows:
        if len(row) < 7:
            continue
        # A (input bed) has at least 3 columns, B (regions) has 4
        # Layout: A_cols..., B_chrom, B_start, B_end, B_name
        b_chrom, b_start_s, b_end_s, b_name = row[-4:]
        region = region_lookup.get(b_name)
        if region is None:
            continue
        feature_chrom = row[0]
        feature_start = int(row[1])
        feature_end = int(row[2])
        if feature_chrom != region.chrom:
            continue
        clipped_start = max(feature_start, region.start)
        clipped_end = min(feature_end, region.end)
        if clipped_end <= clipped_start:
            continue
        new_start = clipped_start - region.start
        new_end = clipped_end - region.start
        new_row = [region.name, str(new_start), str(new_end)] + row[3:-4]
        lifted_rows.append(new_row)

    with open(output_bed, "w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        for row in lifted_rows:
            writer.writerow(row)


def collect_alignments_for_regions(
    alignments_bam: Path, regions: list[Region]
) -> dict[int, list[pysam.AlignedSegment]]:
    dict_alignments: dict[int, list[pysam.AlignedSegment]] = {}
    with pysam.AlignmentFile(str(alignments_bam), "rb") as bam:
        for idx, region in enumerate(regions, start=1):
            alignments = [
                aln
                for aln in bam.fetch(region.chrom, region.start, region.end)
                if not aln.is_secondary
            ]
            dict_alignments[idx] = alignments
    return dict_alignments


def write_reads_fasta(reads: dict[str, SeqRecord], output_fasta: Path) -> None:
    with open(output_fasta, "w") as f:
        SeqIO.write(list(reads.values()), f, "fasta")


def build_reduced_test_data(
    alignments: Path,
    reference: Path,
    trf_bed: Path,
    mononucleotides_bed: Path,
    regions_bed: Path,
    output_dir: Path,
    threads: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    regions = parse_bed_regions(regions_bed)
    reduced_reference = output_dir / "reference.reduced.fasta"
    lifted_trf = output_dir / "trf.reduced.bed"
    lifted_mono = output_dir / "mononucleotides.reduced.bed"
    reads_fasta = output_dir / "reads.reduced.fasta"
    realigned_bam = output_dir / "alignments.reduced.bam"

    create_reduced_reference(reference, regions, reduced_reference)
    lift_bed_annotations(trf_bed, regions, lifted_trf)
    lift_bed_annotations(mononucleotides_bed, regions, lifted_mono)

    dict_alignments = collect_alignments_for_regions(alignments, regions)
    read_sequences = consensus.get_full_read_sequences_of_alignments(
        dict_alignments=dict_alignments, path_alignments=alignments
    )
    write_reads_fasta(read_sequences, reads_fasta)

    util.align_reads_with_minimap(
        reference=reduced_reference,
        reads=reads_fasta,
        bamout=realigned_bam,
        threads=threads,
        tech="map-ont",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a reduced reference genome and lift annotations for testing.",
    )
    parser.add_argument("--alignments", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--trf", required=True, type=Path)
    parser.add_argument("--mononucleotides", required=True, type=Path)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=1)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    build_reduced_test_data(
        alignments=args.alignments,
        reference=args.reference,
        trf_bed=args.trf,
        mononucleotides_bed=args.mononucleotides,
        regions_bed=args.regions,
        output_dir=args.output_dir,
        threads=args.threads,
    )


if __name__ == "__main__":
    main()
