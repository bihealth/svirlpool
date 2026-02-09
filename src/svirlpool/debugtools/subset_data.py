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
from ..scripts import cut_reads_from_alns
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


def collect_cut_reads_from_regions(
    alignments_bam: Path, regions: list[Region]
) -> dict[str, SeqRecord]:
    """Extract and cut reads to the specified regions using cut_reads_from_alns."""
    read_sequences: dict[str, SeqRecord] = {}
    for region in regions:
        region_tuple = (region.chrom, region.start, region.end)
        cut_reads = cut_reads_from_alns.cut_reads_from_alignments(
            alignments=alignments_bam,
            region=region_tuple,
            buffer_clipped_length=1000,
        )
        read_sequences.update(cut_reads)
    return read_sequences


def write_reads_fastq(reads: dict[str, SeqRecord], output_fastq: Path) -> None:
    import gzip
    with gzip.open(output_fastq, "wt") as f:
        SeqIO.write(list(reads.values()), f, "fastq")


def build_reduced_test_data(
    alignments: Path,
    reference: Path,
    trf_bed: Path,
    mononucleotides_bed: Path,
    regions_bed: Path,
    output_dir: Path,
    threads: int,
    subsample: float = 1.0,
    seed: int = 0,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    regions = parse_bed_regions(regions_bed)
    reduced_reference = output_dir / "reference.reduced.fasta.gz"
    lifted_trf = output_dir / "trf.reduced.bed"
    lifted_mono = output_dir / "mononucleotides.reduced.bed"
    reads_fastq = output_dir / "reads.reduced.fastq.gz"
    realigned_bam = output_dir / "alignments.reduced.bam"

    # Create reduced reference in temp file, then bgzip and index
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp_ref:
        tmp_ref_path = Path(tmp_ref.name)
        create_reduced_reference(reference, regions, tmp_ref_path)
    
    # Bgzip the reference
    subprocess.run(
        ["bgzip", "-c", str(tmp_ref_path)],
        stdout=open(reduced_reference, "wb"),
        check=True,
    )
    
    # Remove the temporary uncompressed file
    tmp_ref_path.unlink()
    
    # Index the bgzipped reference with samtools faidx
    subprocess.run(
        ["samtools", "faidx", str(reduced_reference)],
        check=True,
        capture_output=True,
    )
    
    lift_bed_annotations(trf_bed, regions, lifted_trf)
    lift_bed_annotations(mononucleotides_bed, regions, lifted_mono)

    read_sequences = collect_cut_reads_from_regions(alignments, regions)
    
    # Subsample reads if requested
    if subsample < 1.0:
        import random
        random.seed(seed)
        n_reads = len(read_sequences)
        n_subsample = int(n_reads * subsample)
        read_names = list(read_sequences.keys())
        random.shuffle(read_names)
        subsampled_names = set(read_names[:n_subsample])
        read_sequences = {name: seq for name, seq in read_sequences.items() if name in subsampled_names}
    
    write_reads_fastq(read_sequences, reads_fastq)

    util.align_reads_with_minimap(
        reference=reduced_reference,
        reads=reads_fastq,
        bamout=realigned_bam,
        threads=threads,
        tech="map-ont",
    )
    
    # Remove the reads fastq file after alignment to save space
    reads_fastq.unlink()


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
    parser.add_argument(
        "--subsample",
        type=float,
        default=1.0,
        help="Subsample reads to this fraction (0.0-1.0). Default: 1.0 (no subsampling)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for subsampling reproducibility. Default: 0",
    )
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
        subsample=args.subsample,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
