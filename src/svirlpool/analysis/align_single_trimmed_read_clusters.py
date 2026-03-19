import argparse
import logging
import re
import shlex
import subprocess
import tempfile
from pathlib import Path

import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..localassembly import consensus, consensus_class
from ..localassembly.consensus import (
    get_full_read_sequences_of_alignments,
    get_max_extents_of_read_alignments_on_cr,
    get_read_alignment_intervals_in_region,
    trim_reads,
)

log = logging.getLogger(__name__)


def _sanitize_consensus_id(consensus_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", consensus_id)


def _normalize_region(region: tuple[str, int, int]) -> tuple[str, int, int]:
    chrom, start, end = region
    if start <= end:
        return chrom, start, end
    return chrom, end, start


def _find_existing_mmi(reference: Path) -> Path | None:
    exact = reference.parent / f"{reference.name}.mmi"
    if exact.exists() and exact.is_file():
        return exact
    stem_match = reference.parent / f"{reference.stem}.mmi"
    if stem_match.exists() and stem_match.is_file():
        return stem_match
    return None


def _fetch_alignments_from_region(
    alignments_path: Path,
    region: tuple[str, int, int],
    no_secondary: bool = True,
    reference: Path | None = None,
) -> list[pysam.AlignedSegment]:
    result: list[pysam.AlignedSegment] = []
    is_cram = str(alignments_path).endswith(".cram")
    open_mode = "rc" if is_cram else "rb"
    open_kwargs = {}
    if is_cram and reference is not None:
        open_kwargs["reference_filename"] = str(reference)

    with pysam.AlignmentFile(alignments_path, open_mode, **open_kwargs) as aln_file:
        for aln in aln_file.fetch(*region):
            if no_secondary and aln.is_secondary:
                continue
            result.append(aln)
    return result


def _collect_consensus_alignments(
    consensus_object: consensus_class.Consensus,
    alignments_path: Path,
    reference: Path,
) -> dict[int, list[pysam.AlignedSegment]]:
    selected_readnames = consensus_object.get_used_readnames()
    seen = set()
    selected_alignments: list[pysam.AlignedSegment] = []

    for region in consensus_object.original_regions:
        normalized_region = _normalize_region(region)
        for aln in _fetch_alignments_from_region(
            alignments_path=alignments_path,
            region=normalized_region,
            no_secondary=True,
            reference=reference,
        ):
            if aln.query_name not in selected_readnames:
                continue
            key = (
                str(aln.query_name),
                str(aln.reference_name),
                int(aln.reference_start),
                int(aln.reference_end),
                bool(aln.is_reverse),
                str(aln.cigarstring),
            )
            if key in seen:
                continue
            seen.add(key)
            selected_alignments.append(aln)

    return {0: selected_alignments}


def _collect_read_intervals_for_consensus(
    consensus_object: consensus_class.Consensus,
    alignments: list[pysam.AlignedSegment],
    buffer_clipped_length: int,
) -> dict[str, list[tuple[int, int, str, int, int]]]:
    all_intervals: dict[str, list[tuple[int, int, str, int, int]]] = {}
    for region in consensus_object.original_regions:
        chrom, region_start, region_end = _normalize_region(region)
        region_alignments = [
            aln
            for aln in alignments
            if aln.reference_name == chrom
            and aln.reference_start < region_end
            and aln.reference_end > region_start
        ]
        if len(region_alignments) == 0:
            continue
        region_intervals = get_read_alignment_intervals_in_region(
            region_start=region_start,
            regions_end=region_end,
            alignments=region_alignments,
            buffer_clipped_length=buffer_clipped_length,
        )
        for readname, intervals in region_intervals.items():
            all_intervals.setdefault(readname, []).extend(intervals)
    return all_intervals


def _trim_reads_for_consensus(
    consensus_object: consensus_class.Consensus,
    alignments_path: Path,
    reference: Path,
    buffer_clipped_length: int,
) -> dict[str, SeqRecord]:
    dict_alignments = _collect_consensus_alignments(
        consensus_object=consensus_object,
        alignments_path=alignments_path,
        reference=reference,
    )
    flat_alignments = dict_alignments[0]
    if len(flat_alignments) == 0:
        return {}

    intervals = _collect_read_intervals_for_consensus(
        consensus_object=consensus_object,
        alignments=flat_alignments,
        buffer_clipped_length=buffer_clipped_length,
    )
    if len(intervals) == 0:
        return {}

    max_intervals = get_max_extents_of_read_alignments_on_cr(intervals)
    read_records = get_full_read_sequences_of_alignments(
        dict_alignments=dict_alignments,
        path_alignments=alignments_path,
        reference=reference,
    )
    return trim_reads(
        dict_alignments=dict_alignments,
        intervals=max_intervals,
        read_records=read_records,
    )


def _write_trimmed_reads_fastq(
    trimmed_reads: dict[str, SeqRecord], output_fastq: Path
) -> None:
    records: list[SeqRecord] = []
    for record in trimmed_reads.values():
        rec = record[:]
        if not rec.letter_annotations or "phred_quality" not in rec.letter_annotations:
            rec.letter_annotations["phred_quality"] = [30] * len(rec.seq)
        records.append(rec)
    with open(output_fastq, "w") as f:
        SeqIO.write(records, f, "fastq")


def _align_fastq_to_reference(
    reference: Path,
    reads_fastq: Path,
    output_bam: Path,
    threads: int,
) -> None:
    mmi = _find_existing_mmi(reference)
    reference_for_alignment = mmi if mmi is not None else reference
    if mmi is not None:
        log.info(f"Using minimap2 index: {mmi}")
    else:
        log.info("No minimap2 index found in reference parent directory, using FASTA")

    cmd_align = shlex.split(
        f"minimap2 -a -x map-ont -t {threads} {reference_for_alignment} {reads_fastq}"
    )
    cmd_sort = shlex.split(f"samtools sort -O BAM -o {output_bam}")
    cmd_index = shlex.split(f"samtools index {output_bam}")

    p_align = subprocess.Popen(cmd_align, stdout=subprocess.PIPE)
    p_sort = subprocess.Popen(cmd_sort, stdin=p_align.stdout)
    if p_align.stdout:
        p_align.stdout.close()

    sort_returncode = p_sort.wait()
    align_returncode = p_align.wait()
    if align_returncode != 0:
        raise subprocess.CalledProcessError(align_returncode, cmd_align)
    if sort_returncode != 0:
        raise subprocess.CalledProcessError(sort_returncode, cmd_sort)

    subprocess.run(cmd_index, check=True)


def align_single_trimmed_read_clusters(
    input_db: Path,
    reference: Path,
    alignments: Path,
    crids: list[int],
    output_dir: Path,
    buffer_clipped_length: int = 0,
    threads: int = 1,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)

    containers = consensus.load_crs_containers_from_db(path_db=input_db, crIDs=crids)
    if len(containers) == 0:
        raise ValueError(f"No containers found in {input_db} for crIDs={crids}")

    requested_crids = set(containers.keys())
    consensus_objects = consensus_class.read_consensus_from_db(
        database=input_db,
        crIDs=requested_crids,
    )
    selected_consensus_objects = [
        co
        for co in consensus_objects
        if len(requested_crids.intersection(set(co.crIDs))) > 0
    ]

    if len(selected_consensus_objects) == 0:
        raise ValueError(
            f"No consensus objects found in {input_db} for requested crIDs={sorted(requested_crids)}"
        )

    written_bams: list[Path] = []
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        for consensus_object in selected_consensus_objects:
            trimmed_reads = _trim_reads_for_consensus(
                consensus_object=consensus_object,
                alignments_path=alignments,
                reference=reference,
                buffer_clipped_length=buffer_clipped_length,
            )
            if len(trimmed_reads) == 0:
                log.warning(
                    f"No trimmed reads available for consensus {consensus_object.ID}; skipping."
                )
                continue

            safe_consensus_id = _sanitize_consensus_id(consensus_object.ID)
            tmp_fastq = tmp_dir_path / f"{safe_consensus_id}.trimmed.fastq"
            _write_trimmed_reads_fastq(trimmed_reads=trimmed_reads, output_fastq=tmp_fastq)

            output_bam = output_dir / f"{safe_consensus_id}.trimmed.bam"
            log.info(
                f"Aligning {len(trimmed_reads)} trimmed reads of consensus {consensus_object.ID} to reference."
            )
            _align_fastq_to_reference(
                reference=reference,
                reads_fastq=tmp_fastq,
                output_bam=output_bam,
                threads=threads,
            )
            written_bams.append(output_bam)

    return written_bams


def run(args, **kwargs):
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        force=True,
    )

    written_bams = align_single_trimmed_read_clusters(
        input_db=args.input,
        reference=args.reference,
        alignments=args.alignments,
        crids=args.crids,
        output_dir=args.output_dir,
        buffer_clipped_length=args.buffer_clipped_length,
        threads=args.threads,
    )
    log.info(f"Finished writing {len(written_bams)} BAM files to {args.output_dir}")


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to crs_containers.db file.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to reference FASTA.",
    )
    parser.add_argument(
        "-a",
        "--alignments",
        type=Path,
        required=True,
        help="Path to original read alignments (BAM or CRAM).",
    )
    parser.add_argument(
        "--crids",
        type=int,
        nargs="+",
        required=True,
        help="Candidate region IDs to process (at least one).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory where one BAM per consensus is written and indexed.",
    )
    parser.add_argument(
        "--buffer-clipped-length",
        type=int,
        default=0,
        help="Maximum clipped sequence buffer used when trimming reads (default: 0).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for minimap2 alignment (default: 1).",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level.",
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Trim reads for selected consensus clusters and align each cluster to reference."
    )
    add_arguments(parser)
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
