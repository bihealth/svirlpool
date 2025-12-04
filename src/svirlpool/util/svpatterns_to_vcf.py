#!/usr/bin/env python3
"""
Generate VCF debug output from svpatterns.db file.

This module reads SVpatterns from a database file and outputs them as VCF records,
representing INS, DEL, and BND variants from the underlying SVprimitives.
"""

import argparse
import logging
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from shlex import split

import tqdm

from ..localassembly import SVpatterns
from ..util import datatypes

log = logging.getLogger(__name__)


def generate_vcf_header(reference: Path, samplename: str) -> list[str]:
    """Generate VCF header lines."""
    header_lines = [
        "##fileformat=VCFv4.2",
        f"##fileDate={datetime.now().strftime('%Y%m%d')}",
        "##source=svirlpool_svpatterns_to_vcf",
        f"##reference={reference}",
    ]

    # Add contig lines from reference index
    fai_path = Path(str(reference) + ".fai")
    if fai_path.exists():
        with open(fai_path, "r") as fai:
            for line in fai:
                fields = line.strip().split("\t")
                if len(fields) >= 2:
                    contig_name = fields[0]
                    contig_length = fields[1]
                    header_lines.append(
                        f"##contig=<ID={contig_name},length={contig_length}>"
                    )

    # Add INFO fields
    header_lines.extend([
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
        '##INFO=<ID=CONSENSUS_ID,Number=1,Type=String,Description="Consensus ID from which this variant originated">',
        '##INFO=<ID=ALIGNMENT_ID,Number=1,Type=Integer,Description="Alignment ID on the consensus">',
        '##INFO=<ID=SV_ID,Number=1,Type=Integer,Description="SV ID on the consensus alignment">',
        '##INFO=<ID=SVPATTERN_TYPE,Number=1,Type=String,Description="Type of SVpattern (INS/DEL/INV/etc)">',
        '##INFO=<ID=REPEAT_IDS,Number=.,Type=Integer,Description="Tandem repeat IDs overlapping this variant">',
        '##INFO=<ID=IS_REVERSE,Number=0,Type=Flag,Description="Alignment is on reverse strand">',
        '##INFO=<ID=SEQ_COMPLEXITY,Number=1,Type=Float,Description="Mean sequence complexity score">',
        '##ALT=<ID=INS,Description="Insertion">',
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=INV,Description="Inversion">',
        '##ALT=<ID=DUP,Description="Duplication">',
    ])

    # Add FORMAT fields
    header_lines.extend([
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    ])

    # Add column header line
    header_lines.append(
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samplename}"
    )

    return header_lines


def get_ref_base(reference: Path, chrom: str, pos: int) -> str:
    """Get reference base at position using samtools faidx."""
    try:
        cmd = ["samtools", "faidx", str(reference), f"{chrom}:{pos}-{pos}"]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split("\n")
        if len(lines) > 1:
            return lines[1].upper()
        return "N"
    except Exception as e:
        log.warning(f"Failed to get reference base at {chrom}:{pos}: {e}")
        return "N"


def get_ref_sequence(reference: Path, chrom: str, start: int, end: int) -> str:
    """Get reference sequence using samtools faidx."""
    try:
        cmd = ["samtools", "faidx", str(reference), f"{chrom}:{start}-{end}"]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split("\n")
        if len(lines) > 1:
            return "".join(lines[1:]).upper()
        return "N"
    except Exception as e:
        log.warning(f"Failed to get reference sequence at {chrom}:{start}-{end}: {e}")
        return "N"


def svprimitive_to_vcf_record(
    svp: SVpatterns.SVprimitive,
    reference: Path,
    samplename: str,
    svpattern_type: str,
) -> str | None:
    """Convert an SVprimitive to a VCF record line.

    SVprimitive sv_type values:
    0: Insertion
    1: Deletion
    2: Deletion
    3: Breakend (left)
    4: Breakend (right)
    """

    # Map sv_type to string representation
    sv_type_str = datatypes.SV_TYPE_DICT.get(svp.sv_type, "UNK")

    chrom = svp.chr

    # Build INFO field
    info_parts = [
        f"SVTYPE={sv_type_str}",
        f"SVLEN={svp.size if svp.sv_type == 0 else -svp.size}",
        f"END={svp.ref_end}",
        f"CONSENSUS_ID={svp.consensusID}",
        f"ALIGNMENT_ID={svp.alignmentID}",
        f"SV_ID={svp.svID}",
        f"SVPATTERN_TYPE={svpattern_type}",
    ]

    if svp.repeatIDs:
        repeat_ids_str = ",".join(map(str, sorted(svp.repeatIDs)))
        info_parts.append(f"REPEAT_IDS={repeat_ids_str}")

    if svp.aln_is_reverse:
        info_parts.append("IS_REVERSE")

    info = ";".join(info_parts)

    # Build VCF record based on SV type
    if svp.sv_type == 0:  # Insertion
        # For insertions, POS is the position before the insertion
        pos = svp.ref_start
        ref = get_ref_base(reference, chrom, pos)

        # Get alt sequence
        alt_seq = ""
        if svp.original_alt_sequences and len(svp.original_alt_sequences) > 0:
            alt_seq = "".join(svp.original_alt_sequences)

        if alt_seq:
            alt = ref + alt_seq
        else:
            # Symbolic allele
            alt = "<INS>"

        vcf_id = svp.get_vcfID()
        qual = "."
        filter_field = "PASS"
        format_field = "GT:DP"
        sample_field = "./.:."

        return f"{chrom}\t{pos}\t{vcf_id}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}\t{format_field}\t{sample_field}"

    elif svp.sv_type in [1, 2]:  # Deletion (both type 1 and 2)
        # For deletions, POS is the position before the deletion
        pos = svp.ref_start
        ref_seq = get_ref_sequence(reference, chrom, pos, svp.ref_end)
        ref = ref_seq if ref_seq else "N"
        alt = ref[0] if len(ref) > 0 else "N"

        vcf_id = svp.get_vcfID()
        qual = "."
        filter_field = "PASS"
        format_field = "GT:DP"
        sample_field = "./.:."

        return f"{chrom}\t{pos}\t{vcf_id}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}\t{format_field}\t{sample_field}"

    elif svp.sv_type in [3, 4]:  # BND (breakends - left and right)
        pos = svp.ref_start
        ref = get_ref_base(reference, chrom, pos)

        # For BNDs, we need to represent the breakend
        # Format: t[p[ or ]p]t depending on orientation
        if svp.adjacent_bnd:
            alt = svp.adjacent_bnd.tp_str
        else:
            # Generic BND without mate info
            alt = "<BND>"

        vcf_id = svp.get_vcfID()
        qual = "."
        filter_field = "PASS"
        format_field = "GT:DP"
        sample_field = "./.:."

        return f"{chrom}\t{pos}\t{vcf_id}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}\t{format_field}\t{sample_field}"

    else:
        # Should not reach here, but handle unexpected types
        log.warning(f"Unexpected sv_type {svp.sv_type} for SVprimitive")
        return None


def svpatterns_to_vcf(
    input_db: Path,
    output_vcf: Path,
    reference: Path,
    samplename: str | None = None,
) -> None:
    """
    Convert SVpatterns from database to VCF format.

    Args:
        input_db: Path to svpatterns.db file
        output_vcf: Path to output VCF file (can be .vcf or .vcf.gz)
        reference: Path to reference genome FASTA file
        samplename: Sample name to use in VCF (inferred from db if not provided)
    """

    log.info(f"Reading SVpatterns from {input_db}")
    svpatterns_list = SVpatterns.read_svPatterns_from_db(input_db)

    if not svpatterns_list:
        log.warning("No SVpatterns found in database")
        return

    # Infer samplename if not provided
    if samplename is None:
        samplename = svpatterns_list[0].samplename
        log.info(f"Using samplename from database: {samplename}")

    # Determine if we need to compress
    compress = str(output_vcf).endswith(".gz")

    # Write to temporary file if compression is needed
    if compress:
        temp_vcf = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
        output_file = temp_vcf
        temp_path = Path(temp_vcf.name)
    else:
        output_file = open(output_vcf, "w")
        temp_path = None

    try:
        # Write header
        header_lines = generate_vcf_header(reference, samplename)
        for line in header_lines:
            output_file.write(line + "\n")

        # Process each SVpattern
        variant_count = 0
        for svpattern in tqdm.tqdm(svpatterns_list, desc="Writing VCF records"):
            svpattern_type = svpattern.get_sv_type()

            # Write VCF records for each SVprimitive
            for svp in svpattern.SVprimitives:
                vcf_record = svprimitive_to_vcf_record(
                    svp, reference, samplename, svpattern_type
                )
                if vcf_record:
                    output_file.write(vcf_record + "\n")
                    variant_count += 1

        output_file.close()

        log.info(
            f"Wrote {variant_count} variant records from {len(svpatterns_list)} SVpatterns"
        )

        # Compress and index if needed
        if compress and temp_path is not None:
            # Sort the VCF file before compression
            log.info("Sorting VCF file")
            temp_sorted = tempfile.NamedTemporaryFile(
                mode="w", suffix=".vcf", delete=False
            )
            temp_sorted_path = Path(temp_sorted.name)
            temp_sorted.close()

            cmd_sort = f"sort -k1,1 -k2,2n {str(temp_path)} -o {str(temp_sorted_path)}"
            subprocess.check_call(split(cmd_sort))

            # Remove unsorted temp file
            temp_path.unlink()

            log.info(f"Compressing VCF to {output_vcf}")
            with open(output_vcf, "wb") as f:
                cmd_zip = f"bgzip -c {str(temp_sorted_path)}"
                subprocess.check_call(split(cmd_zip), stdout=f)

            # Remove sorted temp file
            temp_sorted_path.unlink()

            # Index with tabix
            log.info("Indexing VCF with tabix")
            try:
                cmd_index = f"tabix -f -p vcf {str(output_vcf)}"
                subprocess.check_call(split(cmd_index))
                log.info(f"Created index: {output_vcf}.tbi")
            except subprocess.CalledProcessError as e:
                log.warning(f"Failed to create tabix index: {e}")

        log.info(f"VCF output written to {output_vcf}")

    except Exception as e:
        log.error(f"Error writing VCF: {e}")
        if temp_path and temp_path.exists():
            temp_path.unlink()
        raise


def run(args):
    """Run the svpatterns_to_vcf conversion."""
    svpatterns_to_vcf(
        input_db=args.input,
        output_vcf=args.output,
        reference=args.reference,
        samplename=args.samplename,
    )


def get_parser():
    """Get argument parser for command-line execution."""
    parser = argparse.ArgumentParser(
        description="Convert SVpatterns database to VCF format for debugging. "
        "Outputs VCF records for each SVprimitive (INS, DEL, BND) in the database."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to svpatterns.db SQLite database file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to output VCF file (.vcf or .vcf.gz)",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to reference genome FASTA file (must be indexed with .fai)",
    )
    parser.add_argument(
        "-s",
        "--samplename",
        type=str,
        default=None,
        help="Sample name to use in VCF header (inferred from database if not provided)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    return parser


def main():
    """Main entry point for command-line execution."""
    parser = get_parser()
    args = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    run(args)


if __name__ == "__main__":
    main()
