# %%
# !/usr/bin/env python
import argparse
import csv
import gzip
import json
import logging
import multiprocessing as mp
import os
import pickle
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path, PosixPath

import cattrs
import numpy as np
import numpy.typing as npt
import psutil
from tqdm import tqdm

from ..svcalling.multisample_sv_calling import cohens_d

# %%
from ..util import datatypes, util

# %%

# Initialize module-level logger
logger = logging.getLogger(__name__)


def setup_logging(log_level: str = "INFO") -> None:
    """Setup logging with the specified log level.

    Args:
        log_level: Logging level as string (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    """
    # Convert string to logging level
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")

    # Configure logging
    logging.basicConfig(
        level=numeric_level,
        format="[%(levelname)s %(asctime)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )

    logger.info(f"Logging level set to {log_level.upper()}")


def log_memory_usage(context: str = "") -> None:
    """Log current memory usage of the process."""
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    mem_mb = mem_info.rss / 1024 / 1024
    logger.info(f"Memory usage {context}: {mem_mb:.2f} MB")


def write_to_bedgraphs(
    signals: list[tuple[str, int, int, int, float]],
    values_signal: npt.NDArray[np.float32],
    values_normalized: npt.NDArray[np.float32],
    values_masked: npt.NDArray[np.float32],
    bedgraph: Path,
) -> None:
    bedgraph = Path(bedgraph)
    logger.info("write bedgraph of signals..")
    with open(bedgraph.with_suffix(".signal.bedgraph"), "w") as f:
        trackdef = "track type=bedGraph name=values_signal description=values_signal \
visibility=display_mode color=150,20,255 altColor=20,150,255 \
priority=10 autoScale=on alwaysZero=off gridDefault=on \
maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper \
windowingFunction=mean smoothingWindow=off"
        print(trackdef, file=f)
        for i, signal in tqdm(enumerate(signals), desc="Writing signal bedgraph"):
            val = float(values_signal[i])
            if not (
                val is None
                or val == float("inf")
                or val == float("-inf")
                or np.isnan(val)
            ):
                print(signal[0], signal[1], signal[2], str(val), file=f)
    with open(bedgraph.with_suffix(".normalized.bedgraph"), "w") as f:
        trackdef = (
            "track type=bedGraph name=values_normalized description=values_normalized \
visibility=display_mode altColor=150,20,255 color=20,150,255 \
priority=10 autoScale=on alwaysZero=off gridDefault=on \
maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper \
windowingFunction=mean smoothingWindow=off"
        )
        print(trackdef, file=f)
        for i, signal in tqdm(enumerate(signals), desc="Writing normalized bedgraph"):
            val = float(values_normalized[i])
            if not (
                val is None
                or val == float("inf")
                or val == float("-inf")
                or np.isnan(val)
            ):
                print(signal[0], signal[1], signal[2], str(val), file=f)
    with open(bedgraph.with_suffix(".masked.bedgraph"), "w") as f:
        trackdef = (
            "track type=bedGraph name=values_signal_masked description=values_signal_masked \
visibility=display_mode altColor=150,50,50 color=250,150,150 \
priority=10 autoScale=on alwaysZero=off gridDefault=on \
maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper \
windowingFunction=mean smoothingWindow=off"
        )
        print(trackdef, file=f)
        for i, signal in tqdm(enumerate(signals), desc="Writing masked bedgraph"):
            val = float(values_masked[i])
            if not (
                val is None
                or val == float("inf")
                or val == float("-inf")
                or np.isnan(val)
            ):
                print(
                    signal[0],
                    signal[1],
                    signal[2],
                    str(float(values_masked[i])),
                    file=f,
                )
    logger.info("Finished bedgraphs.")


def yield_unstructured_crs_from_tsv(tsv_file: Path):
    """Yield unstructured CandidateRegion dicts from a TSV file."""
    # Handle both gzipped and non-gzipped files
    open_func = gzip.open if str(tsv_file).endswith(".gz") else open
    mode = "rt" if str(tsv_file).endswith(".gz") else "r"
    with open_func(tsv_file, mode) as f:
        reader = csv.reader(f, delimiter="\t", quotechar='"')
        for row in reader:
            cr_json = row[3]
            cr_dict = json.loads(cr_json)
            yield cr_dict


def create_crs_db(path_db: Path, timeout: float) -> None:
    # create sqlite3 database with primary key crID and pickled json.dumps(CandidateRegion.unstructure()) as vales named candidate_region
    assert timeout > 0.0, f"timeout must be > 0.0. It is {timeout}"
    assert type(path_db) == PosixPath or type(path_db) == str, (
        f"path_db must be a Path or str. It is {type(path_db)}"
    )

    if path_db.exists():
        logger.warning(f"Database {path_db} exists. Overwriting it.")
        path_db.unlink()

    conn = sqlite3.connect(path_db)
    c = conn.cursor()
    c.execute("""CREATE TABLE IF NOT EXISTS candidate_regions
        (crID INTEGER PRIMARY KEY, candidate_region TEXT)""")
    c.execute(f"pragma busy_timeout={str(int(timeout * 1000))}")
    conn.commit()
    c.close()
    conn.close()


def write_crs_to_db(path_db: Path, crs_iterator) -> None:
    """Write CandidateRegion objects from an iterator to database."""
    assert type(path_db) == PosixPath or type(path_db) == str, (
        f"path_db must be a Path or str. It is {type(path_db)}"
    )
    assert path_db.exists(), f"Database {path_db} does not exist"
    assert path_db.is_file(), f"Database {path_db} is not a file"
    conn = sqlite3.connect(path_db)
    c = conn.cursor()

    # Process in batches to avoid memory issues
    batch_size = 1000
    batch = []
    running_crID: int = 0
    for cr_dict in crs_iterator:
        cr = cattrs.structure(cr_dict, datatypes.CandidateRegion)
        cr.crID = running_crID
        running_crID += 1
        batch.append([cr.crID, pickle.dumps(cr.unstructure())])

        if len(batch) >= batch_size:
            c.executemany("INSERT INTO candidate_regions VALUES (?,?)", batch)
            conn.commit()
            batch = []

    # Insert remaining items
    if batch:
        c.executemany("INSERT INTO candidate_regions VALUES (?,?)", batch)
        conn.commit()

    c.close()
    conn.close()


def load_crs_from_db(
    path_db: Path, crIDs: list[int] | None = None
) -> list[datatypes.CandidateRegion]:
    if crIDs:
        assert type(crIDs) == list, "crIDs is not a list"
        assert all(type(crID) == int for crID in crIDs), (
            "crIDs is not a list of integers"
        )
    assert type(path_db) == PosixPath or type(path_db) == str, (
        f"path_db must be a Path or str. It is {type(path_db)}"
    )
    assert path_db.exists(), f"Database {path_db} does not exist"

    conn = sqlite3.connect("file:" + str(path_db) + "?mode=ro", uri=True)
    c = conn.cursor()
    # use executemany
    if crIDs:
        pickled_crs = c.execute(
            "SELECT candidate_region FROM candidate_regions WHERE crID IN ({crIDs})".format(
                crIDs=",".join([str(crID) for crID in crIDs])
            )
        ).fetchall()
    else:
        # read all pcikled crs objects from db, not the crIDs
        pickled_crs = c.execute(
            "SELECT candidate_region FROM candidate_regions"
        ).fetchall()
    # unpickle and return list of structured crs
    crs = [
        cattrs.structure(pickle.loads(row[0]), datatypes.CandidateRegion)
        for row in pickled_crs
    ]
    c.close()
    conn.close()
    return crs


def has_similar_sv_signals(
    this: datatypes.CandidateRegion,
    other: datatypes.CandidateRegion,
    cohensd: float = 0.5,
    min_abs_signal_size: int = 100,
    min_signals_per_read_fraction: float = 0.2,
) -> bool:
    # hard cutoff: both crs need to have at least min_signals_per_read_fraction * read count SVs of at least min_abs_signal_size
    this_readcount = len(this.get_read_names())
    other_readcount = len(other.get_read_names())
    dels_this = [
        s.size
        for s in this.sv_signals
        if abs(s.size) >= min_abs_signal_size and s.sv_type == 1
    ]
    dels_other = [
        s.size
        for s in other.sv_signals
        if abs(s.size) >= min_abs_signal_size and s.sv_type == 1
    ]
    dels_insufficient = (
        len(dels_this) < min_signals_per_read_fraction * this_readcount
        or len(dels_other) < min_signals_per_read_fraction * other_readcount
    )
    ins_this = [
        s.size
        for s in this.sv_signals
        if abs(s.size) >= min_abs_signal_size and s.sv_type == 0
    ]
    ins_other = [
        s.size
        for s in other.sv_signals
        if abs(s.size) >= min_abs_signal_size and s.sv_type == 0
    ]
    ins_insufficient = (
        len(ins_this) < min_signals_per_read_fraction * this_readcount
        or len(ins_other) < min_signals_per_read_fraction * other_readcount
    )
    if dels_insufficient and ins_insufficient:
        return False
    # both signal sizes distributions of signals greater than min_abs_signal_size of are checkd with cohen's d
    D_dels = (
        cohens_d(np.array(dels_this), np.array(dels_other))
        if len(dels_this) > 0 and len(dels_other) > 0
        else cohensd
    )
    D_ins = (
        cohens_d(np.array(ins_this), np.array(ins_other))
        if len(ins_this) > 0 and len(ins_other) > 0
        else cohensd
    )
    if abs(D_dels) < cohensd or abs(D_ins) < cohensd:
        return True
    return False


def ensure_min_cr_size(
    cr: datatypes.CandidateRegion, min_cr_size: int
) -> datatypes.CandidateRegion:
    size_cr = cr.referenceEnd - cr.referenceStart
    if size_cr < min_cr_size:
        # extend equally on both sides to reach min_cr_size
        difference = min_cr_size - size_cr
        cr.referenceStart = max(0, cr.referenceStart - difference // 2)
        cr.referenceEnd = cr.referenceEnd + (difference - difference // 2)
    return cr


def split_signals_by_chromosome(
    signalstrengths: Path,
    tmp_dir: Path,
    filter_absolute: float,
    filter_normalized: float,
) -> dict[str, Path]:
    """Split signals into chromosome-specific files and filter them.

    Returns dict mapping chromosome -> path to filtered signals file.
    """
    logger.info("Phase 1: Splitting signals by chromosome and filtering...")

    # First pass: collect all signals and filter
    signals_by_chr: dict[str, list[datatypes.ExtendedSVsignal]] = {}

    for svsignal in util.yield_from_extendedSVsignal(
        input=signalstrengths, description="Loading and filtering signals"
    ):
        # Apply filters
        signal_normalized = (
            svsignal.strength / svsignal.coverage if svsignal.coverage > 0 else 0.0
        )
        if (
            svsignal.strength > filter_absolute
            and signal_normalized > filter_normalized
        ):
            if svsignal.chr not in signals_by_chr:
                signals_by_chr[svsignal.chr] = []
            signals_by_chr[svsignal.chr].append(svsignal)

    # Write chromosome-specific files
    chr_signal_files = {}
    for chr_name, signals in signals_by_chr.items():
        chr_file = tmp_dir / f"{chr_name}_signals.tsv"
        with open(chr_file, "w") as f:
            writer = csv.writer(f, delimiter="\t", quotechar='"')
            for signal in signals:
                writer.writerow([
                    signal.chr,
                    signal.ref_start,
                    signal.ref_end,
                    signal.repeatID,
                    signal.coverage,
                    signal.strength,
                    json.dumps(signal.unstructure()),
                ])
        chr_signal_files[chr_name] = chr_file
        logger.info(f"  {chr_name}: {len(signals)} filtered signals")

    return chr_signal_files


def split_tandem_repeats_by_chromosome(
    tandem_repeats: Path, tmp_dir: Path, chromosomes: set[str]
) -> dict[str, Path]:
    """Split tandem repeats into chromosome-specific files.

    Returns dict mapping chromosome -> path to repeats file.
    """
    logger.info("Splitting tandem repeats by chromosome...")

    chr_repeat_files = {}
    chr_file_handles = {}

    # Open file handles for each chromosome
    for chr_name in chromosomes:
        chr_file = tmp_dir / f"{chr_name}_repeats.bed"
        chr_repeat_files[chr_name] = chr_file
        chr_file_handles[chr_name] = open(chr_file, "w")

    # Read and distribute repeats (handle both gzipped and non-gzipped files)
    open_func = gzip.open if str(tandem_repeats).endswith(".gz") else open
    mode = "rt" if str(tandem_repeats).endswith(".gz") else "r"
    with open_func(tandem_repeats, mode) as trf:
        reader = csv.reader(trf, delimiter="\t", quotechar='"')
        for row in reader:
            chr_name = row[0]
            if chr_name in chr_file_handles:
                print("\t".join(row), file=chr_file_handles[chr_name])

    # Close all handles
    for fh in chr_file_handles.values():
        fh.close()

    return chr_repeat_files


def process_chromosome_to_proto_crs(args_tuple) -> tuple[str, Path, dict]:
    """Worker function to process one chromosome and create proto-CRs.

    Args:
        args_tuple: (chr_name, chr_signals_file, chr_repeats_file, tmp_dir, buffer_region_radius)

    Returns:
        (chr_name, proto_crs_file, statistics_dict)
    """
    chr_name, chr_signals_file, chr_repeats_file, tmp_dir, buffer_region_radius = (
        args_tuple
    )

    logger.info(f"Processing chromosome {chr_name}...")

    # Load signals for this chromosome
    signals = []
    with open(chr_signals_file, "r") as f:
        reader = csv.reader(f, delimiter="\t", quotechar='"')
        for i, row in enumerate(reader):
            signal_dict = json.loads(row[6])
            signal = cattrs.structure(signal_dict, datatypes.ExtendedSVsignal)
            signals.append((signal.chr, signal.ref_start, signal.ref_end, i, signal))

    if not signals:
        logger.warning(f"No signals found for chromosome {chr_name}")
        # Return empty results
        proto_crs_file = tmp_dir / f"{chr_name}_proto_crs.tsv"
        proto_crs_file.touch()
        return (
            chr_name,
            proto_crs_file,
            {"chr": chr_name, "n_proto_crs": 0, "read_counts": []},
        )

    # Write signals to BED file with margins
    tmp_signals_bed = tmp_dir / f"{chr_name}_signals.bed"
    with open(tmp_signals_bed, "w") as f:
        for s in signals:
            start = max(0, s[1] - 50)
            end = s[2] + 50
            index = s[3]
            print(s[0], start, end, index, sep="\t", file=f)

    # Add tandem repeats to BED file
    if chr_repeats_file.exists():
        with open(chr_repeats_file, "r") as trf, open(tmp_signals_bed, "a") as f:
            for line in trf:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    print(parts[0], parts[1], parts[2], -1, sep="\t", file=f)

    # Sort the BED file (required for bedtools merge)
    sorted_bed = tmp_dir / f"{chr_name}_signals_sorted.bed"
    cmd_sort = ["sort", "-k1,1", "-k2,2n", str(tmp_signals_bed)]
    with open(sorted_bed, "w") as sf:
        subprocess.check_call(cmd_sort, stdout=sf, stderr=subprocess.DEVNULL)

    # Run bedtools merge
    merged_file = tmp_dir / f"{chr_name}_merged.bed"
    cmd = [
        "bedtools",
        "merge",
        "-d",
        str(buffer_region_radius),
        "-i",
        str(sorted_bed),
        "-c",
        "4",
        "-o",
        "collapse",
    ]
    with open(merged_file, "w") as mf:
        subprocess.check_call(cmd, stdout=mf, stderr=subprocess.DEVNULL)

    # Create proto-CRs from merged regions
    proto_crs_file = tmp_dir / f"{chr_name}_proto_crs.tsv"
    read_counts = []

    with open(merged_file, "r") as mf, open(proto_crs_file, "w") as pcf:
        reader = csv.reader(mf, delimiter="\t")
        writer = csv.writer(pcf, delimiter="\t", quotechar='"')

        for row in tqdm(reader, desc=f"Creating proto-CRs for {chr_name}"):
            start = int(row[1])
            end = int(row[2])
            indices_str = row[3]
            indices = [int(idx) for idx in indices_str.split(",") if int(idx) != -1]

            if not indices:  # Skip if only tandem repeats
                continue

            # Collect signals for this proto-CR
            region_signals = [signals[idx][4] for idx in indices]

            # Create candidate region
            cr = datatypes.CandidateRegion.from_sv_signals(
                sv_signals=region_signals,
                crID=0,  # Will be reassigned later
            )

            read_counts.append(len(cr.get_read_names()))
            writer.writerow([
                cr.chr,
                cr.referenceStart,
                cr.referenceEnd,
                json.dumps(cr.unstructure()),
            ])

    # Compute statistics
    stats = {
        "chr": chr_name,
        "n_proto_crs": len(read_counts),
        "read_counts": read_counts,
    }

    stats_file = tmp_dir / f"{chr_name}_stats.json"
    with open(stats_file, "w") as sf:
        json.dump(stats, sf)

    logger.info(f"Chromosome {chr_name}: created {len(read_counts)} proto-CRs")

    return chr_name, proto_crs_file, stats


def compute_global_statistics(chr_stats: list[dict]) -> dict:
    """Compute global statistics from all chromosome statistics.

    Returns dict with global statistics including median read count.
    """
    logger.info("Phase 3: Computing global statistics...")

    all_read_counts = []
    total_proto_crs = 0

    for stats in chr_stats:
        all_read_counts.extend(stats["read_counts"])
        total_proto_crs += stats["n_proto_crs"]

    median_read_count = float(np.median(all_read_counts)) if all_read_counts else 0.0

    global_stats = {
        "total_proto_crs": total_proto_crs,
        "median_read_count": median_read_count,
        "total_reads_across_all_crs": sum(all_read_counts),
    }

    logger.info(f"Total proto-CRs: {total_proto_crs}")
    logger.info(f"Global median read count: {median_read_count:.1f}")

    return global_stats


def filter_and_merge_chromosome(args_tuple) -> tuple[str, Path, Path, dict]:
    """Worker function to filter and merge CRs for one chromosome.

    Args:
        args_tuple: (chr_name, proto_crs_file, tmp_dir, median_read_count,
                     cutoff_multiplier, min_cr_size)

    Returns:
        (chr_name, final_crs_file, dropped_crs_file, stats_dict)
    """
    (
        chr_name,
        proto_crs_file,
        tmp_dir,
        median_read_count,
        cutoff_multiplier,
        min_cr_size,
    ) = args_tuple

    logger.info(f"Filtering and merging chromosome {chr_name}...")

    max_read_count = cutoff_multiplier * median_read_count

    final_crs_file = tmp_dir / f"{chr_name}_final.tsv"
    dropped_crs_file = tmp_dir / f"{chr_name}_dropped.tsv"

    n_filtered = 0
    n_kept = 0
    n_merged = 0

    # Check if proto_crs file is empty
    if proto_crs_file.stat().st_size == 0:
        logger.info(f"No proto-CRs for chromosome {chr_name}, skipping")
        final_crs_file.touch()
        dropped_crs_file.touch()
        return (
            chr_name,
            final_crs_file,
            dropped_crs_file,
            {"chr": chr_name, "n_kept": 0, "n_filtered": 0, "n_merged": 0},
        )

    # Handle both gzipped and non-gzipped files
    open_func = gzip.open if str(proto_crs_file).endswith(".gz") else open
    mode = "rt" if str(proto_crs_file).endswith(".gz") else "r"

    with (
        open_func(proto_crs_file, mode) as pcf,
        open(final_crs_file, "w") as fcf,
        open(dropped_crs_file, "w") as dcf,
    ):
        reader = csv.reader(pcf, delimiter="\t", quotechar='"')
        writer_final = csv.writer(fcf, delimiter="\t", quotechar='"')
        writer_dropped = csv.writer(dcf, delimiter="\t", quotechar='"')

        previous_cr: datatypes.CandidateRegion | None = None

        for row in reader:
            cr_dict = json.loads(row[3])
            cr = cattrs.structure(cr_dict, datatypes.CandidateRegion)

            # Filter by read count
            read_count = len(cr.get_read_names())
            if read_count > max_read_count:
                writer_dropped.writerow(row)
                n_filtered += 1
                continue

            # Try to merge with previous CR
            if previous_cr is not None:
                # Check if should merge (same chromosome, close proximity, similar signals)
                if (
                    previous_cr.chr == cr.chr
                    and cr.referenceStart <= previous_cr.referenceEnd + 100_000
                    and has_similar_sv_signals(previous_cr, cr)
                ):
                    # Merge
                    merged_signals = previous_cr.sv_signals + cr.sv_signals
                    merged_start = min(previous_cr.referenceStart, cr.referenceStart)
                    merged_end = max(previous_cr.referenceEnd, cr.referenceEnd)
                    previous_cr = datatypes.CandidateRegion.from_sv_signals(
                        crID=previous_cr.crID, sv_signals=merged_signals
                    )
                    n_merged += 1
                    logger.debug(
                        f"Merged CRs at {previous_cr.chr}:{merged_start}-{merged_end}"
                    )
                    continue
                else:
                    # Write previous CR
                    previous_cr = ensure_min_cr_size(previous_cr, min_cr_size)
                    writer_final.writerow([
                        previous_cr.chr,
                        previous_cr.referenceStart,
                        previous_cr.referenceEnd,
                        json.dumps(previous_cr.unstructure()),
                    ])
                    n_kept += 1
                    logger.debug(
                        f"Wrote CR at {previous_cr.chr}:{previous_cr.referenceStart}-{previous_cr.referenceEnd}"
                    )

            previous_cr = cr

        # Write last CR
        if previous_cr is not None:
            previous_cr = ensure_min_cr_size(previous_cr, min_cr_size)
            writer_final.writerow([
                previous_cr.chr,
                previous_cr.referenceStart,
                previous_cr.referenceEnd,
                json.dumps(previous_cr.unstructure()),
            ])
            n_kept += 1
            logger.debug(
                f"Wrote CR at {previous_cr.chr}:{previous_cr.referenceStart}-{previous_cr.referenceEnd}"
            )
    stats = {
        "chr": chr_name,
        "n_kept": n_kept,
        "n_filtered": n_filtered,
        "n_merged": n_merged,
    }

    logger.info(
        f"Chromosome {chr_name}: kept={n_kept}, filtered={n_filtered}, merged={n_merged}"
    )

    return chr_name, final_crs_file, dropped_crs_file, stats


def concatenate_and_write_to_db(
    chr_final_files: list[Path],
    chr_dropped_files: list[Path],
    output_db: Path,
    dropped_db: Path | None,
    chromosomes_order: list[str],
) -> None:
    """Phase 5: Concatenate chromosome results and write to databases."""
    logger.info("Phase 5: Concatenating results and writing to databases...")

    # Create mapping from chromosome to files
    chr_to_final = {f.stem.replace("_final", ""): f for f in chr_final_files}
    chr_to_dropped = {f.stem.replace("_dropped", ""): f for f in chr_dropped_files}

    # Create temporary concatenated files
    tmp_all_final = output_db.parent / "tmp_all_final.tsv"
    tmp_all_dropped = output_db.parent / "tmp_all_dropped.tsv" if dropped_db else None

    # Concatenate final CRs in chromosome order
    n_final_crs = 0
    with open(tmp_all_final, "w") as outf:
        writer = csv.writer(outf, delimiter="\t", quotechar='"')
        for chr_name in chromosomes_order:
            # Try to find the file for this chromosome
            chr_file = chr_to_final.get(chr_name)
            if chr_file and chr_file.exists() and chr_file.stat().st_size > 0:
                with open(chr_file, "r") as inf:
                    reader = csv.reader(inf, delimiter="\t", quotechar='"')
                    for row in reader:
                        writer.writerow(row)
                        n_final_crs += 1

    logger.info(f"Concatenated {n_final_crs} final CRs from all chromosomes")

    # Write to output database
    if n_final_crs > 0:
        logger.info("Writing final CRs to output database...")
        create_crs_db(path_db=output_db, timeout=240.0)
        write_crs_to_db(
            path_db=output_db,
            crs_iterator=yield_unstructured_crs_from_tsv(tmp_all_final),
        )
    else:
        logger.warning("No final CRs to write to database")

    # Handle dropped CRs if requested
    if dropped_db is not None and tmp_all_dropped is not None:
        n_dropped_crs = 0
        with open(tmp_all_dropped, "w") as outf:
            writer = csv.writer(outf, delimiter="\t", quotechar='"')
            for chr_name in chromosomes_order:
                chr_file = chr_to_dropped.get(chr_name)
                if chr_file and chr_file.exists() and chr_file.stat().st_size > 0:
                    with open(chr_file, "r") as inf:
                        reader = csv.reader(inf, delimiter="\t", quotechar='"')
                        for row in reader:
                            writer.writerow(row)
                            n_dropped_crs += 1

        if n_dropped_crs > 0:
            logger.info(f"Writing {n_dropped_crs} dropped CRs to database...")
            create_crs_db(path_db=dropped_db, timeout=240.0)
            write_crs_to_db(
                path_db=dropped_db,
                crs_iterator=yield_unstructured_crs_from_tsv(tmp_all_dropped),
            )
        else:
            logger.info("No dropped CRs to write")

    # Cleanup temporary files
    tmp_all_final.unlink()
    if tmp_all_dropped and tmp_all_dropped.exists():
        tmp_all_dropped.unlink()

    logger.info("Database writing complete.")


def create_candidate_regions(
    reference: Path,
    signalstrengths: Path,
    output: Path,
    tandem_repeats: Path,
    threads: int,
    buffer_region_radius: int,
    min_cr_size: int,
    filter_absolute: float,
    filter_normalized: float,
    cutoff_median_readcount_per_region: float = 6.0,
    dropped: Path | None = None,
    bedgraph: Path | None = None,
    tmp_dir_path: Path | None = None,
) -> None:
    csv.field_size_limit(sys.maxsize)

    with tempfile.TemporaryDirectory(dir=tmp_dir_path) as tdir:
        tmp_dir = Path(tdir)

        # Read chromosome names from fasta.fai
        fai_file = Path(str(reference) + ".fai")
        chromosomes: dict[str, int] = {}
        with open(fai_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                chromosomes[parts[0]] = int(parts[1])
        chromosomes_order = list(chromosomes.keys())

        # PHASE 1: Split inputs by chromosome
        chr_signal_files = split_signals_by_chromosome(
            signalstrengths=signalstrengths,
            tmp_dir=tmp_dir,
            filter_absolute=filter_absolute,
            filter_normalized=filter_normalized,
        )

        chr_repeat_files = split_tandem_repeats_by_chromosome(
            tandem_repeats=tandem_repeats,
            tmp_dir=tmp_dir,
            chromosomes=set(chr_signal_files.keys()),
        )

        # Optional: Generate bedgraph from filtered signals in parallel
        bedgraph_process = None
        if bedgraph is not None:
            logger.info("Starting bedgraph generation in background...")
            # Collect all filtered signals for bedgraph
            all_signals = []
            all_strengths = []
            all_depths = []

            for chr_name in chromosomes_order:
                if chr_name not in chr_signal_files:
                    continue
                with open(chr_signal_files[chr_name], "r") as f:
                    reader = csv.reader(f, delimiter="\t", quotechar='"')
                    for row in reader:
                        chr_name_sig = row[0]
                        start = int(row[1])
                        end = int(row[2])
                        coverage = int(row[4])
                        strength = float(row[5])
                        all_signals.append((
                            chr_name_sig,
                            start,
                            end,
                            coverage,
                            strength,
                        ))
                        all_strengths.append(strength)
                        all_depths.append(coverage)

            values_strength = np.array(all_strengths, dtype=float)
            values_depth = np.array(all_depths, dtype=int)

            with np.errstate(divide="ignore", invalid="ignore"):
                signal_normalized = np.divide(values_strength, values_depth)
                signal_normalized[~np.isfinite(signal_normalized)] = 0.0

            mask_final = np.ones(len(all_signals), dtype=bool)

            # Start bedgraph writing in a separate process
            bedgraph_process = mp.Process(
                target=write_to_bedgraphs,
                args=(
                    all_signals,
                    values_strength,
                    signal_normalized,
                    mask_final.astype(float),
                    bedgraph,
                ),
            )
            bedgraph_process.start()
            logger.info("Bedgraph generation running in background process")

        # PHASE 2: Process each chromosome in parallel
        logger.info(
            f"Phase 2: Processing {len(chr_signal_files)} chromosomes in parallel..."
        )

        process_args = []
        for chr_name in chromosomes_order:
            if chr_name not in chr_signal_files:
                continue
            process_args.append((
                chr_name,
                chr_signal_files[chr_name],
                chr_repeat_files.get(chr_name, tmp_dir / f"{chr_name}_repeats.bed"),
                tmp_dir,
                buffer_region_radius,
            ))

        chr_results = []
        with mp.Pool(processes=threads) as pool:
            chr_results = pool.map(process_chromosome_to_proto_crs, process_args)

        logger.info("Chromosome processing complete.")

        # Extract statistics
        chr_stats = [result[2] for result in chr_results]
        chr_proto_files = {result[0]: result[1] for result in chr_results}

        # PHASE 3: Compute global statistics
        logger.info("Phase 3: Computing global statistics...")
        global_stats = compute_global_statistics(chr_stats)
        median_read_count = global_stats["median_read_count"]

        # PHASE 4: Filter and merge each chromosome in parallel
        logger.info(
            f"Phase 4: Filtering and merging {len(chr_proto_files)} chromosomes in parallel..."
        )

        filter_args = []
        for chr_name in chromosomes_order:
            if chr_name not in chr_proto_files:
                continue
            filter_args.append((
                chr_name,
                chr_proto_files[chr_name],
                tmp_dir,
                median_read_count,
                cutoff_median_readcount_per_region,
                min_cr_size,
            ))

        logger.info("Filtering and merging candidate regions...")
        filter_results = []
        with mp.Pool(processes=threads) as pool:
            filter_results = list(
                tqdm(
                    pool.imap(filter_and_merge_chromosome, filter_args),
                    total=len(filter_args),
                    desc="Filtering and merging chromosomes",
                )
            )

        # Extract results
        chr_final_files = [result[1] for result in filter_results]
        chr_dropped_files = [result[2] for result in filter_results]
        filter_stats = [result[3] for result in filter_results]

        # Log summary statistics
        total_kept = sum(s["n_kept"] for s in filter_stats)
        total_filtered = sum(s["n_filtered"] for s in filter_stats)
        total_merged = sum(s["n_merged"] for s in filter_stats)

        logger.info(
            f"Summary: kept={total_kept}, filtered={total_filtered}, merged={total_merged}"
        )

        if total_kept == 0:
            raise ValueError("No candidate regions remaining after filtering. Exiting.")

        # PHASE 5: Concatenate and write to databases
        concatenate_and_write_to_db(
            chr_final_files=chr_final_files,
            chr_dropped_files=chr_dropped_files,
            output_db=output,
            dropped_db=dropped,
            chromosomes_order=chromosomes_order,
        )

        logger.info(f"Done. Created and stored {total_kept} candidate regions.")
        log_memory_usage("final")

        # Wait for bedgraph process to complete if it's running
        if bedgraph_process is not None:
            logger.info("Waiting for bedgraph generation to complete...")
            bedgraph_process.join()
            logger.info("Bedgraph generation completed.")


def run(args, **kwargs):
    # Setup logging with the specified level
    setup_logging(args.log_level)

    create_candidate_regions(
        signalstrengths=args.input,
        reference=args.reference,
        output=args.output,
        tandem_repeats=args.tandem_repeats,
        buffer_region_radius=args.buffer_region_radius,
        bedgraph=args.bedgraph,
        filter_absolute=args.filter_absolute,
        filter_normalized=args.filter_normalized,
        cutoff_median_readcount_per_region=args.cutoff_median_readcount_per_region,
        threads=args.threads,
        min_cr_size=args.min_cr_size,
        dropped=args.dropped,
        tmp_dir_path=getattr(args, "tmp_dir", None),
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Writes crs file from s signals, repeats and a referene fasta file."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to signalstrengths file that was signaldepths_to_signalstrength.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to reference in uncompressed .fa(sta) format. Only the fasta index (.fai) is required to exist.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to database with schema [crID,CandidateRegion].",
    )
    parser.add_argument(
        "-u",
        "--tandem-repeats",
        type=Path,
        required=True,
        help="Path to bed file with tandem repeats. Can be generated with 'tandem_repeats'.",
    )
    parser.add_argument(
        "--dropped",
        type=Path,
        required=False,
        default=None,
        help="Path to database with schema [crID,CandidateRegion].",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=8,
        help="Number of threads to use.",
    )
    parser.add_argument(
        "-b",
        "--buffer-region-radius",
        type=int,
        required=False,
        default=300,
        help="Radius of the initial seeds' sizes around signal that can constitute a candidate region seed.",
    )
    parser.add_argument(
        "--cutoff-median-readcount-per-region",
        type=float,
        required=False,
        default=6.0,
        help="Cutoff multiplier for median read count per candidate region to filter out regions with excessive read counts.",
    )
    parser.add_argument(
        "--bedgraph",
        required=False,
        default=None,
        help="Path to bedgraph file to write signal values to.",
    )
    parser.add_argument(
        "--filter-absolute",
        type=float,
        required=False,
        default=1.2,
        help="Filter signal values below this value.",
    )
    parser.add_argument(
        "--filter-normalized",
        type=float,
        required=False,
        default=0.12,
        help="Filter signal values below this value.",
    )
    parser.add_argument(
        "--min-cr-size",
        type=int,
        required=False,
        default=1200,
        help="Minimum size of candidate region.",
    )
    parser.add_argument(
        "--tmp-dir",
        type=Path,
        required=False,
        default=None,
        help="Directory for temporary files.",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        required=False,
        default="DEBUG",
        help="Logging level. One of DEBUG, INFO, WARNING, ERROR, CRITICAL.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
