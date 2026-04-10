# %%


import logging
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from intervaltree import Interval
from pysam import AlignedSegment

from ..localassembly import consensus
from ..util.datatypes import MergedSVSignal, SVsignal
from . import consensus_class

log = logging.getLogger(__name__)

# %%


# %%


def ref_and_read_start_from_merging_sv_signals(
    sv_signals: list[MergedSVSignal],
) -> tuple[int, int]:
    """returns the ref_start and read_start of the first signal in the list. Can be more elaborate in the future"""
    sv_signals = sorted(sv_signals, key=lambda x: x.ref_start)
    merged_ref_start = sv_signals[0].ref_start
    merged_read_start = sv_signals[0].read_start
    return (merged_ref_start, merged_read_start)


def alt_sequence_for_MergedSVSignal(
    merged_signal: MergedSVSignal,
    consensus_sequence: str,
    reverse: bool,
    interval_core: tuple[int, int],
) -> str:
    assert len(consensus_sequence) > 0, "consensus_sequence must have a length > 0"
    # there are different strategies for insertions and deletions. The simplest one is to
    # cut the consensus sequence at the ref_start position and add the merged signal size
    # iterate SV signals and add an alt sequence for each one. '' for each DEL and the sequence for each INS
    seq = ""
    readstart = merged_signal.read_start - interval_core[0]
    readend = merged_signal.read_end - interval_core[0]
    if merged_signal.sv_type == 0:
        alt_seq = Seq(consensus_sequence[readstart:readend])
        seq = str(alt_seq.reverse_complement()) if reverse else str(alt_seq)
    elif merged_signal.sv_type == 3:  # BND
        if reverse:
            seq = str(Seq(consensus_sequence[readstart]).reverse_complement())
        else:
            seq = str(consensus_sequence[readstart])
        # if reverse:
        #     seq = str(Seq(consensus_sequence[readstart:]).reverse_complement())
        # else:
        #     seq = str(consensus_sequence[:readstart])
    elif merged_signal.sv_type == 4:  # BND
        if reverse:
            seq = str(Seq(consensus_sequence[readstart]).reverse_complement())
        else:
            seq = str(consensus_sequence[readstart])
        # if reverse:
        #     seq = str(Seq(consensus_sequence[:readstart]).reverse_complement())
        # else:
        #     seq = str(consensus_sequence[readstart:])
    elif merged_signal.sv_type == 1:  # DEL
        seq = ""
    else:
        raise ValueError(
            f"merged signal type must be 0 (INS), 1 (DEL), 3 (BND) or 4 (BND). Got {merged_signal.sv_type}"
        )
    return seq


# to debug and test merge_svs_in_dict_alignments, we need to save all input to a pickled file
def save_merge_svs_in_dict_alignments_input(input, path: Path | str):
    """Saves the input for merge_svs_in_dict_alignments to a file."""
    if isinstance(path, str):
        path = Path(path)
    with open(path, "wb") as f:
        import pickle

        pickle.dump(input, f)
    log.info(f"Saved input for merge_svs_in_dict_alignments to {path}")


def assign_repeat_ids(
    sv_signals: list[MergedSVSignal],
    trf_intervals: list[tuple[int, int, int]],
) -> None:
    # trf_intervals: list of (start, end, repeatID)
    # Precompute Interval objects for all trf_intervals
    precomputed_trf_intervals = [
        (Interval(start, end), repeatID) for start, end, repeatID in trf_intervals
    ]
    for sv in sv_signals:
        sv_interval = Interval(sv.ref_start, sv.ref_end)
        sv.repeatIDs = []
        for trf_interval, repeatID in precomputed_trf_intervals:
            if sv_interval.overlaps(trf_interval):
                sv.repeatIDs.append(repeatID)
        sv.repeatIDs = list(set(sv.repeatIDs))


def parse_sv_signals_from_consensus(
    samplename: str,
    reference_name: str,
    consensus_sequence: str,
    pysam_alignment: AlignedSegment,
    interval_core: tuple[int, int],
    trf_intervals: list[tuple[int, int, int]],
    min_signal_size: int = 5,
    min_bnd_size: int = 50,
) -> list[MergedSVSignal]:
    """
    Parse SV signals from consensus alignment. Returned items are sorted by their position on the consensus.
    """
    # parse sv signals from alignments
    sv_signals: list[SVsignal] = consensus.parse_ReadAlignmentSignals_from_alignment(
        samplename=samplename,
        alignment=pysam_alignment,
        min_signal_size=min_signal_size,
        min_bnd_size=min_bnd_size,
    ).SV_signals

    if interval_core is not None:
        _filtered_outside_core = [
            sv_signal
            for sv_signal in sv_signals
            if not (
                sv_signal.read_start > interval_core[0]
                and sv_signal.read_end < interval_core[1]
            )
        ]
        sv_signals = [
            sv_signal
            for sv_signal in sv_signals
            if sv_signal.read_start > interval_core[0]
            and sv_signal.read_end < interval_core[1]
        ]
        # also log the signals that were taken
        for signal in sv_signals:
            _consensusID = str(pysam_alignment.query_name)
            _crID = _consensusID.split(".")[0]
            log.debug(
                f"TRANSFORMED::parse_sv_signals_from_consensus:(keeping signals inside of interval_core {interval_core}): consensusID={_consensusID}, crID={_crID}, description={signal._get_description(reference_name)}"
            )

        if _filtered_outside_core:
            _consensusID = str(pysam_alignment.query_name)
            _crID = _consensusID.split(".")[0]
            log.debug(
                f"DROPPED::parse_sv_signals_from_consensus:(dropping signals outside of interval_core {interval_core}): consensusID={_consensusID}, crID={_crID}, description={[sv_signal._get_description(reference_name) for sv_signal in _filtered_outside_core]}"
            )

    # transform all SVsignals to MergedSVSignals
    merged_sv_signals = [
        MergedSVSignal(
            chr=reference_name,
            ref_start=sv_signal.ref_start,
            ref_end=sv_signal.ref_end,
            read_start=sv_signal.read_start,
            read_end=sv_signal.read_end,
            size=sv_signal.size,
            sv_type=sv_signal.sv_type,
            repeatIDs=[],
            original_alt_sequences=[],
            original_ref_sequences=[],
        )
        for sv_signal in sv_signals
    ]

    assign_repeat_ids(merged_sv_signals, trf_intervals)

    # add alt sequences to all MergedSVSignals
    for merged_signal in merged_sv_signals:
        merged_signal.original_alt_sequences = [
            alt_sequence_for_MergedSVSignal(
                merged_signal=merged_signal,
                consensus_sequence=consensus_sequence,
                reverse=pysam_alignment.is_reverse,
                interval_core=interval_core,
            )
        ]
        # the ALT sequence of the SV signal is now in the forward direction and must not be reversed again.
        if merged_signal.sv_type >= 3 and len(merged_signal.get_alt_sequence()) == 0:
            raise ValueError(
                f"merged signal {merged_signal} has no alt_sequence. All breakends need to have an alt_sequence!"
            )

    return sorted(merged_sv_signals, key=lambda x: x.read_start)


# re-implement consensusAlignment.proto_svs=merged_svs. This is a simpler version, ignoring original_ref_sequences.
# it only parses the ref sequence corresponding to the merged signal chr,ref_start,ref_end
def add_ref_sequences_to_dict_alignments(
    dict_alignments: dict[str, list[consensus_class.ConsensusAlignment]],
    path_reference: str | Path,
    tmp_dir_path: Path | str | None = None,
) -> dict[str, list[consensus_class.ConsensusAlignment]]:
    regions: list[tuple[str, int, int, str]] = []
    for consensusAlignments in dict_alignments.values():
        for consensusAlignment in consensusAlignments:
            for merged_signal in consensusAlignment.proto_svs:
                if merged_signal.sv_type == 1:  # only consider deletions
                    region_string = f"{merged_signal.chr}:{merged_signal.ref_start}-{merged_signal.ref_end}"
                    regions.append((
                        merged_signal.chr,
                        merged_signal.ref_start,
                        merged_signal.ref_end,
                        region_string,
                    ))
                elif merged_signal.sv_type == 3:  # BND left
                    region_string = f"{merged_signal.chr}:{merged_signal.ref_start}-{merged_signal.ref_start + 1}"
                    regions.append((
                        merged_signal.chr,
                        merged_signal.ref_start,
                        merged_signal.ref_start + 1,
                        region_string,
                    ))
                elif merged_signal.sv_type == 4:  # BND right
                    region_string = f"{merged_signal.chr}:{merged_signal.ref_start - 1}-{merged_signal.ref_start}"
                    regions.append((
                        merged_signal.chr,
                        merged_signal.ref_start - 1,
                        merged_signal.ref_start,
                        region_string,
                    ))
    # make regions unique and sorted
    regions = sorted(set(regions), key=lambda x: (x[0], x[1]))
    # log.info(f"{len(regions)} unique regions to extract ref sequences from")
    # build a set of region tuples that are sorted and then written as region strings to a tmp bed file
    # use getfasta to get the sequences
    # build dict of interval:sequence from the resulting fasta file
    tmp_bed = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path,
        mode="w",
        delete=False if tmp_dir_path else True,
        suffix=".bed",
    )
    with open(tmp_bed.name, "wt") as f:
        for chr, start, end, region in regions:
            if end - start > 0:
                print(chr, start, end, region, sep="\t", file=f)
            else:
                raise ValueError(
                    f"end-start must be > 0. Got start={start}, end={end}, end-start={end - start}"
                )
    cmd_getfasta = [
        "bedtools",
        "getfasta",
        "-fi",
        str(path_reference),
        "-bed",
        tmp_bed.name,
    ]
    tmp_alt_fasta = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta")
    with open(tmp_alt_fasta.name, "wt") as f:
        # log.info(f"extracting ref sequences into {tmp_alt_fasta.name}...")
        subprocess.check_call(cmd_getfasta, stdout=f)
    # load fasta file # record name is the region string. ref_sequences is a dict of region:sequence
    ref_sequences = {
        record.name: str(record.seq)
        for record in SeqIO.parse(tmp_alt_fasta.name, "fasta")
    }

    # iterate over all merged signals and add the sequences.
    for consensusAlignments in dict_alignments.values():
        for consensusAlignment in consensusAlignments:
            for mss in consensusAlignment.proto_svs:
                # if the merged signal was merged from more than one original signal, add the original sequences as well
                if mss.sv_type == 0:
                    mss.original_ref_sequences = []
                    continue
                elif mss.sv_type == 1:
                    # add ref sequences (ignoring original signals, just add for each proto_sv the corresponding ref sequence)
                    mss.original_ref_sequences = [
                        ref_sequences[
                            f"{mss.chr}:{mss.ref_start}-{mss.ref_end}"
                        ].upper()
                    ]
                elif mss.sv_type == 3:
                    mss.original_ref_sequences = [
                        ref_sequences[
                            f"{mss.chr}:{mss.ref_start}-{mss.ref_start + 1}"
                        ].upper()
                    ]
                elif mss.sv_type == 4:
                    mss.original_ref_sequences = [
                        ref_sequences[
                            f"{mss.chr}:{mss.ref_start - 1}-{mss.ref_start}"
                        ].upper()
                    ]
                else:
                    raise ValueError(
                        f"merged signal type must be 0 (INS), 1 (DEL), 3 (BND) or 4 (BND). Got {mss.sv_type}"
                    )
    # return the augmented list of merged signals
    return dict_alignments
