# %%

import cProfile
import io
import pstats
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from intervaltree import Interval, IntervalTree
from logzero import logger as log
from pysam import AlignedSegment

from . import consensus, consensus_class, datastructures, datatypes

# %%


def profile_parse_sv_signals_from_consensus(
    samplename: str,
    consensusAlignment: consensus_class.ConsensusAlignment,
    consensus_sequence: str,
    interval_core: tuple[int, int] | None = None,
    min_signal_size: int = 5,
    min_bnd_size: int = 200,
    profile_output_path: Path | str | None = None,
) -> tuple[list[datatypes.MergedSVSignal], str]:
    """
    Profile the parse_sv_signals_from_consensus function and return results with profiling stats.

    Args:
        Same as parse_sv_signals_from_consensus
        profile_output_path: Optional path to save detailed profiling stats

    Returns:
        tuple: (results from function, profiling stats as string)
    """
    # Create a profiler
    profiler = cProfile.Profile()

    # Run the function under profiling
    profiler.enable()
    results = parse_sv_signals_from_consensus(
        samplename=samplename,
        consensusAlignment=consensusAlignment,
        consensus_sequence=consensus_sequence,
        interval_core=interval_core,
        min_signal_size=min_signal_size,
        min_bnd_size=min_bnd_size,
    )
    profiler.disable()

    # Generate stats
    stats_stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stats_stream)
    stats.sort_stats("cumulative")
    stats.print_stats(20)  # Show top 20 functions

    profiling_output = stats_stream.getvalue()

    # Optionally save detailed stats to file
    if profile_output_path:
        with open(profile_output_path, "w") as f:
            f.write(profiling_output)
        log.info(f"Detailed profiling stats saved to {profile_output_path}")

    return results, profiling_output


# %%


def ref_and_read_start_from_merging_sv_signals(
    sv_signals: list[datatypes.MergedSVSignal],
) -> tuple[int, int]:
    """returns the ref_start and read_start of the first signal in the list. Can be more elaborate in the future"""
    sv_signals = sorted(sv_signals, key=lambda x: x.ref_start)
    merged_ref_start = sv_signals[0].ref_start
    merged_read_start = sv_signals[0].read_start
    return (merged_ref_start, merged_read_start)


def merge_sv_signals(
    chr: str, sv_signals: list[datatypes.SVsignal], repeatIDs: list[int]
) -> datatypes.MergedSVSignal:
    """Merges sv_signals. Does not accept svs with type 3 or 4 (break ends).

    Args:
        sv_signals (list[datatypes.SVsignal]): list of sv_signals to merge
        repeatID (list[int]|None): IDs of the repeat where the signals are located.

    Returns:
        datatypes.MergedSVSignal
    """
    # assert that repeatIDs are all int and >= 0
    if len(repeatIDs) == 0:
        assert all(
            isinstance(repeatID, int) and repeatID >= 0 for repeatID in repeatIDs
        ), f"all repeatIDs must be int >= 0. Got {repeatIDs}"
    assert len(sv_signals) > 0, "sv_signals must have at least 1 signal"
    assert all(sv_signal.sv_type <= 1 for sv_signal in sv_signals), (
        f"all sv_signal types must be either 0 (INS) or 1 (DEL). Got {[sv_signal.sv_type for sv_signal in sv_signals]}"
    )
    # sort sv_signals by ref_start position
    # get new length of the merged signal. add if the sv signal is an insertion, subtract if it is a deletion
    sum_ins = sum([
        sv_signal.size for sv_signal in sv_signals if sv_signal.sv_type == 0
    ])
    sum_del = sum([
        sv_signal.size for sv_signal in sv_signals if sv_signal.sv_type == 1
    ])
    total_length = sum_ins - sum_del

    merged_type = 0 if total_length > 0 else 1
    merged_ref_start, merged_read_start = ref_and_read_start_from_merging_sv_signals(
        sv_signals
    )
    merged_ref_end = (
        merged_ref_start + total_length if merged_type == 0 else merged_ref_start
    )
    merged_read_end = (
        merged_read_start + total_length if merged_type == 0 else merged_read_start
    )
    return datatypes.MergedSVSignal(
        chr=chr,
        ref_start=merged_ref_start,
        ref_end=merged_ref_end,
        read_start=merged_read_start,
        read_end=merged_read_end,
        size=abs(total_length),
        sv_type=merged_type,
        original_signals=[svs.unstructure() for svs in sv_signals],
        repeatIDs=list(set(repeatIDs)) if len(repeatIDs) > 0 else None,
    )


# def add_alt_sequences_to_MergedSVSignal(merged_signal:datatypes.MergedSVSignal, consensus_sequence:str, reverse:bool) -> tuple[str,list[str]|None]:
#     assert len(consensus_sequence) > 0, "consensus_sequence must have a length > 0"
#     # there are different strategies for insertions and deletions. The simplest one is to
#     # cut the consensus sequence at the ref_start position and add the merged signal size
#     # iterate SV signals and add an alt sequence for each one. '' for each DEL and the sequence for each INS
#     if len(merged_signal.original_signals) > 0:
#         original_alt_sequences:list[str] = []
#         for sv_signal_unstructured in merged_signal.original_signals:
#             sv_signal = cattrs.structure(sv_signal_unstructured, datatypes.SVsignal)
#             original_alt_sequences.append('')
#             if sv_signal.sv_type == 0:
#                 alt_seq = Seq(consensus_sequence[sv_signal.read_start:sv_signal.read_end])
#                 original_alt_sequences.append(str(alt_seq.reverse_complement()) if reverse else str(alt_seq))
#         merged_alt_sequence:str = ''
#         if merged_signal.sv_type == 0:
#             alt_seq = Seq(consensus_sequence[merged_signal.read_start:merged_signal.read_start+merged_signal.size])
#             merged_alt_sequence:str = str(alt_seq.reverse_complement()) if reverse else str(alt_seq)
#         return merged_alt_sequence, original_alt_sequences
#     else:
#         if merged_signal.sv_type == 0:
#             alt_seq = Seq(consensus_sequence[merged_signal.read_start:merged_signal.read_end])
#             return str(alt_seq.reverse_complement()) if reverse else str(alt_seq), []
#         else:
#             return '', []


def alt_sequence_for_MergedSVSignal(
    merged_signal: datatypes.MergedSVSignal,
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


# # horizontal merge
# def merge_svs_in_dict_alignments(
#         dict_alignments:dict[str,list[datatypes.ConsensusAlignment]],
#         kmer_letter_dict:dict[str,int],
#         verbose:bool=False,
#         max_low_complexity_ignored_size:int=15,
#         low_complexity_max_2mer_count:int=4,
#         intersection_size_kmer_sketches:float=0.5,
#         merge_large_dels_with_small_gaps:tuple[int,int]=(50,500)) -> None:
#     # merges proto SVs (MergedSVSignal objects) in all consensusAlignments in dict_alignments
#     # following some criteria
#     # iterate all ConsensusAlignments and for each one, build union find data structures
#     # and connect all SV signals that share a repeatID
#     # then test each connected component, if all members fulfill the criteria for merging
#     # and merge them. If not, add them as singletons
#     for consensusID,consensusAlignments in dict_alignments.items():
#         for consensusAlignment in consensusAlignments:
#             svs:list[datatypes.MergedSVSignal] = consensusAlignment.proto_svs
#             # make all sv.repeatIDs sets
#             for sv in svs:
#                 if sv.sv_type >= 3: # don't horizontally merge BNDs
#                     continue
#                 if sv.repeatIDs is None or len(sv.repeatIDs) == 0:
#                     sv.repeatIDs = []
#             # count all 2-mers and save their counts in a dict index:dict[2-mer:count]
#             dict_2mer_sets: dict[int,set[int]] = dict()
#             for i,sv in enumerate(svs):
#                 if sv.sv_type == 0: # insertion
#                     seq = sv.get_alt_sequence()
#                 elif sv.sv_type == 1: # deletion
#                     seq = sv.get_ref_sequence()
#                 else:
#                     continue
#                 dict_2mer_sets[i] = util.kmer_set_from_string(string=seq,letter_dict=kmer_letter_dict,k=2)
#             # remove all sv signals if their 2-mer counts are <= 4. also re-index the dict_2mer_sets
#             filtered_sv_indeces = [i for i,sv in enumerate(svs) if (i not in dict_2mer_sets or len(dict_2mer_sets[i]) <= 4) and sv.size < max_low_complexity_ignored_size]
#             dict_2mer_sets_remapped = dict()
#             new_index = 0
#             for i,sv in enumerate(svs):
#                 if sv.sv_type >= 3:
#                     continue
#                 if i in filtered_sv_indeces:
#                     continue
#                 dict_2mer_sets_remapped[new_index] = dict_2mer_sets[i]
#                 new_index += 1
#             dict_2mer_sets = dict_2mer_sets_remapped
#             svs = [sv for i,sv in enumerate(svs) if i not in filtered_sv_indeces]

#             uf = datastructures.UnionFind(range(len(svs)))
#             # connect all sv_signals that share a repeatID and that either have only up to 4 2mers or both have more 2-mers
#             for i in range(len(svs)):
#                 if svs[i].sv_type >= 3:
#                     continue # skip break ends for now
#                 if i not in dict_2mer_sets:
#                     if verbose:
#                         print(f"no 2mers for signal {i}, size={svs[i].size}, ref_pos={svs[i].ref_start}")
#                     continue
#                 for j in range(i+1,len(svs)):
#                     if svs[j].sv_type >= 3:
#                         continue # skip break ends for now
#                     if j not in dict_2mer_sets or len(dict_2mer_sets[j]) == 0:
#                         if verbose:
#                             print(f"no 2mers for signal {j}, size={svs[j].size}, ref_pos={svs[j].ref_start}")
#                         continue
#                     # if the SVs to be merged are very large deletions and are very close (within close bp), then merge them
#                     if svs[i].sv_type == 1 and svs[j].sv_type == 1 and abs(svs[i].size) > merge_large_dels_with_small_gaps[1] and abs(svs[j].size) > merge_large_dels_with_small_gaps[1]:
#                         if Interval(svs[i].ref_start,svs[i].ref_end).distance_to(Interval(svs[j].ref_start,svs[j].ref_end)) <= merge_large_dels_with_small_gaps[0]:
#                             if verbose:
#                                 print(f"union by distance {i} and {j} because both are large evenets (>1kb) and are close")
#                             uf.union_by_name(i,j)
#                             continue
#                     if len(set(svs[i].repeatIDs).intersection(set(svs[j].repeatIDs))) > 0:
#                         # if both i and j have max low_complexity_max_2mer_count 2mers
#                         if len(dict_2mer_sets[j]) <= low_complexity_max_2mer_count and len(dict_2mer_sets[i]) <= low_complexity_max_2mer_count \
#                                 or len(dict_2mer_sets[j]) > 3 and len(dict_2mer_sets[i]) > low_complexity_max_2mer_count:
#                             # check similarity of 2mers:
#                             intersection_size = util.intersection_ratio_of_smaller_set(dict_2mer_sets[i],dict_2mer_sets[j])
#                             if intersection_size > intersection_size_kmer_sketches:
#                                 if verbose:
#                                     print(f"union by name {i} and {j}")
#                                 uf.union_by_name(i,j)
#                             else:
#                                 if verbose:
#                                     print(f"no matching 2mers overlap. Len of 2mers of {i}={len(dict_2mer_sets[i])}, {j}={len(dict_2mer_sets[j])}")
#                         else:
#                             if verbose:
#                                 print(f"no matching 2mers overlap. Len of 2mers of {i}={len(dict_2mer_sets[i])}, {j}={len(dict_2mer_sets[j])}")
#                     else:
#                         if verbose:
#                             print(f"no repeatID overlap between {i} and {j}")
#             # get connected components
#             # merge each connected component of svs
#             # save the newly merged signals in consensusAlignment.proto_svs
#             cc = uf.get_connected_components(allow_singletons=True)
#             # start new SVs with all bnds
#             # This doesn't make sense right now!
#             new_svs = []
#             for indices in cc:
#                 signals_to_merge = [svs[i] for i in indices if svs[i].sv_type < 3]
#                 if len(signals_to_merge) > 0:
#                     new_svs.append(merge_merged_svs(signals_to_merge))
#             bnds = [sv for sv in svs if sv.sv_type >= 3]
#             consensusAlignment.proto_svs = sorted(new_svs + bnds, key=lambda x: x.ref_start)


# to debug and test merge_svs_in_dict_alignments, we need to save all input to a pickled file
def save_merge_svs_in_dict_alignments_input(input, path: Path | str):
    """Saves the input for merge_svs_in_dict_alignments to a file."""
    if isinstance(path, str):
        path = Path(path)
    with open(path, "wb") as f:
        import pickle

        pickle.dump(input, f)
    log.info(f"Saved input for merge_svs_in_dict_alignments to {path}")


# --- horizontal merge --- #
# test this!
def merge_svs_in_dict_alignments(
    dict_alignments: dict[str, list[consensus_class.ConsensusAlignment]],
    max_low_complexity_ignored_size: int,
    merge_large_dels_with_small_gaps: tuple[int, int] = (50, 500),
    verbose: bool = False,
    minmax_GC_tolerance: float = 0.0,
) -> None:
    # merges proto SVs (MergedSVSignal objects) in all consensusAlignments in dict_alignments
    # following some criteria
    # iterate all ConsensusAlignments and for each one, build union find data structures
    # and connect all SV signals that share a repeatID
    # then test each connected component, if all members fulfill the criteria for merging
    # and merge them. If not, add them as singletons
    for _consensusID, consensusAlignments in dict_alignments.items():
        for consensusAlignment in consensusAlignments:
            svs: list[datatypes.MergedSVSignal] = consensusAlignment.proto_svs
            # make all sv.repeatIDs sets
            for sv in svs:
                if sv.sv_type >= 3:  # don't horizontally merge BNDs
                    continue
                if sv.repeatIDs is None or len(sv.repeatIDs) == 0:
                    sv.repeatIDs = []
            # count all 2-mers and save their counts in a dict index:dict[2-mer:count]
            dict_gc_content: dict[int, float] = {}
            for i, sv in enumerate(svs):
                if sv.sv_type == 0:  # insertion
                    dict_gc_content[i] = SeqUtils.gc_fraction(sv.get_alt_sequence())
                elif sv.sv_type == 1:  # deletion
                    dict_gc_content[i] = SeqUtils.gc_fraction(sv.get_ref_sequence())
                else:
                    continue  # break end
            # remove all sv signals if their 2-mer counts are <= 4. also re-index the dict_2mer_sets
            filtered_sv_indeces = [
                i
                for i, sv in enumerate(svs)
                if sv.sv_type >= 3
                or i not in dict_gc_content
                or (
                    (
                        dict_gc_content[i] > (1.0 - minmax_GC_tolerance)
                        or dict_gc_content[i] < minmax_GC_tolerance
                    )
                    and sv.size <= max_low_complexity_ignored_size
                )
            ]
            uf = datastructures.UnionFind(range(len(svs)))
            for i in range(len(svs)):
                if i in filtered_sv_indeces:
                    continue
                for j in range(i + 1, len(svs)):
                    if j in filtered_sv_indeces:
                        continue
                    if (
                        svs[i].sv_type == 1
                        and svs[j].sv_type == 1
                        and abs(svs[i].size) > merge_large_dels_with_small_gaps[1]
                        and abs(svs[j].size) > merge_large_dels_with_small_gaps[1]
                    ):
                        if (
                            Interval(svs[i].ref_start, svs[i].ref_end).distance_to(
                                Interval(svs[j].ref_start, svs[j].ref_end)
                            )
                            <= merge_large_dels_with_small_gaps[0]
                        ):
                            if verbose:
                                print(
                                    f"union by distance {i} and {j} because both are large evenets (>1kb) and are close"
                                )
                            uf.union_by_name(i, j)
                    elif (
                        len(set(svs[i].repeatIDs).intersection(set(svs[j].repeatIDs)))
                        > 0
                    ):
                        if verbose:
                            print(
                                f"union by name {i} and {j} with shared repeatIDs: {set(svs[i].repeatIDs).intersection(set(svs[j].repeatIDs))}"
                            )
                        uf.union_by_name(i, j)
                        if svs[i].chr != svs[j].chr:
                            raise ValueError(
                                f"SVs {i} and {j} have different chromosomes: {svs[i].chr} and {svs[j].chr}. This is not allowed.\nsv a: {svs[i]}\nsv b: {svs[j]}"
                            )
                    else:
                        if verbose:
                            print(f"no union for {i} and {j}")
            # get connected components
            # merge each connected component of svs
            # save the newly merged signals in consensusAlignment.proto_svs
            cc = uf.get_connected_components(allow_singletons=True)
            # start new SVs with all bnds
            # This doesn't make sense right now!
            new_svs = []
            for indices in cc:
                signals_to_merge = [
                    svs[i] for i in indices if i not in filtered_sv_indeces
                ]
                if len(signals_to_merge) > 0:
                    new_svs.append(merge_merged_svs(signals_to_merge))
            bnds = [sv for sv in svs if sv.sv_type >= 3]
            consensusAlignment.proto_svs = sorted(
                new_svs + bnds, key=lambda x: x.ref_start
            )


# TODO: needs testing
def merge_merged_svs(
    signals_to_merge: list[datatypes.MergedSVSignal],
) -> datatypes.MergedSVSignal:
    # merging signals: all input signals must be represented in the new original signals
    # new total length is the sum of sum(insertions) - sum(deletions)
    # ref_start is ref_start of first signal, ref_end is ref_start + new total length
    # repeat IDs are all repeat IDs of the input signals (then cast to list)
    # alt_sequence or ref_sequence is a concatenation of all alt_sequences or ref_sequences,
    # subsetted to the new total length
    if len(signals_to_merge) == 0:
        raise ValueError("no signals to merge")
    if len(signals_to_merge) == 1:
        return signals_to_merge[0]
    for merged_signal in signals_to_merge:
        if merged_signal.sv_type >= 3:
            raise ValueError("merged signal must not be a BND signal")
    repeatIDs = list({
        repeatID
        for merged_signal in signals_to_merge
        for repeatID in merged_signal.repeatIDs
    })
    chr = signals_to_merge[0].chr
    total_size = sum([
        merged_signal.size
        for merged_signal in signals_to_merge
        if merged_signal.sv_type == 0
    ]) - sum([
        merged_signal.size
        for merged_signal in signals_to_merge
        if merged_signal.sv_type == 1
    ])
    original_alt_sequences = [
        seq
        for merged_signal in signals_to_merge
        for seq in merged_signal.original_alt_sequences
        if seq != ""
    ]
    original_ref_sequences = [
        seq
        for merged_signal in signals_to_merge
        for seq in merged_signal.original_ref_sequences
        if seq != ""
    ]

    sv_type = 0 if total_size > 0 else 1
    size = abs(total_size)
    # to define read_start and ref_start, find the largest contributing signal and take its position
    largest_signal = max(signals_to_merge, key=lambda x: x.size)
    ref_start = largest_signal.ref_start
    ref_end = ref_start + size if sv_type == 1 else ref_start + 1
    read_start = largest_signal.read_start
    read_end = read_start + size if sv_type == 0 else read_start + 1

    return datatypes.MergedSVSignal(
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=read_start,
        read_end=read_end,
        size=size,
        sv_type=sv_type,
        chr=chr,
        repeatIDs=repeatIDs,
        original_alt_sequences=original_alt_sequences,
        original_ref_sequences=original_ref_sequences,
    )


# horizontal merging
def merge_to_proto_svs(
    chr: str,
    it_repeats: IntervalTree,
    sv_signals: list[datatypes.SVsignal],
    dict_repeatIDs: dict[int, set[int]],
) -> list[datatypes.MergedSVSignal]:
    UF = datastructures.UnionFind(range(len(sv_signals)))
    # connect all sv_signals that share a repeatID
    for indices in dict_repeatIDs.values():
        indices = list(indices)
        # union_by_name all others with the first in the list
        for i in indices[1:]:
            # don't connect if the signal is a deletion composed of only 1 or 2 2-mers
            if sv_signals[i].sv_type == 2:  # deletion
                pass
            UF.union_by_name(indices[0], i)
        # get connected components
    connected_components = UF.get_connected_components(allow_singletons=True)
    merged_svs: list[datatypes.MergedSVSignal] = []
    for cc in connected_components:
        if len(cc) == 1:
            # just the corresponding SV signal to results
            signal: datatypes.SVsignal = sv_signals[cc.pop()]
            merged_sv_signal = datatypes.MergedSVSignal(
                ref_start=signal.ref_start,
                ref_end=signal.ref_end,
                read_start=signal.read_start,
                read_end=signal.read_end,
                size=signal.size,
                sv_type=signal.sv_type,
                chr=chr,
                original_signals=[signal.unstructure()],
            )
            repeatIDs = list({
                repeat.data for repeat in it_repeats[signal.ref_start : signal.ref_end]
            })
            if len(repeatIDs) > 0:
                merged_sv_signal.repeatIDs = repeatIDs
            merged_svs.append(merged_sv_signal)
        else:
            # bnds might overlap repeats and should not be merged. They are added to the list of merged signals as singletons
            bnd_signals = [sv_signals[i] for i in cc if sv_signals[i].sv_type > 2]
            for bnd_signal in bnd_signals:
                merged_sv_signal = datatypes.MergedSVSignal(
                    ref_start=bnd_signal.ref_start,
                    ref_end=bnd_signal.ref_end,
                    read_start=bnd_signal.read_start,
                    read_end=bnd_signal.read_end,
                    size=bnd_signal.size,
                    sv_type=bnd_signal.sv_type,
                    chr=chr,
                    original_signals=[bnd_signal.unstructure()],
                )
                repeatIDs = list({
                    repeat.data
                    for repeat in it_repeats[
                        bnd_signal.ref_start : bnd_signal.ref_end + 1
                    ]
                })
                if len(repeatIDs) > 0:
                    merged_sv_signal.repeatIDs = repeatIDs
                merged_svs.append(merged_sv_signal)
            selected_sv_signals = [
                sv_signals[i] for i in cc if sv_signals[i].sv_type <= 2
            ]
            if len(selected_sv_signals) == 0:
                continue
                # get all repeatIDs that are overlapping any of the original signals
            repeatIDs = set()
            for sv_signal in selected_sv_signals:
                for repeat in it_repeats[
                    sv_signal.ref_start : sv_signal.ref_end + 1
                ]:  # add +1 to end to avoid empty intervals
                    repeatIDs.add(repeat.data)
            merged_sv_signal = merge_sv_signals(
                chr=chr, sv_signals=selected_sv_signals, repeatIDs=list(repeatIDs)
            )
            merged_svs.append(merged_sv_signal)
    return merged_svs


# def proto_svs(consensusAlignment:datatypes.ConsensusAlignment,
#               consensus_sequence:str,
#               min_signal_size:int=8,
#               min_bnd_size:int=200) -> list[datatypes.MergedSVSignal]:
#         # horizontal merging
#         # build an Intervaltree of the trf intervals
#         it_repeats = IntervalTree()
#         dict_repeatIDs = dict()
#         for start,end,repeatID in consensusAlignment.trf_intervals:
#             it_repeats.addi(start,end,repeatID)
#             dict_repeatIDs[repeatID] = set()
#         # parse sv signals from alignments
#         pysam_alignment = consensusAlignment.alignment.to_pysam()
#         sv_signals:list[datatypes.SVsignal] = consensus_lib.parse_ReadAlignmentSignals_from_alignment(
#             alignment=pysam_alignment,
#             sampleID=0,
#             min_signal_size=min_signal_size,
#             min_bnd_size=min_bnd_size).SV_signals
#         # iterate i,sv_signals and create a dictionary {repeatID:{i}} for each signal. Add all repeatIDs to the dictionary
#         for i,sv_signal in enumerate(sv_signals):
#             for interval in it_repeats[sv_signal.ref_start:sv_signal.ref_end+1]: # add +1 to end to avoid empty intervals
#                 dict_repeatIDs[interval.data].add(i)
#         proto_svs = merge_to_proto_svs(
#             chr=consensusAlignment.alignment.reference_name,
#             it_repeats=it_repeats,
#             sv_signals=sv_signals,
#             dict_repeatIDs=dict_repeatIDs)
#         # add alt sequences to all MergedSVSignals
#         for merged_signal in proto_svs:
#             alt_sequence, original_alt_sequences = \
#                 add_alt_sequences_to_MergedSVSignal(merged_signal=merged_signal, consensus_sequence=consensus_sequence, reverse=pysam_alignment.is_reverse)
#             # check if the alt_sequence is empty. if also sv_type is 0, then raise an error
#             if merged_signal.sv_type == 0 and len(alt_sequence) == 0:
#                 raise ValueError(f"merged signal {merged_signal} has sv_type 0 (INS) but no alt_sequence. Any insertion needs to have an alt_sequence.")
#             merged_signal.alt_sequence = alt_sequence
#             merged_signal.original_alt_sequences = original_alt_sequences
#         return sorted(proto_svs)


def assign_repeat_ids(
    sv_signals: list[datatypes.MergedSVSignal],
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
    consensusAlignment: consensus_class.ConsensusAlignment,
    consensus_sequence: str,
    interval_core: tuple[int, int] | None = None,
    min_signal_size: int = 5,
    min_bnd_size: int = 200,
    _pysam_cache: dict[int, AlignedSegment] = None,
) -> list[datatypes.MergedSVSignal]:
    """
    Parse SV signals from consensus alignment.

    Args:
        _pysam_cache: Optional cache for pysam alignments to avoid repeated to_pysam() calls
    """
    # Cache pysam alignment to avoid repeated expensive conversions
    if _pysam_cache is None:
        _pysam_cache = {}

    if consensusAlignment.uid not in _pysam_cache:
        _pysam_cache[consensusAlignment.uid] = consensusAlignment.alignment.to_pysam()
    pysam_alignment = _pysam_cache[consensusAlignment.uid]

    # parse sv signals from alignments
    sv_signals: list[datatypes.SVsignal] = (
        consensus.parse_ReadAlignmentSignals_from_alignment(
            samplename=samplename,
            alignment=pysam_alignment,
            min_signal_size=min_signal_size,
            min_bnd_size=min_bnd_size,
        ).SV_signals
    )

    # if interval_core is given: this describes the actual start and end of the consensus sequence
    # filter any signals whose read_start and read_end are not fully inside interval_core
    if interval_core is not None:
        sv_signals = [
            sv_signal
            for sv_signal in sv_signals
            if sv_signal.read_start > interval_core[0]
            and sv_signal.read_end < interval_core[1]
        ]

    # transform all SVsignals to MergedSVSignals
    merged_sv_signals = [
        datatypes.MergedSVSignal(
            chr=consensusAlignment.alignment.reference_name,
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

    assign_repeat_ids(merged_sv_signals, consensusAlignment.trf_intervals)

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

    return sorted(merged_sv_signals, key=lambda x: x.ref_start)


# # def add_ref_sequences_to_MergedSVSignals(merged_svs:list[datatypes.MergedSVSignal], path_reference:str) -> list[datatypes.MergedSVSignal]:
# def add_ref_sequences_to_dict_alignments(
#             dict_alignments:dict[str,list[datatypes.ConsensusAlignment]],
#             path_reference:str,
#             tmp_dir_path:Path|None=None) -> dict[str,list[datatypes.ConsensusAlignment]]:
#     """Adds ref_sequence to all merged sv signals. If the signal was merged from more than one original sv signal,
# then original_ref_sequences is also added."""
#     # procedure:
#     # for each signal, generate the region string like chr:start-end and the region tuple (chr,start,end)
#     regions:list[tuple[str,int,int,str]] = []
#     for consensusID,consensusAlignments in tqdm(dict_alignments.items()):
#         for consensusAlignment in consensusAlignments:
#             for merged_signal in consensusAlignment.proto_svs:
#                 region_string = f"{merged_signal.chr}:{merged_signal.ref_start}-{merged_signal.ref_end}"
#                 regions.append((merged_signal.chr, merged_signal.ref_start, merged_signal.ref_end, region_string))
#                 for sv_signal_unstructured in merged_signal.original_signals:
#                     sv_signal = cattrs.structure(sv_signal_unstructured, datatypes.SVsignal)
#                     region_string = f"{merged_signal.chr}:{sv_signal.ref_start}-{sv_signal.ref_end}"
#                     regions.append((merged_signal.chr, sv_signal.ref_start, sv_signal.ref_end, region_string))
#     # make regions unique and sorted
#     regions = list(set(regions))
#     log.info(f"{len(regions)} unique regions to extract ref sequences from")

#     # build a set of region tuples that are sorted and then written as region strings to a tmp bed file
#     # use getfasta to get the sequences
#     # build dict of interval:sequence from the resulting fasta file
#     tmp_bed = tempfile.NamedTemporaryFile(dir=tmp_dir_path, mode='w',delete=False if tmp_dir_path else True, suffix='.bed')
#     with open(tmp_bed.name, 'wt') as f:
#         for chr,start,end,region in regions:
#             if end-start > 0:
#                 print(chr, start, end, region, sep='\t', file=f)
#     cmd_getfasta = ['bedtools','getfasta','-fi',str(path_reference),'-bed',tmp_bed.name]
#     tmp_alt_fasta = tempfile.NamedTemporaryFile(mode='w',delete=False,suffix='.fasta')
#     with open(tmp_alt_fasta.name,'wt') as f:
#         log.info(f"extracting ref sequences into {tmp_alt_fasta.name}...")
#         subprocess.check_call(cmd_getfasta,stdout=f)
#     # load fasta file # record name is the region string. ref_sequences is a dict of region:sequence
#     ref_sequences = {record.name:str(record.seq) for record in SeqIO.parse(tmp_alt_fasta.name, "fasta")}

#     # iterate over all merged signals and add the sequences.
#     for consensusID,consensusAlignments in tqdm(dict_alignments.items()):
#         for consensusAlignment in consensusAlignments:
#             for mss in consensusAlignment.proto_svs:
#                 # if the merged signal was merged from more than one original signal, add the original sequences as well
#                 if mss.sv_type >= 3:
#                     mss.ref_sequence = ''
#                     continue
#                 # actually add ref sequences
#                 assert len(mss.original_signals) > 0, "merged signal must have at least one original signal"
#                 mss.original_ref_sequences = []
#                 for sv_signal_unstructured in mss.original_signals:
#                     sv_signal = cattrs.structure(sv_signal_unstructured, datatypes.SVsignal)
#                     if sv_signal.ref_start == sv_signal.ref_end:
#                         mss.original_ref_sequences.append('')
#                     else:
#                         mss.original_ref_sequences.append(ref_sequences[f"{mss.chr}:{sv_signal.ref_start}-{sv_signal.ref_end}"])
#                 # create one chimeric sequence for the merged signal if there are several original signals
#                 # the sequence chimera can only be as long as the merged signal size
#                 # first concat all sequences and then cut the sequence to the merged signal size
#                 chimera_sequence = ''.join(mss.original_ref_sequences)
#                 mss.ref_sequence = chimera_sequence[:mss.size]
#     # return the augmented list of merged signals
#     return dict_alignments


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
