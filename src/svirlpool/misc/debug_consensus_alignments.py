# %%
from copy import deepcopy
from pathlib import Path

from tqdm import tqdm

from ..scripts import consensus_align_lib, consensuses_to_svprimitives, datatypes, util

# %%

consensus_to_reference_alignments = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/HG002/svirlpool/problems/consensus.bam"
)
path_trf = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/human_hs37d5.trf.bed"
)
name_reference = "hs37d5"
samplename = "HG002"
path_db = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/HG002/svirlpool/problems/consensuses.db"
)
path_reference = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/hs37d5.fa"
)

# %%

sequences: dict[str, str] = util.load_consensus_sequences_from_db(path_db)

consensus_alignments: dict[str, list[datatypes.Alignment]] = (
    consensuses_to_svprimitives.load_alignments(
        consensus_to_reference_alignments, "HG002"
    )
)

dict_alignments: dict[str, list[datatypes.ConsensusAlignment]] = (
    consensuses_to_svprimitives.add_trf_annotations_to_alignments(
        alignments=consensus_alignments,
        trf_intervals=consensuses_to_svprimitives.trf_to_interval_tree(path_trf),
        reference_name=name_reference,
    )
)

consensus_alignments = dict()  # clear memory
unaligned_consensusIDs = set(sequences.keys()) - set(dict_alignments.keys())

for consensusID, consensusAlignments in tqdm(dict_alignments.items()):
    sequence = sequences[consensusID]
    for consensusAlignment in consensusAlignments:
        # add proto_svs (MergedSVSignal)s to all ConsensusAlignments in dict_alignments
        # contains the alt sequences
        merged_svs: list[datatypes.MergedSVSignal] = (
            consensus_align_lib.parse_sv_signals_from_consensus(
                samplename=samplename,
                consensusAlignment=consensusAlignment,
                consensus_sequence=sequence,
            )
        )
        # TODO: chech if the merged_svs are correct
        consensusAlignment.proto_svs = merged_svs
# %%
dict_alignments = consensus_align_lib.add_ref_sequences_to_dict_alignments(
    dict_alignments=dict_alignments, path_reference=path_reference
)
# %%

dict_alignments_merged = deepcopy(dict_alignments)

kmer_letter_dict = util.compute_letter_dict("ACGTN")
# BNDs need to be processed beforehand. They are not merged with other SVs, but are reported separately
consensus_align_lib.merge_svs_in_dict_alignments(
    dict_alignments=dict_alignments_merged, kmer_letter_dict=kmer_letter_dict
)

# %%
# merge_svs_in_dict_alignments needs debugging!
dict_alignments_merged_test = deepcopy(dict_alignments)
dict_alignments_merged_test = {"10.0": dict_alignments_merged_test["10.0"]}
consensus_align_lib.merge_svs_in_dict_alignments(
    dict_alignments=dict_alignments_merged_test,
    kmer_letter_dict=kmer_letter_dict,
    intersection_size_kmer_sketches=0.5,
)

# %%
mmm = consensus_align_lib.merge_merged_svs(
    dict_alignments_merged_test["10.0"][0].proto_svs
)
# %%
