# %%
from pathlib import Path

from . import (
    consensus_containers_to_db,
    consensuses_to_svprimitives,
    sv_calling,
)

# %%

path_container = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG004/mini500/consensus/0/consensus.252.txt"
)
path_consensus_alignments_db = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG004/mini500/HG004.consensus_alignments.db"
)

# load the consensus container
container = next(
    consensus_containers_to_db.parse_crs_container_results(path=path_container)
)
consensusIDs = list(container.consensus_dicts.keys())
consensuses = container.consensus_dicts
all_consensus_alignments = (
    consensuses_to_svprimitives.load_consensus_alignments_from_database(
        consensusIDs=consensusIDs,
        name_reference="hs37d5",
        path_db=path_consensus_alignments_db,
    )
)

svPrimitives = []
for consensusID in consensusIDs:
    consensus = consensuses[consensusID]
    consensus_alignments = all_consensus_alignments[consensusID]
    svPrimitives.extend(
        sv_calling.svPrimitives_from_consensus(
            consensus=consensus, consensus_alignments=consensus_alignments, verbose=True
        )
    )

assert len(svPrimitives[1].genotypeMeasurements[0].supporting_reads) > 10
# DEBUG! consensus.intervals_cutread_alignments is somehow empty or bullshit. Check this!
# TODO: make this a test!
# %%
