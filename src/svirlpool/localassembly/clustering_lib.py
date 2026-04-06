
from ..util import datatypes

def cluster_reads_by_reference_alignments_signals(
    crs:list[datatypes.CandidateRegion]
)
    # rework summed indels a bit.
    # construct a series of boxes, where a box is normal sequence or a repeat.
    # within each box, the ins and dels are summed for each read.
    # two reads are then compared by overlap of the boxes and similarity of the summed ins and dels within the boxes.
    # each box is weighted by the mean abs(ins-del) sizes per read within it.
    # either normalize all boxes or each box. include a switch to change logic here.
    # This results in a vector of scores for each read. Since reads can span several crs, the boxes need to be constructed
    # for all crs and sorted, e.g. by crID.
    # the alignment scores can be extracted from the crs.
    # if there is no signal in a box for a given read but the read has an alignment interval that spans the whole box, then
    # the read has no SV signal in that box.
    # If the read has no signals in a box and does not overlap with its alignment over the box, then treat this as 0.
    