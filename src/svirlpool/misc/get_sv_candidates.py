# load reads
# scan for deletion candidates
# 1) intra alignment
# 2) inter alignments
# %%
import csv
import json
from pathlib import Path

import cattrs
import numpy as np
import pandas as pd

from . import datatypes, util


# %%
def reads_to_informative_loci(input: Path, reference: Path) -> list:
    L = []
    with open(input, "r") as f:
        ref_dict = util.create_ref_dict(reference=reference)
        reader = csv.reader(f, delimiter="\t", quotechar='"')
        for line in reader:
            L.append(
                cattrs.structure(json.loads(line[2]), datatypes.ReadAlignmentsSequence)
            )
    return L


# %%

rpath = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/testcases/reads_sequences/pb-clr.32x.del_be.reads.tsv"
)
refpath = Path(
    "/fast/work/projects/cubit/18.12/static_data/reference/GRCh37/hs37d5/hs37d5.fa"
)

reads = reads_to_informative_loci(rpath, refpath)

# %%
r = reads[4]
# %%
# intra-deletion
L = []
for r in reads:
    for a in r.ReadAlignmentFragments:
        for d in a.deletions:
            L.append((a.referenceID, d.ref_pos, d.ref_pos + d.ref_length, d.ref_length))
# %%

# clustering approaches:
#   distance metric (c. SVIM)
#   vertex cover

df = pd.DataFrame(L)

# %%
import matplotlib.pyplot as plt


# %%
def draw_read(read: datatypes.ReadAlignmentsSequence):
    plt.figure()
    plt.title(read.ReadAlignmentFragments[0].read_name)
    A = sorted(read.ReadAlignmentFragments, key=lambda x: x.readAlignmentStart)
    for j, aln in enumerate(read.ReadAlignmentFragments):
        color = "red" if aln.alignmentForward else "blue"
        y = j  # int(aln.alignmentSecondary)
        x = (aln.readAlignmentStart, aln.readAlignmentEnd)
        plt.plot(x, (y, y), color=color)
        if aln.breakendLeft:
            plt.plot(
                (x[0] - aln.breakendLeft.clippedLength, x[0]), (y, y), color="black"
            )
        if aln.breakendRight:
            plt.plot(
                (x[1], x[1] + aln.breakendRight.clippedLength), (y, y), color="grey"
            )
        for i, ins in enumerate(aln.insertions):
            h = y + 0.1 * i
            # plt.plot((ins.read_pos,ins.read_pos+ins.read_length),(h,h),color='orange')
            plt.scatter(ins.read_pos + ins.read_length / 2, h, color="orange")
        for i, dl in enumerate(aln.deletions):
            h = y + 0.1 * i
            plt.scatter(dl.read_pos, h, color="green")


# %%

factor_supplementary = 1.0
factor_secondary = 1.0

# inter-deletions
# 1) find all that have left breakend
# 2) find all that have right breakend

deletion_candidates = []
for r in reads:
    L = []
    A = sorted(r.ReadAlignmentFragments, key=lambda x: x.readAlignmentStart)
    for i, a in enumerate(A):
        if not a.breakendRight:
            continue
        s = None
        for j in range(i, len(A)):
            b = A[j]
            if (
                a.referenceID != b.referenceID
                or a.alignmentForward != b.alignmentForward
            ):
                continue
            if a.breakendRight and b.breakendLeft:
                s = a.breakendRight
                e = b.breakendLeft
                delsize = s.referencePosition - e.referencePosition
                cscorea = (
                    factor_supplementary * a.alignmentSupplementary
                    + factor_secondary * a.alignmentSecondary
                )
                cscoreb = (
                    factor_supplementary * b.alignmentSupplementary
                    + factor_secondary * b.alignmentSecondary
                )
                if delsize:
                    print(delsize, s, e, a.alignmentForward)
                cscore = (cscorea + cscoreb) / 2 + (
                    (1 / np.log(delsize)) if delsize > 1 else 0
                )
                L.append(
                    datatypes.DeletionCandidate(
                        a.referenceID, s.referencePosition, delsize, cscore
                    )
                )
    if len(L) > 0:
        deletion_candidates.append(L)

# %%
