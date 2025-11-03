# %%
# iterate signaldepths
# find easy-to-solve SV candidates
# filter signals, if coverage is above provided threshold (per chromosome)
# parallelize over chromosomes (?)
# %%
import typing
from pathlib import Path

import numpy as np

# %%
window_size: int = 1000
sv_gathering_range: int = 30
path_signaldepths = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d01/HG002.signaldepth"
)
# %%
iter_lines = open(path_signaldepths, "r")
# %%


def parse_line(line) -> list:
    l = line.rstrip().split("\t")
    chrID = int(l[0])
    ref_start = int(l[1])
    ref_end = int(l[2])
    svtype = l[3]
    read_start = int(l[4])
    read_end = int(l[5])
    svsize = int(l[6])
    readname = l[7]
    sampleID = int(l[8])
    forward = int(l[9])
    chrName = l[10]
    depth = int(l[11])
    return [
        chrID,
        ref_start,
        ref_end,
        svtype,
        read_start,
        read_end,
        svsize,
        readname,
        sampleID,
        forward,
        chrName,
        depth,
    ]


def in_range(a: list, b: list, w: int) -> bool:
    if a[0] != b[0]:
        return False
    return abs(a[1] - b[1]) <= w


def svs_similar_size(l0, l1, precision=0.01) -> bool:
    return (
        abs(l0[6] - l1[6]) / max(l0[6], l1[6]) <= precision or abs(l0[6] - l1[6]) <= 3
    )


def similar_sv(a, b, w, precision=0.01) -> bool:
    return a[3] == b[3] and in_range(a, b, w) and svs_similar_size(a, b, precision)


# connects SV signals with similar size and position and of same type (INS,DELL)
def get_connected_components(lines, window_size, precision=0.01) -> typing.List[set]:
    UF = datastructures.UnionFind(vertex_names=range(len(lines)))
    for i in range(len(lines)):
        for j in range(len(lines)):
            if similar_sv(lines[i], lines[j], window_size, precision):
                UF.union(i, j)
    # get connected components
    cc = UF.get_connected_components()
    return cc


# %%
# min_lines is the fraction of lines compared to their avg depth
def process_buffer(
    lines,
    sv_gathering_range: int,
    min_lines: float = 0.3,
    max_relative_size_difference: float = 0.01,
    max_relative_std: float = 0.05,
):
    # append index to each line
    lines = [[*line, i] for i, line in enumerate(lines)]
    cc = get_connected_components(lines=lines, window_size=sv_gathering_range)
    # for each connected component, get the lines
    # and create a new SVcandidate
    used_lines = set()
    for cluster in cc:
        clines = [lines[i] for i in cluster]
        # m is an array of arrays of some numericals of the lines in the cluster
        # m is: ref_start, svsize, depth, sampleID, local_index
        m = np.array(
            [[line[1], line[6], line[11], line[8], line[-1]] for line in clines]
        )
        # process the cluster.
        # 1. calc std of sv sizes
        # 2. skip a cluster if its less than min_lines lines
        median_size = np.median(m[:, 1])
        relative_size_differences = np.abs((m[:, 1] - median_size) / m[:, 1])
        if np.median(relative_size_differences) > max_relative_size_difference:
            continue
        # std should not be too high
        if np.std(m[:, 1]) / median_size > max_relative_std:
            continue
        # there should be at least one sample which has min_lines lines relative to the median coverage
        for sampleID in np.unique(m[:, 3]):
            if np.sum(m[:, 3] == sampleID) / len(m) < min_lines:
                continue
        # if the cluster is not skipped, create a SVcandidate and save the used lines to a set
        svtypes_dict = {"DELL": "DEL", "INS": "INS"}
        readnames = [line[7] for line in clines]
        # depth is calculated for each sampleID.
        # e.g. {""ref_cov_l"": 32, ""ref_cov_r"": 32, ""alt_cov_l"": 1, ""alt_cov_r"": 1, ""score_l"": 0.0, ""score_r"": 0.0}
        # just score_r can't be solved here, so its just the same as score_l.
        # ref_cov_l is the coverage in m. ref_cov_r will be the same. alt_cov_l is the number of lines of this sampleID / median coverage.
        depths = {sampleID: {} for sampleID in np.unique(m[:, 3])}
        for sampleID in np.unique(m[:, 3]):
            mm = m[m[:, 3] == sampleID]
            depths[sampleID]["ref_cov_l"]: int(np.median(mm[:, 2]))
            depths[sampleID]["ref_cov_r"]: depths[sampleID]["ref_cov_l"]
            depths[sampleID]["alt_cov_l"]: mm.shape[0]
            depths[sampleID]["alt_cov_r"]: depths[sampleID]["alt_cov_l"]
            depths[sampleID]["score_l"]: 1.0
            depths[sampleID]["score_r"]: 1.0
        attributes = {
            "mcrID": -1,
            "assembly": 0,
            "readnames": readnames,
            "depths": depths,
        }
        # NOOOOOO, this can only work if refA and refB are the same! skip this for now...


#        svc = datatypes.SVcandidate(
#            SVType=svtypes_dict[clines[0][3]],
#            rafs=[0],
#            attributes

# %%

account = window_size
# find the first line that is not a DELR, BNDL, BNDR
while True:
    current = parse_line(iter_lines.readline())
    if current[3] in ["DELR", "BNDL", "BNDR"] or current[6] < 30:
        continue
    break
# initialize buffer
buffer = [current]
while True:
    try:
        current = parse_line(iter_lines.readline())
    except:
        # end of file
        process_buffer(buffer, sv_gathering_range)
        break
    # 'DELR','BNDL','BNDR' are just ignored
    if current[3] in ["DELR", "BNDL", "BNDR"]:
        continue
    gain = sum([similar_sv(current, item) for item in buffer])
    # add current to buffer
    buffer.append(current)
    # update account
    account -= 1 + int(round(np.log2(gain)))
    print("debug: account=", account)
    # if account is empty, process buffer and reset account
    if account <= 0:
        process_buffer(buffer, sv_gathering_range)
        account = window_size
        buffer = [current]

# %%


def test_in_range__various():
    objects = [
        ([0, 10], [0, 20], 20, True),
        ([0, 10, 12, 8], [0, 20, 55, 5], 20, True),
        ([0, 10], [0, 20], 5, False),
        ([0, 200], [0, 100], 100, True),
    ]
    results = [in_range(a, b, w) for a, b, w, _ in objects]
    expected = [e for _, _, _, e in objects]
    assert results == expected


# %%


def test_svs_similar_size__various():
    objects = [
        (
            [0, 100, 200, "DEL", 1010, 1200, 150, "readname_00", 0, 1, "chr1", 30],
            [0, 100, 200, "DEL", 1010, 1200, 150, "readname_01", 0, 1, "chr1", 30],
        ),
        (
            [0, 100, 200, "INS", 1010, 1200, 150, "readname_00", 0, 1, "chr1", 30],
            [0, 100, 200, "INS", 1010, 1200, 150, "readname_01", 0, 1, "chr1", 30],
        ),
        (
            [0, 100, 200, "DEL", 1010, 1200, 100, "readname_00", 0, 1, "chr1", 30],
            [0, 100, 200, "DEL", 1010, 1200, 150, "readname_01", 0, 1, "chr1", 30],
        ),
        (
            [0, 100, 200, "INS", 1010, 1200, 150, "readname_00", 0, 1, "chr1", 30],
            [0, 100, 200, "INS", 1010, 1200, 180, "readname_01", 0, 1, "chr1", 30],
        ),
    ]
    expected = [True, True, False, False]
    results = [svs_similar_size(a, b) for a, b in objects]
    assert results == expected


def test_similar_sv__simple():
    objects = [
        (
            [0, 100, 200, "DEL", 1010, 1200, 150, "readname_00", 0, 1, "chr1", 30],
            [0, 100, 200, "DEL", 1010, 1200, 150, "readname_01", 0, 1, "chr1", 30],
        ),
        (
            [0, 100, 200, "INS", 1010, 1200, 150, "readname_00", 0, 1, "chr1", 30],
            [0, 100, 200, "INS", 1010, 1200, 150, "readname_01", 0, 1, "chr1", 30],
        ),
    ]
    expected = [True, True]
    results = [similar_sv(a, b, 100) for a, b in objects]
    assert results == expected


def test_similar_sv__dels():
    w = 20
    precision = 0.01
    objects = [
        (
            [0, 100, 200, "DEL", 1010, 1200, 150, "readname_00", 0, 1, "chr1", 30],
            [0, 101, 200, "DEL", 1010, 1200, 151, "readname_01", 0, 1, "chr1", 30],
        ),
        (
            [0, 100, 200, "DEL", 1010, 1200, 150, "readname_00", 0, 1, "chr1", 30],
            [0, 130, 200, "DEL", 1010, 1200, 151, "readname_01", 0, 1, "chr1", 30],
        ),
        (
            [0, 100, 200, "DEL", 1010, 1200, 110, "readname_00", 0, 1, "chr1", 30],
            [0, 100, 200, "DEL", 1010, 1200, 151, "readname_01", 0, 1, "chr1", 30],
        ),
        (
            [0, 100, 200, "DEL", 1010, 1200, 150, "readname_00", 0, 1, "chr1", 30],
            [0, 100, 200, "INS", 1010, 1200, 150, "readname_01", 0, 1, "chr1", 30],
        ),
    ]
    expected = [True, False, False, False]
    results = [similar_sv(a, b, w, precision) for a, b in objects]
    assert results == expected


def test_similar_sv__real_false():
    a = parse_line(
        "\t".join(
            "17      1219409 1219410 DELL    19438   19439   61      f06f3067-e4da-4ddc-84f3-08f6a81f8e45    0       1       18      32".split()
        )
    )
    b = parse_line(
        "\t".join(
            "17      1219390 1219391 DELL    12596   12597   52      d5aa1e4d-ab7c-4246-be28-7df1a7978e7b    0       0       18      32".split()
        )
    )
    w = 30
    precision = 0.01
    expected = False
    result = similar_sv(a, b, w, precision)
    assert result == expected


# %%
l = """
1680 17      1219389 1219390 DELL    11975   11976   55      7e28f17d-1f20-4f94-a8e7-13711e8d1eb6    0       0       18      32
1681 17      1219389 1219390 DELL    15271   15272   55      c6198a01-f2a1-446e-8ca5-8519365109ad    0       1       18      32
1682 17      1219389 1219390 DELL    1988    1989    55      285cf5f6-ca4e-4c8f-90e1-dbf692d07d8e    0       0       18      32
1683 17      1219389 1219390 DELL    2072    2073    55      a5e34912-ce0c-44fe-914b-111329ab1e3c    0       1       18      32
1684 17      1219389 1219390 DELL    39062   39063   55      ae9cfeb5-6cfa-4432-8a2e-c32cec2c0858    0       1       18      32
1685 17      1219389 1219390 DELL    4942    4943    55      b292a462-986d-4196-9673-39b3b76f99c0    0       1       18      32
1686 17      1219389 1219390 DELL    62541   62542   55      b849eed0-c425-410c-852a-3067ea364147    0       0       18      32
1687 17      1219389 1219390 DELL    83381   83382   55      f2512597-1a12-4ab0-a1c9-2ff183af4cd8    0       0       18      32
1688 17      1219390 1219391 DELL    12596   12597   52      d5aa1e4d-ab7c-4246-be28-7df1a7978e7b    0       0       18      32
1689 17      1219391 1219392 DELL    477     478     55      48172a6b-906d-451c-99fb-4260a5e3df35    0       1       18      32
1690 17      1219393 1219394 DELL    14937   14938   56      df55b318-77ef-4bc8-847e-f5527fcc49e7    0       1       18      32
1691 17      1219397 1219398 DELL    52344   52345   56      c696c72a-dcc7-4991-bf28-80105741cc14    0       1       18      32
1692 17      1219398 1219399 DELL    5032    5033    54      a493c88a-d158-46f5-8c51-0a4900323757    0       0       18      32
1693 17      1219404 1219405 DELL    8531    8532    56      e09ffc68-6ccf-4b89-adec-826088f92e70    0       0       18      32
1694 17      1219408 1219409 DELL    29962   29963   56      889f2123-428c-4d33-8b7d-3b550a806666    0       1       18      32
1695 17      1219408 1219409 DELL    3510    3511    59      e9205a19-62ee-4067-984c-7af69ad71467    0       1       18      32
1696 17      1219408 1219409 DELL    52902   52903   54      fab8d8e2-8a7b-42d3-ba47-0f4c461710bf    0       0       18      32
1697 17      1219409 1219410 DELL    11534   11535   60      ab4bdada-5855-4a1f-9ac8-85eab6464e4f    0       1       18      32
1698 17      1219409 1219410 DELL    19438   19439   61      f06f3067-e4da-4ddc-84f3-08f6a81f8e45    0       1       18      32
1699 17      1219410 1219411 DELL    36305   36306   59      b548ab82-4eeb-412d-9c3e-191d58e6c144    0       1       18      32
1700 17      1219413 1219414 DELL    40144   40145   56      50736577-e244-4e51-bdf7-6a3f3fb32154    0       1       18      32
1701 17      1219413 1219414 DELL    5744    5745    58      48f6db22-6259-41c6-bd87-2c4732f03b28    0       0       18      32
1702 17      1219414 1219415 DELL    11081   11082   54      738a1738-f59b-4129-8dfa-edbd650cd636    0       1       18      32
1703 17      1219415 1219416 DELL    35735   35736   59      c533df82-0330-41f2-93af-f34644dc6189    0       1       18      32
1704 17      1219421 1219422 DELL    80204   80205   54      f241b015-07ff-4ccf-95b0-b064f109adfa    0       1       18      32
1705 17      1219430 1219431 DELL    45416   45417   56      35c10d83-bad6-4d40-86f4-d6ca0aab4df8    0       1       18      32
1706 17      1219431 1219432 DELL    18415   18416   57      9d00eb02-693e-498a-a019-d6a68409584c    0       0       18      32
1707 17      1219431 1219432 DELL    54597   54598   54      ff3030f3-aab0-4c5c-bf17-f333293e21b5    0       1       18      32
""".split(
    "\n"
)
lines = [parse_line("\t".join(line.split()[1:])) for line in l if line != ""]
# %%
