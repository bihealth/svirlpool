# %%
#!/usr/bin/env python
import argparse
import csv
import multiprocessing as mp
import pickle
import sqlite3
import sys
from pathlib import Path, PosixPath

import cattrs
import numpy as np
import numpy.typing as npt
from intervaltree import IntervalTree
from logzero import logger
from tqdm import tqdm

# %%
from . import datatypes, util

# %%

# # find stretches, accepts a np array of ints, returns a list of tuples
# def find_stretches(arr:np.ndarray) -> typing.List[typing.Tuple[int,int]]:
#     if len(arr) == 0:
#         return []
#     mask = arr
#     a = np.concatenate([mask,np.array([mask.max()+1])])
#     b = np.concatenate([np.array([mask.min()-1]),mask])
#     c = a != b
#     starts = np.where(c)[0]
#     ends = starts -1
#     return [list(x) for x in zip(starts[:-1],ends[1:],arr[starts[:-1]])]


def write_to_bedgraphs(
    signals: list[datatypes.ExtendedSVsignal],
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
        for i, signal in tqdm(enumerate(signals)):
            val = float(values_signal[i])
            if not (
                val is None
                or val == float("inf")
                or val == float("-inf")
                or np.isnan(val)
            ):
                print(signal.chr, signal.ref_start, signal.ref_end, str(val), file=f)
    with open(bedgraph.with_suffix(".normalized.bedgraph"), "w") as f:
        trackdef = "track type=bedGraph name=values_normalized description=values_normalized \
visibility=display_mode altColor=150,20,255 color=20,150,255 \
priority=10 autoScale=on alwaysZero=off gridDefault=on \
maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper \
windowingFunction=mean smoothingWindow=off"
        print(trackdef, file=f)
        for i, signal in tqdm(enumerate(signals)):
            val = float(values_normalized[i])
            if not (
                val is None
                or val == float("inf")
                or val == float("-inf")
                or np.isnan(val)
            ):
                print(signal.chr, signal.ref_start, signal.ref_end, str(val), file=f)
    with open(bedgraph.with_suffix(".masked.bedgraph"), "w") as f:
        trackdef = "track type=bedGraph name=values_signal_masked description=values_signal_masked \
visibility=display_mode altColor=150,50,50 color=250,150,150 \
priority=10 autoScale=on alwaysZero=off gridDefault=on \
maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper \
windowingFunction=mean smoothingWindow=off"
        print(trackdef, file=f)
        for i, signal in tqdm(enumerate(signals)):
            val = float(values_masked[i])
            if not (
                val is None
                or val == float("inf")
                or val == float("-inf")
                or np.isnan(val)
            ):
                print(
                    signal.chr,
                    signal.ref_start,
                    signal.ref_end,
                    str(float(values_masked[i])),
                    file=f,
                )
    logger.info("Finished bedgraphs.")


# def effective_intervals_intervaltrees(rafs:list[Path], readnames_dict, chromosomes) -> dict[str,IntervalTree]:
#     # Pre-allocate per-chromosome buffers
#     per_chr: dict[str, list[tuple[int, int, int]]] = {chr: [] for chr in chromosomes}
#     get_read_id = readnames_dict.get

#     for path_raf in rafs:
#         logger.info(f"Load all rafs from {path_raf}..")
#         all_rafs = list(util.yield_from_raf(input=path_raf))
#         for raf in all_rafs:
#             rid = get_read_id(raf.read_name)
#             if rid is None:
#                 continue  # skip reads without significant signal
#             chr_, start, end = raf.effective_interval  # expect tuple (chr, start, end)
#             # Make sure chromosome exists (skip unexpected contigs)
#             if chr_ in per_chr:
#                 per_chr[chr_].append((start, end, rid))

#     # Build interval trees in bulk
#     raf_trees: dict[str, IntervalTree] = {}
#     for chr_ in chromosomes:
#         tuples_ = per_chr[chr_]
#         if tuples_:
#             # from_tuples expects (begin, end, data)
#             raf_trees[chr_] = IntervalTree.from_tuples(tuples_)
#         else:
#             raf_trees[chr_] = IntervalTree()
#     return raf_trees


# def unique_region_intervaltrees(unique_regions, chromosomes) -> dict[str,IntervalTree]:
#     unique_regions_trees = {chr:IntervalTree() for chr in chromosomes}
#     with open(unique_regions,'r') as f:
#         reader = csv.reader(f,delimiter='\t',quotechar='"')
#         for i,(chr,start,end) in tqdm(enumerate(reader)):
#             unique_regions_trees[chr][int(start):int(end)] = i
#     return unique_regions_trees


# # TODO: needs testing
# def get_valid_intervals(
#         seed_interval:tuple[int,int],
#         unique_regions:set[Interval],
#         max_extents:tuple[int,int],
#         min_cr_size:int,
#         min_anchor_size:int) -> list[Interval]:
#     """"computes a set of valid sub-intervals from unique regions that \
# are within the max_extens & don't overlap with the seed start-end-interval & \
# have a minimum anchor size"""
#     # compute a set of intervals that match criteria
#     # first, create an intervaltree of the unique regions
#     it_unique_regions = IntervalTree()
#     for ur in unique_regions:
#         it_unique_regions[ur[0]:ur[1]] = ur[2]
#     # chop all before and after maximum extents
#     it_unique_regions.chop(0,max_extents[0])
#     it_unique_regions.chop(max_extents[1],2**31-1)
#     # chop the core region
#     core_interval = (seed_interval[0]-min_cr_size//2,seed_interval[1]+min_cr_size//2)
#     it_unique_regions.chop(core_interval[0],core_interval[1])
#     # from all unique regions that are before the seed region, subtract min_anchor_size from end
#     # - get all unique regions that are before the core region
#     # - get the intervals that are to be subtracted
#     # - chop from it_unique_regions each last min_anchor_size bp
#     for interval in it_unique_regions[0:core_interval[0]]:
#         it_unique_regions.chop(interval[1]-min_anchor_size,interval[1])
#     # from all unique regions that are after the core region, add min_anchor_size to start
#     for interval in it_unique_regions[core_interval[1]:]:
#         it_unique_regions.chop(interval[0],interval[0]+min_anchor_size)
#     return sorted(it_unique_regions)


# def find_optimal_extent(crseed:datatypes.CrSeed, optimal_cr_size:int) -> tuple[int,int]:
#     left_optimum = find_crseed_closest_valid_point(crseed=crseed,target=crseed.start-optimal_cr_size//2)
#     right_optimum = find_crseed_closest_valid_point(crseed=crseed,target=crseed.end+optimal_cr_size//2)
#     return (left_optimum,right_optimum)


# def find_crseed_closest_valid_point(
#             crseed:datatypes.CrSeed,
#             target:int) -> int|None:
#     # target is a point on the reference
#     # min and max expansion are total bp from the start and end of CrSeed
#     # find closest point that is in a valid_interval of crseed
#     # if no such point exists, return None
#     # 1) check if the point is < or > crseed.start or crseed.end. If it is within crseed.start to crseed.end, then raise error
#     if crseed.start < target < crseed.end:
#         raise ValueError(f"target {target} is within crseed {crseed.start} to {crseed.end}")
#     valid_intervaltree = IntervalTree(crseed.valid_intervals)
#     # if target is left of crseed.start, then chop the right half of the valid_intervaltree
#     if target <= crseed.start:
#         valid_intervaltree.chop(crseed.start,2**31-1)
#     else: # target overlapping with (crseed.start:crseed.end) is already handled
#         valid_intervaltree.chop(0,crseed.end)
#     overlapping_valid_intervals:set[Interval] = valid_intervaltree[target]
#     # if target is on one oft the valid intervals, then set expansion to the distance of the target to the closest end of the interval
#     if len(overlapping_valid_intervals) > 0:
#         return target
#     else: # no direct overlap, find the closest valid point
#         sorted_intervals = sorted(valid_intervaltree, key=lambda x: min(abs(target-x.end),abs(x.begin-target)))
#         # the first interval is the closest. pick either start or end of the interval, depending on what is closer
#         if len(sorted_intervals) == 0:
#             # return either start or end of crseed, depending on what is closer to target
#             if abs(target-crseed.start) < abs(target-crseed.end):
#                 return crseed.start-50
#             else:
#                 return crseed.end+50
#         closest_interval = sorted_intervals[0]
#         return closest_interval.begin \
#             if abs(target-closest_interval.begin) < abs(target-closest_interval.end) \
#             else closest_interval.end


# def merge_CrSeeds(a:datatypes.CrSeed,b:datatypes.CrSeed) -> datatypes.CrSeed:
#     if a.chr != b.chr:
#         raise(f"cannot merge CrSeeds of different chromosomes: {a.chr} and {b.chr}")
#     if None in a.max_extents or None in b.max_extents:
#         raise ValueError(f"max_extents of CrSeed is None: {a.max_extents} or {b.max_extents}, but all must be integers.")
#     # check if all values in max_extents for both a and b are 0<=x<=2**31-1
#     if not all([0 <= x <= 2**31-1 for x in a.max_extents]) or not all([0 <= x <= 2**31-1 for x in b.max_extents]):
#         raise ValueError(f"max_extents of CrSeed is not within 0 to 2**31-1: {a.max_extents} or {b.max_extents}")
#     # transform to intervaltree and merge, then chop new start,end.
#     new_start = min(a.start,b.start)
#     new_end = max(a.end,b.end)
#     new_valid_intervals_it = IntervalTree(a.valid_intervals+b.valid_intervals)
#     new_valid_intervals_it.chop(new_start,new_end)
#     new_valid_intervals = sorted(new_valid_intervals_it)
#     return datatypes.CrSeed(
#         chr = a.chr,
#         chrID = a.chrID,
#         start = min(a.start,b.start),
#         end = max(a.end,b.end),
#         readIDs=set(a.readIDs).union(b.readIDs),
#         max_extents=(min(a.max_extents[0],b.max_extents[0]),max(a.max_extents[1],b.max_extents[1])),
#         actual_extents=(min(a.actual_extents[0],b.actual_extents[0]), max(a.actual_extents[1],b.actual_extents[1])),
#         signals=a.signals+b.signals,
#         valid_intervals=new_valid_intervals)

# # test if two crseeds need to be merged
# # def crSeed_overlapping_cores(
# #         cra:datatypes.CrSeed,
# #         crb:datatypes.CrSeed,
# #         min_cr_size:int) -> bool:
# #     half_min_cr_size = min_cr_size//2
# #     a:Interval = Interval(cra.start-half_min_cr_size,cra.end+half_min_cr_size)
# #     b:Interval = Interval(crb.start-half_min_cr_size,crb.end+half_min_cr_size)
# #     return a.overlaps(b)

# # test if two crseeds need to be separated
# def crSeed_overlapping_actual_extents(
#         cra:datatypes.CrSeed,
#         crb:datatypes.CrSeed) -> bool:
#     a:Interval = Interval(*cra.actual_extents)
#     b:Interval = Interval(*crb.actual_extents)
#     return a.overlaps(b)

# test_a = list(zip([10,20,50,70,90],[10,20,50,70,90]))
# test_b = list(zip([15,30,45,70,90,100],[10,20,40,65,80,95]))


# # used in find_CrSeed_separation
# def find_points_separation(
#             points_a:list[tuple[int,int]],
#             points_b:list[tuple[int,int]],
#             dist_func:callable=lambda x: x**2) -> tuple[tuple[int,int,int],tuple[int,int,int]]:
#     # build triples of points of form (position,distance,owner(0|1))
#     points_a = [(*x,0) for x in points_a]
#     points_b = [(*x,1) for x in points_b]
#     points = sorted(points_a+points_b,key=lambda x: x[0])
#     # find a pair x,y of adjacend elements in points, where
#     # x[0] < y[0] and x[2] != y[2] and x[1] + y[1] is minimal
#     # report both x[0] and y[0]
#     min_distance = 2**31-1
#     min_distance_points = None
#     for i in range(len(points)-1):
#         if points[i][2] == 0 and points[i+1][2] == 1:
#             distance = dist_func(points[i][1]) + dist_func(points[i+1][1])
#             if distance < min_distance:
#                 min_distance = distance
#                 min_distance_points = (points[i],points[i+1])
#     return min_distance_points


# # used in find_CrSeed_separation
# def intervals_to_points(
#         intervals:list[Interval],
#         point_of_distance:int,
#         stepsize:int) -> list[tuple[int,int]]:
#     points = []
#     for interval in intervals:
#         points.append((interval.begin,abs(interval.begin-point_of_distance)))
#         points.append((interval.end,abs(interval.end-point_of_distance)))
#         # add a point every stepsize (e.g. 100bp) for each interval
#         for point in range(interval.begin+stepsize,interval.end,stepsize):
#             points.append((point,abs(point-point_of_distance)))
#     return points


# def get_valid_intervals_of_a_in_a_to_b(
#         a:datatypes.CrSeed,
#         b:datatypes.CrSeed) -> list[Interval]:
#     a_intervalTree = IntervalTree(a.valid_intervals)
#     a,b = (a,b) if a.start < b.start else (b,a)
#     a_intervalTree.chop(0,a.min_extents[1])
#     a_intervalTree.chop(b.min_extents[0],2**31-1)
#     return sorted(a_intervalTree)


# def find_CrSeed_separation(
#         cra:datatypes.CrSeed,
#         crb:datatypes.CrSeed) -> tuple[int,int]|None:
#     # find two points (that can be the same) that are between cra and crb (start-end interval)
#     # both points must be in the valid_intervals of both cra and crb
#     cra,crb = (cra,crb) if cra.start < crb.start else (crb,cra)
#     cra_valid_intervals_between_cra_and_crb = get_valid_intervals_of_a_in_a_to_b(a=cra,b=crb)
#     crb_valid_intervals_between_cra_and_crb = get_valid_intervals_of_a_in_a_to_b(a=crb,b=cra)
#     cra_valid_intervals_tree = IntervalTree(cra_valid_intervals_between_cra_and_crb)
#     crb_valid_intervals_tree = IntervalTree(crb_valid_intervals_between_cra_and_crb)
#     # check if the point between cra and crb is in both valid intervals
#     midpoint = (cra.end+crb.start)//2
#     if cra_valid_intervals_tree[midpoint] and crb_valid_intervals_tree[midpoint]:
#         return (midpoint,midpoint)
#     # if not, find the two closest points that are in both valid intervals
#     stepsize_a = max(100,
#         int(sum([interval.end - interval.begin
#             for interval in cra_valid_intervals_between_cra_and_crb])/100))
#     points_a = intervals_to_points(
#         intervals=cra_valid_intervals_between_cra_and_crb,
#         point_of_distance=cra.end,
#         stepsize=stepsize_a)
#     stepsize_b = max(100,
#         int(sum([interval.end - interval.begin
#             for interval in crb_valid_intervals_between_cra_and_crb])/100))
#     points_b = intervals_to_points(
#         intervals=crb_valid_intervals_between_cra_and_crb,
#         point_of_distance=crb.start,
#         stepsize=stepsize_b)
#     # all intervals between cra and crb
#     points_separation = find_points_separation(points_a=points_a,points_b=points_b)
#     if points_separation is None:
#         return None
#     else:
#         return (points_separation[0][0], points_separation[1][0])


# def compute_min_extents(crseed:datatypes.CrSeed) -> tuple[int,int]:
#     min_start = crseed.max_extents[0]
#     min_end = crseed.max_extents[1]
#     # find the next point to min_start
#     # crseed has valid_intervals, which are in the two intervals
#     #   crseed.max_extents[0] <= x <= crseed.start
#     #   crseed.end <= x <= crseed.max_extents[1]
#     # thus, the valid_interval x closest to crseed.start that is x <= crseed.start
#     # is the left min_extent
#     # and the valid_interval y closest to crseed.end that is x >= crseed.end
#     # is the right min_extent
#     # 1) find the interval that is closest to crseed.start
#     for interval in sorted(crseed.valid_intervals,key=lambda x: abs(crseed.start - x.end)):
#         if interval.end <= crseed.start:
#             min_start = interval.end
#             break
#     # 2) find the interval that is closest to crseed.end
#     for interval in sorted(crseed.valid_intervals,key=lambda x: abs(x.begin - crseed.end)):
#         if interval.begin >= crseed.end:
#             min_end = interval.begin
#             break
#     return (min_start,min_end)


# def crSeed_to_candidate_region(crseed:datatypes.CrSeed,crID:int) -> datatypes.CandidateRegion:
#     # check type of crseed
#     if not isinstance(crseed,datatypes.CrSeed):
#         raise ValueError(f"crseed is not a CrSeed, but a {type(crseed)}")
#     sv_signals:list[datatypes.ExtendedSVsignal] = crseed.signals
#     # for sv_signal in sv_signals:
#     #     sv_signal.readname = sv_signal.readname.split(".")[0]
#     return datatypes.CandidateRegion(
#         crID=crID,
#         chr=crseed.chr,
#         referenceID=crseed.chrID,
#         referenceStart=crseed.actual_extents[0],
#         referenceEnd=crseed.actual_extents[1],
#         sv_signals=sv_signals)


# def create_containers_db(path_db:Path,timeout:float):
#     assert timeout > 0.0, f"timeout must be > 0.0. It is {timeout}"
#     assert type(path_db) == PosixPath or type(path_db) == str, f"path_db must be a Path or str. It is {type(path_db)}"

#     if path_db.exists():
#         log.warning(f"Database {path_db} exists. Overwriting it.")
#         path_db.unlink()

#     conn = sqlite3.connect(path_db)
#     curs = conn.cursor()
#     curs.execute("""CREATE TABLE IF NOT EXISTS containers
#         (crID INTEGER PRIMARY KEY, data TEXT)""")
#     conn.execute(f'pragma busy_timeout={str(int(timeout*1000))}')
#     conn.commit()
#     conn.close()

# def write_containers_to_db(path_db:Path,crs_containers:list[dict[str,object]]):
#     assert type(path_db) == PosixPath or type(path_db) == str, f"path_db must be a Path or str. It is {type(path_db)}"
#     assert path_db.exists(), f"Database {path_db} does not exist"
#     assert path_db.is_file(), f"Database {path_db} is not a file"
#     assert type(crs_containers) == list, f"crs_containers is not a list"
#     assert all([type(crs_container) == dict for crs_container in crs_containers]), f"crs_containers is not a list of dicts"

#     conn = sqlite3.connect(path_db)
#     curs = conn.cursor()
#     pairs = [(int(crs_container['crID']),serialize_crs_container(crs_container)) for crs_container in crs_containers]
#     curs.executemany("INSERT INTO containers VALUES (?,?)",pairs)
#     conn.commit()
#     conn.close()


def create_crs_db(path_db: Path, timeout: float) -> None:
    # create sqlite3 database with primary key crID and pickled json.dumps(CandidateRegion.unstructure()) as vales named candidate_region
    assert timeout > 0.0, f"timeout must be > 0.0. It is {timeout}"
    assert (
        type(path_db) == PosixPath or type(path_db) == str
    ), f"path_db must be a Path or str. It is {type(path_db)}"

    if path_db.exists():
        logger.warning(f"Database {path_db} exists. Overwriting it.")
        path_db.unlink()

    conn = sqlite3.connect(path_db)
    c = conn.cursor()
    c.execute(
        """CREATE TABLE IF NOT EXISTS candidate_regions
        (crID INTEGER PRIMARY KEY, candidate_region TEXT)"""
    )
    c.execute(f"pragma busy_timeout={str(int(timeout*1000))}")
    conn.commit()
    c.close()
    conn.close()


def write_crs_to_db(path_db: Path, crs: list[datatypes.CandidateRegion]) -> None:
    assert (
        type(path_db) == PosixPath or type(path_db) == str
    ), f"path_db must be a Path or str. It is {type(path_db)}"
    assert path_db.exists(), f"Database {path_db} does not exist"
    assert path_db.is_file(), f"Database {path_db} is not a file"
    assert type(crs) == list, "crs is not a list"
    assert all(
        [type(cr) == datatypes.CandidateRegion for cr in crs]
    ), "crs is not a list of datatypes.CandidateRegion"

    conn = sqlite3.connect(path_db)
    c = conn.cursor()
    pairs = [[cr.crID, pickle.dumps(cr.unstructure())] for cr in crs]
    c.executemany("INSERT INTO candidate_regions VALUES (?,?)", pairs)
    conn.commit()
    c.close()
    conn.close()


def load_crs_from_db(
    path_db: Path, crIDs: list[int] | None = None
) -> list[datatypes.CandidateRegion]:
    if crIDs:
        assert type(crIDs) == list, "crIDs is not a list"
        assert all(
            [type(crID) == int for crID in crIDs]
        ), "crIDs is not a list of integers"
    assert (
        type(path_db) == PosixPath or type(path_db) == str
    ), f"path_db must be a Path or str. It is {type(path_db)}"
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


# def crSeeds_to_candidate_regions(
#         output:Path,
#         threads:int,
#         crseeds:list[datatypes.CrSeed],
#         tmp_dir_path:Path|None=None) -> typing.Generator[datatypes.CandidateRegion,None,None]:
#     for crseed in crseeds:
#         cr:datatypes.CandidateRegion = crSeed_to_candidate_region(crseed=crseed,crID=0)
#         yield cr


#     # write all crs to file
#     tmp_output = tempfile.NamedTemporaryFile(suffix='.tmp.crs.tsv',dir=tmp_dir_path,delete=False if tmp_dir_path else True)
#     with open(tmp_output.name,'w') as f:
#         writer = csv.writer(f,delimiter='\t',quotechar='"')
#         for crseed in crseeds:
#             cr:datatypes.CandidateRegion = crSeed_to_candidate_region(crseed=crseed,crID=0)
#             writer.writerow([cr.referenceID,cr.referenceStart,cr.referenceEnd,json.dumps(cr.unstructure())])
#     # compress and index
#     tmp_out_re_index = tempfile.NamedTemporaryFile(suffix='.tmp.crs.tsv.reindexed',dir=tmp_dir_path,delete=False if tmp_dir_path else True)
#     alignments_to_rafs.compress_and_index_bedlike(
#         input=tmp_output.name,
#         output=tmp_out_re_index.name,
#         threads=threads,
#         sort_numerically=True)
#     with open(output,'w') as f:
#         writer = csv.writer(f,delimiter='\t',quotechar='"')
#         for crID,cr in enumerate(util.yield_from_crs(input=tmp_out_re_index.name)):
#             cr.crID = crID
#             writer.writerow([cr.referenceID,cr.referenceStart,cr.referenceEnd,json.dumps(cr.unstructure())])
#     # read candidate regions again and give new crIDs accoring to the order in the file


# def actual_crseed_extents(
#             min_cr_size:int,
#             optimal_cr_size:int,
#             tandem_repeat_tree:IntervalTree,
#             tree:IntervalTree) -> IntervalTree:
#     dropped_intervals = IntervalTree()
#     seeds = []
#     for interval in sorted(tree):
#         cr:datatypes.CrSeed = interval.data
#         # extend the seed to min_anchor_size
#         cr.

#         cr.valid_intervals = get_valid_intervals(
#                 seed_interval=(interval.begin,interval.end),
#                 unique_regions=unique_region_tree[cr.max_extents[0]:cr.max_extents[1]],
#                 max_extents=cr.max_extents,
#                 min_cr_size=min_cr_size,
#                 min_anchor_size=300)
#             # if a cr has no valid intervals, left and right, skip it
#         if len(cr.valid_intervals) == 0:
#             logger.debug(f"removed CrSeed on interval {cr.chr}:{cr.start}-{cr.end} because it has no valid intervals left or right.")
#             dropped_intervals[interval.begin:interval.end] = interval.data
#             tree.remove(interval)
#             continue
#             # check if valid intervals are present left and righht of cr.start,cr.end
#         valid_intervals_tree = IntervalTree(cr.valid_intervals)
#         if len(valid_intervals_tree[0:cr.start]) == 0 and len(valid_intervals_tree[cr.end:]) == 0:
#             logger.debug(f"removed CrSeed on interval {cr.chr}:{cr.start}-{cr.end} because it has no valid intervals left and right.")
#             dropped_intervals[interval.begin:interval.end] = cr
#             interval.data.actual_extents = (cr.start,cr.end)
#             tree.remove(interval)
#             continue
#         actual_extents = (find_crseed_closest_valid_point(crseed=cr,target=interval.begin-int(optimal_cr_size/2)),
#                               find_crseed_closest_valid_point(crseed=cr,target=interval.end+int(optimal_cr_size/2)))
#         if actual_extents[0] is None or actual_extents[1] is None:
#             raise ValueError(f"expansion_left or expansion_right is None, but should be integers.")
#         cr.actual_extents = actual_extents
#     return tree, dropped_intervals

# # test this!
# def separate_actual_extents(tree:IntervalTree) -> list[datatypes.CrSeed]:
#     if len(tree) == 0:
#         return []
#     if len(tree) == 1:
#         return [next(tree.iter()).data]
#     final_crseeds:list[datatypes.CrSeed] = []
#     crseeds:list[datatypes.CrSeed] = sorted([x.data for x in tree],key=lambda x: (x.min_extents))
#     for i in range(len(crseeds)-1):
#         cra = crseeds[i]
#         crb = crseeds[i+1]
#         # check if both overlap
#         if cra.chrID == crb.chrID and Interval(*cra.actual_extents).overlaps(Interval(*crb.actual_extents)):
#                 # adjust the actual extents of cra and crb so that they don't overlap anymore.
#             separation:tuple[int,int]|None = find_CrSeed_separation(cra=cra,crb=crb)
#             if separation is None:
#                     # no good separation could be found. merge both crseeds and make crseeds[i+1] the new one
#                 merged_cr = merge_CrSeeds(a=cra,b=crb)
#                 crseeds[i+1] = merged_cr
#             else:
#                 cra.actual_extents = (cra.actual_extents[0],separation[0])
#                 crb.actual_extents = (separation[1],crb.actual_extents[1])
#                 # check if cra or crb have illegitimate actual extents
#                 if cra.actual_extents[0] >= cra.actual_extents[1] or crb.actual_extents[0] >= crb.actual_extents[1]:
#                     import pdb; pdb.set_trace()
#                     # only add cra if it was not merged
#                 final_crseeds.append(cra)
#         else:
#             final_crseeds.append(cra)
#     final_crseeds.append(crseeds[-1])
#     return final_crseeds

# def crSeeds_add_max_extents(tree, rafs_tree) -> None:
#     for interval in tqdm(tree):
#             # in-place processing of CrSeeds
#         cr:datatypes.CrSeed = interval.data
#         max_extents = (interval.begin,interval.end)
#         for raf_interval in rafs_tree[interval.begin:interval.end]:
#             max_extents = (min(max_extents[0],raf_interval.begin),max(max_extents[1],raf_interval.end))
#         cr.max_extents = max_extents


# def crSeeds_from_cr_proto_regions(chr, tree) -> IntervalTree:
#     cr_seeds_tree = IntervalTree()
#     for interval in tree:
#         local_svsignals:list[datatypes.ExtendedSVsignal] = interval.data
#         cr_seeds_tree[interval.begin:interval.end] = datatypes.CrSeed(
#                 chr=chr,
#                 chrID=local_svsignals[0].chrID,
#                 start=interval.begin,
#                 end=interval.end,
#                 readIDs=list(set([x.readname for x in local_svsignals])),
#                 signals=local_svsignals)
#     return cr_seeds_tree


# def process_cr_seeds(
#         chromosome:str,
#         cr_proto_regions_tree:IntervalTree,
#         tandem_repeat_tree:IntervalTree,
#         min_cr_size:int,
#         optimal_cr_size:int) -> list[datatypes.CandidateRegion]:
#     # ==== merge cr_proto_regions ==== #
#     # merge down the cr_proto_regions_tree
#     # each interval in cr_proto_regions_tree has a list of svsignals
#     # merge overlapping intervals and concatenate the lists of svsignals in the data reducer
#     # extend each cr_proto_region to the max extents of overlapping tandem repets
#     # extend each cr_proto_region to have a minimal size of min_cr_size
#     # finally, create candidate regions from lists of datatypes.ExtendedSVsignal (the merged cr_proto_regions_tree intervals)


def process_cr_seeds(
    chromosome: str,
    cr_proto_regions_tree: IntervalTree,
    tandem_repeat_tree: IntervalTree,
    min_cr_size: int,
) -> list[datatypes.CandidateRegion]:
    """
    Process candidate region seeds into final CandidateRegion objects.

    Steps:
    1. Merge overlapping intervals in cr_proto_regions_tree
    2. Concatenate their sv_signal lists
    3. Extend to tandem repeat boundaries
    4. Ensure minimum size
    5. Create CandidateRegion objects
    """
    if len(cr_proto_regions_tree) == 0:
        return []

    # Step 1: Merge overlapping intervals and concatenate sv_signals
    cr_proto_regions_tree.merge_overlaps(
        data_reducer=lambda x, y: x + y  # concatenate lists of ExtendedSVsignal
    )

    # Step 2: Extend by overlapping tandem repeats
    # Collect new intervals with extended boundaries
    extended_intervals = []
    for interval in cr_proto_regions_tree:
        start = interval.begin
        end = interval.end
        sv_signals = interval.data

        # Find overlapping tandem repeats and extend boundaries
        overlapping_repeats = tandem_repeat_tree[start:end]
        if overlapping_repeats:
            for repeat in overlapping_repeats:
                start = min(start, repeat.begin)
                end = max(end, repeat.end)

        extended_intervals.append((start, end, sv_signals))

    # Rebuild the tree with extended intervals
    cr_proto_regions_tree = IntervalTree.from_tuples(extended_intervals)

    # Step 3: Extend to minimum size
    size_extended_intervals = []
    for interval in cr_proto_regions_tree:
        start = interval.begin
        end = interval.end
        sv_signals = interval.data

        current_size = end - start
        if current_size < min_cr_size:
            extension = (min_cr_size - current_size) // 2
            start = max(0, start - extension)
            end += extension
            # If min_cr_size is odd, extend the end by 1 more
            if (min_cr_size - current_size) % 2 != 0:
                end += 1

        size_extended_intervals.append((start, end, sv_signals))

    # Rebuild tree again
    cr_proto_regions_tree = IntervalTree.from_tuples(size_extended_intervals)

    # Step 4: Merge again after extensions
    cr_proto_regions_tree.merge_overlaps(
        data_reducer=lambda x, y: x + y  # concatenate lists of ExtendedSVsignal
    )

    # Step 5: Create CandidateRegion objects
    candidate_regions = []
    for i, interval in enumerate(sorted(cr_proto_regions_tree)):
        sv_signals = interval.data  # list of ExtendedSVsignal
        cr = datatypes.CandidateRegion(
            crID=i,  # Placeholder, should be set later
            chr=chromosome,
            referenceID=sv_signals[0].chrID,
            referenceStart=interval.begin,
            referenceEnd=interval.end,
            sv_signals=sv_signals,
        )
        candidate_regions.append(cr)

    return candidate_regions


def mp_process_cr_seeds(args: dict) -> list[datatypes.CandidateRegion]:
    return process_cr_seeds(**args)


def tandem_repeat_intervaltrees(
    tandem_repeats: Path, chromosomes: list[str]
) -> dict[str, IntervalTree]:
    """Load tandem repeats from bed file into interval trees per chromosome"""
    tandem_repeat_trees = {c: IntervalTree() for c in chromosomes}

    with open(tandem_repeats, "r") as f:
        reader = csv.reader(f, delimiter="\t", quotechar='"')
        previous_chromosome: str | None = None
        cache = []

        for row in reader:
            chr = str(row[0])
            if chr not in tandem_repeat_trees:
                continue

            if previous_chromosome == chr or previous_chromosome is None:
                cache.append((int(row[1]), int(row[2])))
                previous_chromosome = chr
            else:
                # Construct the interval tree for the previous chromosome
                if previous_chromosome and cache:
                    tandem_repeat_trees[previous_chromosome] = IntervalTree.from_tuples(
                        cache
                    )
                cache = [(int(row[1]), int(row[2]))]
                previous_chromosome = chr

        # Don't forget the last chromosome!
        if previous_chromosome and cache:
            tandem_repeat_trees[previous_chromosome] = IntervalTree.from_tuples(cache)

    return tandem_repeat_trees


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
    dropped: Path | None = None,
    bedgraph: Path | None = None,
    tmp_dir_path: Path | None = None,
) -> None:
    csv.field_size_limit(sys.maxsize)

    # Load chromosome names
    chromosomes = list(util.create_ref_dict(reference).values())
    logger.info("Loading tandem repeats to interval trees...")
    tandem_repeat_trees = tandem_repeat_intervaltrees(
        tandem_repeats=tandem_repeats, chromosomes=chromosomes
    )

    logger.info("Parsing signalstrengths...")
    svsignals: list[datatypes.ExtendedSVsignal] = list(
        util.yield_from_extendedSVsignal(input=signalstrengths)
    )

    # Compute signal statistics
    values_signal = np.array([x.strength for x in svsignals], dtype=float)
    values_depth = np.array([x.coverage for x in svsignals], dtype=int)

    # Handle division by zero
    with np.errstate(divide="ignore", invalid="ignore"):
        signal_normalized = np.divide(values_signal, values_depth)
        signal_normalized[~np.isfinite(signal_normalized)] = 0.0

    mask_absolute_signal = values_signal > filter_absolute
    mask_normalized_signal = signal_normalized > filter_normalized
    mask_final = mask_absolute_signal & mask_normalized_signal

    logger.info(f"Signals passing filters: {np.sum(mask_final)}/{len(svsignals)}")

    # Optional bedgraph output
    if bedgraph is not None:
        write_to_bedgraphs(
            signals=svsignals,
            values_signal=values_signal,
            values_normalized=signal_normalized,
            values_masked=mask_final.astype(float),
            bedgraph=bedgraph,
        )

    # Create proto regions per chromosome
    logger.info("Creating cr_proto_regions (batched)...")
    tmp_per_chr: dict[str, list[tuple[int, int, list[datatypes.ExtendedSVsignal]]]] = {
        c: [] for c in chromosomes
    }

    idxs = np.flatnonzero(mask_final)
    for i in tqdm(idxs, total=idxs.size, desc="Building proto regions"):
        s = svsignals[i]
        c = s.chr
        if c not in tmp_per_chr:
            logger.warning(f"Signal chromosome {c} not in reference, skipping")
            continue
        start = max(0, s.ref_start - buffer_region_radius)
        end = s.ref_end + buffer_region_radius
        tmp_per_chr[c].append((start, end, [s]))

    # Build interval trees
    cr_proto_regions_trees = {}
    for c, tups in tmp_per_chr.items():
        cr_proto_regions_trees[c] = (
            IntervalTree.from_tuples(tups) if tups else IntervalTree()
        )

    # Prepare multiprocessing jobs
    logger.info("Creating jobs for parallel processing of cr_seeds...")
    jobs_cr_seeds = [
        {
            "chromosome": chr,
            "cr_proto_regions_tree": tree,
            "tandem_repeat_tree": tandem_repeat_trees[chr],
            "min_cr_size": min_cr_size,
        }
        for chr, tree in cr_proto_regions_trees.items()
        if len(tree) > 0  # Skip empty chromosomes
    ]

    logger.info(
        f"Processing {len(jobs_cr_seeds)} chromosomes with {threads} threads..."
    )

    # Process in parallel - results is a list of lists
    results: list[list[datatypes.CandidateRegion]] = []
    with mp.Pool(min(threads, len(jobs_cr_seeds))) as pool:
        # Use imap for progress tracking
        for result in tqdm(
            pool.imap(mp_process_cr_seeds, jobs_cr_seeds),
            total=len(jobs_cr_seeds),
            desc="Processing chromosomes",
        ):
            results.extend(result)

    logger.info(f"Generated {len(results)} candidate regions")

    # Filter by read count
    if len(results) > 0:
        read_counts = [len(cr.get_read_names()) for cr in results]
        median_read_count = np.median(read_counts)
        max_read_count = 4 * median_read_count

        logger.info(
            f"Median read count: {median_read_count:.1f}, "
            f"max threshold: {max_read_count:.1f}"
        )

        filtered_results = []
        filtered_regions = []

        for cr in results:
            read_count = len(cr.get_read_names())
            if read_count <= max_read_count:
                filtered_results.append(cr)
            else:
                filtered_regions.append(cr)

        logger.info(
            f"Filtered {len(filtered_regions)} candidate regions with "
            f"excessive read counts (>{max_read_count:.1f})"
        )

        results = filtered_results

    # Reassign sequential crIDs
    for i, cr in enumerate(results):
        cr.crID = i

    # Write to database
    logger.info(f"Writing {len(results)} candidate regions to database...")
    create_crs_db(path_db=output, timeout=100.0)
    if len(results) > 0:
        write_crs_to_db(path_db=output, crs=results)

    # Write dropped regions if requested
    if dropped is not None and len(filtered_regions) > 0:
        # assign sequential crIDs to filtered regions
        for i, cr in enumerate(filtered_regions):
            cr.crID = i

        logger.info(f"Writing {len(filtered_regions)} dropped regions to database...")
        create_crs_db(path_db=dropped, timeout=100.0)
        write_crs_to_db(path_db=dropped, crs=filtered_regions)

    logger.info(f"Done. Created {len(results)} candidate regions.")


def run(args, **kwargs):
    create_candidate_regions(
        signalstrengths=args.input,
        reference=args.reference,
        output=args.output,
        tandem_repeats=args.tandem_repeats,
        buffer_region_radius=args.buffer_region_radius,
        bedgraph=args.bedgraph,
        filter_absolute=args.filter_absolute,
        filter_normalized=args.filter_normalized,
        threads=args.threads,
        min_cr_size=args.min_cr_size,
        dropped=args.dropped,
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
        help="Path to reference in uncompressed .fa(sta) format. Only the fasta index (.fai) is really required to exist.",
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
        "--buffer_region_radius",
        type=int,
        required=False,
        default=300,
        help="Radius of the initial seeds' sizes around signal that can constitute a candidate region seed.",
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
        default=500,
        help="Minimum size of candidate region.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
