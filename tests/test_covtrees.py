#%%
from pathlib import Path

import numpy as np
from intervaltree import Interval, IntervalTree

from svirlpool.util import covtree

#%%

# ================================================================
# ----------100           -> 1,0,100,r0
# 20----------120         -> 1,20,120,r1
#       80----------180   -> 1,80,180,r2
# -> (0,20,1), (20,80,2), (80,100,3), (100,120,2), (120,180,1)
# ================================================================

def test_parallel_coverage_computation__simple():
    all_positions = {'1': [0, 20, 80, 100, 120, 180]}
    intervall_trees = {'1': IntervalTree([Interval(0, 100, 'r0'), Interval(20, 120, 'r1'), Interval(80, 180, 'r2')])}
    num_workers=2
    result = covtree.parallel_coverage_computation(all_positions=all_positions, intervall_trees=intervall_trees, num_workers=num_workers)
    expected = {'1': IntervalTree([Interval(0, 20, 1), Interval(20, 80, 2), Interval(80, 100, 3), Interval(100, 120, 2), Interval(120, 180, 1)])}
    assert result == expected
    
def test_construct_interval_trees_simple():
    data = np.array([
        ['1', '0', '100', 'r0'],
        ['1', '20', '120', 'r1'],
        ['1', '80', '180', 'r2']])
    result = covtree.construct_interval_trees(data=data)
    expected = {'1': IntervalTree([Interval(0, 100, 'r0'), Interval(20, 120, 'r1'), Interval(80, 180, 'r2')])}, {'1': [0, 20, 80, 100, 120, 180]}
    assert result[0] == expected[0]
    assert result[1] == expected[1]
