
import numpy as np
import pysam

from svirlpool.util.util import Direction, get_read_position_on_ref


def test_deletion_none_direction():
    # Create a dummy alignment with a deletion
    # Ref:  0123456789
    #       MMDDMM
    # Read: MM--MM
    # CIGAR: 2M2D2M
    
    a = pysam.AlignedSegment()
    a.query_name = "test"
    a.query_sequence = "ACGT" # 4 bases
    a.reference_start = 100
    a.cigartuples = [(0, 2), (2, 2), (0, 2)]
    
    # Read coordinates:
    # Block 0 (2M): starts 0, ends 2. Ref: 100-102.
    # Block 1 (2D): starts 2, ends 2. Ref: 102-104.
    # Block 2 (2M): starts 2, ends 4. Ref: 104-106.
    
    # If we query position 2.
    # It matches Block 1 (Deletion) because it's the first block where 2 <= 2 <= 2.
    
    print("Testing position 2 (at deletion) with Direction.NONE")
    try:
        pos = get_read_position_on_ref(a, 2, Direction.NONE)
        print(f"Result: {pos}")
    except Exception as e:
        print(f"Error: {e}")

    # Also check what happens if we use Direction.LEFT or RIGHT
    print("Testing position 2 with Direction.LEFT")
    pos_left = get_read_position_on_ref(a, 2, Direction.LEFT)
    print(f"Result LEFT: {pos_left}")

    print("Testing position 2 with Direction.RIGHT")
    pos_right = get_read_position_on_ref(a, 2, Direction.RIGHT)
    print(f"Result RIGHT: {pos_right}")

if __name__ == "__main__":
    test_deletion_none_direction()
