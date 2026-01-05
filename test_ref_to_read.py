
import pysam
import numpy as np
from svirlpool.util.util import get_ref_position_on_read, Direction

def test_insertion_direction():
    # Create a dummy alignment with an insertion
    # Ref:  0123456789
    #       MM--MM
    # Read: MMIIMM
    # CIGAR: 2M2I2M
    
    a = pysam.AlignedSegment()
    a.query_name = "test_ins"
    a.query_sequence = "ACGTAC" # 6 bases
    a.reference_start = 100
    a.cigartuples = [(0, 2), (1, 2), (0, 2)]
    
    # Ref coordinates:
    # Block 0 (2M): starts 100, ends 102. Read: 0-2.
    # Block 1 (2I): starts 102, ends 102. Read: 2-4. (Ref length 0)
    # Block 2 (2M): starts 102, ends 104. Read: 4-6.
    
    print("Testing position 102 (at insertion site)")
    
    # LEFT should pick Block 0 (ending at 102). Read pos should be 2.
    print("Testing position 102 with Direction.LEFT")
    try:
        pos_left = get_ref_position_on_read(a, 102, Direction.LEFT)
        print(f"Result LEFT: {pos_left}")
    except Exception as e:
        print(f"Error LEFT: {e}")

    # RIGHT should pick Block 2 (starting at 102). Read pos should be 4.
    print("Testing position 102 with Direction.RIGHT")
    try:
        pos_right = get_ref_position_on_read(a, 102, Direction.RIGHT)
        print(f"Result RIGHT: {pos_right}")
    except Exception as e:
        print(f"Error RIGHT: {e}")

    # NONE should pick Block 0 (first match). Read pos should be 2.
    print("Testing position 102 with Direction.NONE")
    try:
        pos_none = get_ref_position_on_read(a, 102, Direction.NONE)
        print(f"Result NONE: {pos_none}")
    except Exception as e:
        print(f"Error NONE: {e}")

if __name__ == "__main__":
    test_insertion_direction()
