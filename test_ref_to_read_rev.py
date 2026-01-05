
import pysam
import numpy as np
from svirlpool.util.util import get_ref_position_on_read, Direction

def test_insertion_direction_reverse():
    # Create a dummy alignment with an insertion, REVERSE strand
    # Ref:  0123456789
    #       MM--MM
    # Read: MMIIMM
    # CIGAR: 2M2I2M
    
    a = pysam.AlignedSegment()
    a.query_name = "test_ins_rev"
    a.query_sequence = "ACGTAC" # 6 bases
    a.reference_start = 100
    a.cigartuples = [(0, 2), (1, 2), (0, 2)]
    a.is_reverse = True
    
    # Ref coordinates:
    # Block 0 (2M): starts 100, ends 102. 
    # Block 1 (2I): starts 102, ends 102.
    # Block 2 (2M): starts 102, ends 104.
    
    # Read coordinates (Original Fragment):
    # If is_reverse, get_starts_ends flips coordinates.
    # Length = 6.
    # Block 0 (2M): Read BAM 0-2. Original: 6-2 = 4 to 6-0 = 6? (4,5)
    # Block 1 (2I): Read BAM 2-4. Original: 6-4 = 2 to 6-2 = 4? (2,3)
    # Block 2 (2M): Read BAM 4-6. Original: 6-6 = 0 to 6-4 = 2? (0,1)
    
    # So:
    # Block 0 (Ref 100-102) -> Read Orig 4-6
    # Block 1 (Ref 102-102) -> Read Orig 2-4
    # Block 2 (Ref 102-104) -> Read Orig 0-2
    
    # Note: Read Orig 0 is the 5' end of the original fragment.
    # If aligned to reverse strand, 5' end of original fragment aligns to High Ref Coord (104).
    # 3' end of original fragment aligns to Low Ref Coord (100).
    
    print("Testing position 102 (at insertion site) REVERSE")
    
    # Ref 102 is at the boundary of Block 0 and Block 2 (separated by Block 1).
    
    # LEFT on Ref (towards 100).
    # Should pick Block 0 (Ref 100-102).
    # Read pos should be 4 (start of Block 0 on Read Orig? or end?).
    # Block 0 corresponds to Read Orig 4-6.
    # Ref 102 is the END of Block 0.
    # So Read pos should be 4 (if 4-6 is 4,5).
    
    print("Testing position 102 with Direction.LEFT")
    try:
        pos_left = get_ref_position_on_read(a, 102, Direction.LEFT)
        print(f"Result LEFT: {pos_left}")
    except Exception as e:
        print(f"Error LEFT: {e}")

    # RIGHT on Ref (towards 104).
    # Should pick Block 2 (Ref 102-104).
    # Read pos should be 2 (start of Block 2 on Read Orig? or end?).
    # Block 2 corresponds to Read Orig 0-2.
    # Ref 102 is the START of Block 2.
    # So Read pos should be 2.
    
    print("Testing position 102 with Direction.RIGHT")
    try:
        pos_right = get_ref_position_on_read(a, 102, Direction.RIGHT)
        print(f"Result RIGHT: {pos_right}")
    except Exception as e:
        print(f"Error RIGHT: {e}")

if __name__ == "__main__":
    test_insertion_direction_reverse()
