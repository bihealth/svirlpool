from libcpp.vector cimport vector
from cython.operator cimport dereference as deref
from typing import List

def compute_complexity(int wsize, kmers: List[int], dna: str) -> List[float]:
    """
    Compute the sequence complexity.

    Parameters:
    wsize (int): Window size.
    kmers (List[int]): List of k-mer values.
    dna (str): DNA sequence.

    Returns:
    List[float]: Computed complexity values.
    """
    cdef vector[uint8_t] c_kmers = kmers  # Convert Python list to C++ vector
    cdef vector[float] result = run_program(wsize, c_kmers, dna)  # Call C++ function
    return [deref(val) for val in result]  # Convert C++ vector to Python list
