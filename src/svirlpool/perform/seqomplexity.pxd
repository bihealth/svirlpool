from libcpp.vector cimport vector

cdef extern from "sequence_complexity.cpp":
    cdef vector[float] run_program(size_t wsize, vector[uint8_t] kmers, str dna)
