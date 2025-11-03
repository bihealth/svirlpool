#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "sequence_complexity.cpp"

namespace py = pybind11;

std::vector<float> compute_complexity(int wsize, std::vector<uint8_t> kmers, std::string dna) {
    return run_program(wsize, kmers, dna);
}

PYBIND11_MODULE(seqomplexity, m) {
    m.def("compute_complexity", &compute_complexity, "Compute sequence complexity");
}
