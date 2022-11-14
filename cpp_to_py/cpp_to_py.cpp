#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include "..\scripts_cpp\causets_cpp\functions.h"
#include "..\scripts_cpp\causets_cpp\vecfunctions.h"
#include "..\scripts_cpp\sprinkle.cpp"
/*
Produced module is in build/Release and can be imported as
import module_name inside a python file (with relevant).

To do it.
Go to build folder.
In the command line, write:

cmake ..
cmake --build . --target ALL_BUILD --config Release

*/

namespace py = pybind11;

int add(int i, int j) {
    return i + j;


}

PYBIND11_MODULE(causets_cppmodule, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");
    m.def("MM_drelation", &MM_drelation, "Calculates MM dimension");
    m.def("distinct_randint2", &distinct_randint2, "Create random array of size 'size'\
                                            with numbers between [0,N-1]");
    m.def("multiply", &multiply, "Multiplies two ints");
    m.def("get_coords", &get_coords, "gets coords");
}
