#include <pybind11/pybind11.h>
#include <iostream>

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

PYBIND11_MODULE(module_name, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");
}
