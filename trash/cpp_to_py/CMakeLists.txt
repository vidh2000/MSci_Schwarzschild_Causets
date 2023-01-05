cmake_minimum_required(VERSION 3.4...3.18)
project(example)
add_subdirectory(pybind11)
pybind11_add_module(causets_cppmodule cpp_to_py.cpp)