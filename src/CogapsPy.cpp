#include "Rpackage/src/GapsRunner.h"
#include "Rpackage/src/math/Random.h"

#include <iostream>
#include <unistd.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

//int CoGAPS(const std::string &dataPath, unsigned numPatterns,
//unsigned maxIterations)
//{
//    return 1;
//}

int add(int i, int j)
{
    return i + j;
}

PYBIND11_MODULE(cogaps, m)
{
    m.doc() = "CoGAPS Python Package";
    //m.def("CoGAPS", &CoGAPS, "main function");
    m.def("add", &add, "A function which adds two numbers");
}