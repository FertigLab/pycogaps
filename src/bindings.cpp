#include "Rpackage/src/GapsRunner.h"
#include "Rpackage/src/utils/GlobalConfig.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int i, int j)
{
    return i + j;
}

int multiply(float x, float y)
{
    return x * y;
}

int getRows(const std::string &path)
{
    FileParser fp(path);
    return fp.nRow();
}

int getCols(const std::string &path)
{
    FileParser fp(path);
    return fp.nCol();
}

float runCogaps(const std::string &path)
{
    GapsParameters params(path);
    GapsRandomState randState(params.seed);
    GapsResult result(gaps::run(path, params, std::string(), &randState));
    return result.meanChiSq;
}

PYBIND11_MODULE(cogaps, m)
{
    m.doc() = "CoGAPS Python Package";
    m.def("add", &add, "A function which adds two numbers");
    m.def("multiply", &multiply, "A function which multiplies two numbers");
    m.def("get_rows", &getRows, "Get number of rows in file");
    m.def("get_cols", &getCols, "Get number of columns in file");
    m.def("runCogaps", &runCogaps, "Run CoGAPS Algorithm");
}