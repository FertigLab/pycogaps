#include "Rpackage/src/GapsRunner.h"
#include "Rpackage/src/utils/GlobalConfig.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

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
    m.def("runCogaps", &runCogaps, "Run CoGAPS Algorithm");
}