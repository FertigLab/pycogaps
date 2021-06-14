#include "Rpackage/src/GapsRunner.h"
#include "Rpackage/src/utils/GlobalConfig.h"
#include "Rpackage/src/GapsParameters.h"
#include "Rpackage/src/GapsResult.h"
#include "Rpackage/src/math/Random.h"

#include <pybind11/pybind11.h>
#include <iostream>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
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