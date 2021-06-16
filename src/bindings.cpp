//#include "CoGAPS/src/GapsRunner.h"
//#include "CoGAPS/src/utils/GlobalConfig.h"
//#include "CoGAPS/src/GapsParameters.h"
//#include "CoGAPS/src/GapsResult.h"
//#include "CoGAPS/src/math/Random.h"

#include <pybind11/pybind11.h>
#include <iostream>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
namespace py = pybind11;

float runCogaps(const std::string &path)
{
//    GapsParameters params(path);
//    GapsRandomState randState(params.seed);
//    GapsResult result(gaps::run(path, params, std::string(), &randState));
//    return result.meanChiSq;
      std::cout<<"in the runcogaps fxn";
      return 0.0;
}

PYBIND11_MODULE(pycogaps, m)
{
    m.doc() = "CoGAPS Python Package";
    m.def("runCogaps", &runCogaps, "Run CoGAPS Algorithm");
}