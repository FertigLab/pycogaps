#include "CoGAPS/src/GapsRunner.h"
#include "CoGAPS/src/utils/GlobalConfig.h"
#include "CoGAPS/src/GapsParameters.h"
#include "CoGAPS/src/GapsResult.h"
#include "CoGAPS/src/math/Random.h"
#include "CoGAPS/src/cpp_tests/catch.h"
#include <pybind11/pybind11.h>
#include <iostream>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
namespace py = pybind11;


float runCogaps(const std::string &path, GapsParameters params)
{
    std::cout<<"1";
    std::cout<<"2";
    GapsRandomState randState(params.seed);
    GapsResult result(gaps::run(path, params, std::string(), &randState));
    return result.meanChiSq;
    std::cout<<"in the runcogaps fxn";
    return 0.0;
}

int runCPPTests()
{
    std::cout<<"in the cpp testing function";
    return 0;
}

PYBIND11_MODULE(pycogaps, m)
{
    m.doc() = "CoGAPS Python Package";
    m.def("runCogaps", &runCogaps, "Run CoGAPS Algorithm");
    m.def("runCPPTests", &runCPPTests, "Run CoGAPS C++ Tests");
    py::class_<GapsParameters>(m, "GapsParameters")
        .def(py::init<const std::string &>())
        .def("print", &GapsParameters::print)
        .def_readwrite("checkpointOutFile", &GapsParameters::checkpointOutFile)
        .def_readwrite("checkpointFile", &GapsParameters::checkpointFile)
        .def_readwrite("seed", &GapsParameters::seed)
        .def_readwrite("nGenes", &GapsParameters::nGenes)
        .def_readwrite("nSamples", &GapsParameters::nSamples)
        .def_readwrite("nPatterns", &GapsParameters::nPatterns)
        .def_readwrite("nIterations", &GapsParameters::nIterations)
        .def_readwrite("maxThreads", &GapsParameters::maxThreads)
        .def_readwrite("outputFrequency", &GapsParameters::outputFrequency)
        .def_readwrite("checkpointInterval", &GapsParameters::checkpointInterval)
        .def_readwrite("snapshotFrequency", &GapsParameters::snapshotFrequency)
        .def_readwrite("alphaA", &GapsParameters::alphaA)
        .def_readwrite("alphaP", &GapsParameters::alphaP)
        .def_readwrite("maxGibbsMassA", &GapsParameters::maxGibbsMassA)
        .def_readwrite("pumpThreshold", &GapsParameters::pumpThreshold)
        .def_readwrite("snapshotPhase", &GapsParameters::snapshotPhase)
        .def_readwrite("useFixedPatterns", &GapsParameters::useFixedPatterns)
        .def_readwrite("subsetData", &GapsParameters::subsetData)
        .def_readwrite("useCheckPoint", &GapsParameters::useCheckPoint)
        .def_readwrite("transposeData", &GapsParameters::transposeData)
        .def_readwrite("printMessages", &GapsParameters::printMessages)
        .def_readwrite("subsetGenes", &GapsParameters::subsetGenes)
        .def_readwrite("printThreadUsage", &GapsParameters::printThreadUsage)
        .def_readwrite("useSparseOptimization", &GapsParameters::useSparseOptimization)
        .def_readwrite("takePumpSamples", &GapsParameters::takePumpSamples)
        .def_readwrite("asynchronousUpdates", &GapsParameters::asynchronousUpdates)
        .def_readwrite("whichMatrixFixed", &GapsParameters::whichMatrixFixed)
        .def_readwrite("workerID", &GapsParameters::workerID)
        .def_readwrite("runningDistributed", &GapsParameters::runningDistributed);
}
