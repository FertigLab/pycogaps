#include "CoGAPS/src/GapsRunner.h"
#include "CoGAPS/src/utils/GlobalConfig.h"
#include "CoGAPS/src/GapsParameters.h"
#include "CoGAPS/src/GapsResult.h"
#include "CoGAPS/src/math/Random.h"
#include "CoGAPS/src/cpp_tests/catch.h"
#include "CoGAPS/src/file_parser/FileParser.h"
#include "CoGAPS/src/include/boost/algorithm/string/join.hpp"
#include "CoGAPS/src/GapsStatistics.h"
#include "CoGAPS/src/data_structures/Matrix.h"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>


#include <iostream>
#include <vector>
#include <iterator>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
namespace py = pybind11;


// overload, with given params
GapsResult runCogaps(const std::string &path, GapsParameters params)
{
    GapsRandomState randState(params.seed);
    GapsResult result(gaps::run(path, params, std::string(), &randState));
    return result;
}

//GapsResult runCogapsFromMatrix(Matrix mat, GapsParameters params)
//{
//    // TODO: implement
//}

std::string getBuildReport()
{
    return buildReport();
}

bool isCheckpointsEnabled()
{
#ifdef GAPS_DISABLE_CHECKPOINTS
    return false;
#else
    return true;
#endif
}

bool isCompiledWithOpenMPSupport()
{
#ifdef __GAPS_OPENMP__
    return true;
#else
    return false;
#endif
}

std::string getFileInfo(const std::string &path)
{
    FileParser fp(path);
    return "dimensions: " + std::to_string(fp.nRow()) + ", " + std::to_string(fp.nCol())
    + "\nrowNames: " + boost::algorithm::join(fp.rowNames(), " ") + "\ncolNames: " + boost::algorithm::join(fp.colNames(), " ");
    return 0;
}

void runCPPTests()
{
    std::cout << "running CPPTests";
}

PYBIND11_MODULE(pycogaps, m)
{
    m.doc() = "CoGAPS Python Package";
    m.def("runCogaps", &runCogaps, "Run CoGAPS Algorithm");
    m.def("runCPPTests", &runCPPTests, "Run CoGAPS C++ Tests");
    py::enum_<GapsAlgorithmPhase>(m, "GapsAlgorithmPhase")
        .value("GAPS_EQUILIBRATION_PHASE", GAPS_EQUILIBRATION_PHASE)
        .value("GAPS_SAMPLING_PHASE", GAPS_SAMPLING_PHASE)
        .value("GAPS_ALL_PHASES", GAPS_ALL_PHASES)
        .export_values();
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
        .def_readwrite("maxGibbsMassP", &GapsParameters::maxGibbsMassP)
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
    m.def("getBuildReport", &getBuildReport, "Return build report.");
    m.def("isCheckpointsEnabled", &isCheckpointsEnabled, "Return whether checkpoints enabled.");
    m.def("isCompiledWithOpenMPSupport", &isCompiledWithOpenMPSupport, "Return whether compiled with Open MP Support.");
    m.def("getFileInfo", &getFileInfo, "Get info of inputted file.");

    py::class_<GapsResult>(m, "GapsResult")
        .def(py::init<const GapsStatistics &>())
        .def("writeToFile", &GapsResult::writeToFile)
        .def_readwrite("Amean", &GapsResult::Amean)
        .def_readwrite("Asd", &GapsResult::Asd)
        .def_readwrite("Pmean", &GapsResult::Pmean)
        .def_readwrite("Psd", &GapsResult::Psd)
        .def_readwrite("pumpMatrix", &GapsResult::pumpMatrix)
        .def_readwrite("meanPatternAssignment", &GapsResult::meanPatternAssignment) // Matrix
        .def_readwrite("equilibrationSnapshotsA", &GapsResult::equilibrationSnapshotsA)
        .def_readwrite("equilibrationSnapshotsP", &GapsResult::equilibrationSnapshotsP)
        .def_readwrite("samplingSnapshotsA", &GapsResult::samplingSnapshotsA)
        .def_readwrite("samplingSnapshotsP", &GapsResult::samplingSnapshotsP) // std::vector<Matrix>
        .def_readwrite("chisqHistory", &GapsResult::chisqHistory)
        .def_readwrite("atomHistoryA", &GapsResult::atomHistoryA)
        .def_readwrite("atomHistoryP", &GapsResult::atomHistoryP)
        .def_readwrite("totalUpdates", &GapsResult::totalUpdates)
        .def_readwrite("seed", &GapsResult::seed)
        .def_readwrite("totalRunningTime", &GapsResult::totalRunningTime)
        .def_readwrite("meanChiSq", &GapsResult::meanChiSq)
        .def_readwrite("averageQueueLengthA", &GapsResult::averageQueueLengthA)
        .def_readwrite("averageQueueLengthP", &GapsResult::averageQueueLengthP);

    py::class_<Matrix>(m, "Matrix")
        .def(py::init<>())
        .def(py::init<unsigned &, unsigned &>())
        .def(py::init<const Matrix &, bool &, bool &,
        std::vector<unsigned> &>())
        .def(py::init<const std::string &, bool &, bool &,
        std::vector<unsigned> &>());
}
