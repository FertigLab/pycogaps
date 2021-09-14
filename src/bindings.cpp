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
#include "CoGAPS/src/data_structures/Vector.h"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>


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

GapsResult runCogapsFromMatrix(Matrix mat, GapsParameters params, Matrix unc)
{
    GapsRandomState randState(params.seed);
    GapsResult result(gaps::run(mat, params, unc, &randState));
    return result;
}

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

float getElement(Vector v, unsigned i) {
    return v[i];
}


PYBIND11_MODULE(pycogaps, m)
{
    m.doc() = "CoGAPS Python Package";
    m.def("runCogaps", &runCogaps, "Run CoGAPS Algorithm");
    m.def("runCogapsFromMatrix", &runCogapsFromMatrix, "Run CoGAPS Algorithm");
    m.def("runCPPTests", &runCPPTests, "Run CoGAPS C++ Tests");
    m.def("getElement", &getElement, "Get an element of a Vector");
    py::enum_<GapsAlgorithmPhase>(m, "GapsAlgorithmPhase")
        .value("GAPS_EQUILIBRATION_PHASE", GAPS_EQUILIBRATION_PHASE)
        .value("GAPS_SAMPLING_PHASE", GAPS_SAMPLING_PHASE)
        .value("GAPS_ALL_PHASES", GAPS_ALL_PHASES)
        .export_values();
    py::enum_<PumpThreshold>(m, "PumpThreshold")
        .value("PUMP_UNIQUE", PUMP_UNIQUE)
        .value("PUMP_CUT", PUMP_CUT)
        .export_values();
    py::class_<GapsParameters>(m, "GapsParameters")
        .def(py::init<>())
        .def(py::init<const std::string &>())
        .def(py::init<const Matrix&>())
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
        .def_readwrite("runningDistributed", &GapsParameters::runningDistributed)
        .def_readwrite("dataIndicesSubset", &GapsParameters::dataIndicesSubset)
        .def_readwrite("fixedPatterns", &GapsParameters::fixedPatterns)
        .def(py::pickle(
            [](const GapsParameters &prm) {
                return py::make_tuple(prm.checkpointOutFile, prm.checkpointFile, prm.seed, prm.nGenes,
                prm.nSamples, prm.nPatterns, prm.nIterations, prm.maxThreads, prm.outputFrequency,
                prm.checkpointInterval, prm.snapshotFrequency, prm.alphaA, prm.alphaP, prm.maxGibbsMassA,
                prm.maxGibbsMassP, prm.pumpThreshold, prm.snapshotPhase, prm.useFixedPatterns,
                prm.subsetData, prm.useCheckPoint, prm.transposeData, prm.printMessages, prm.subsetGenes,
                prm.printThreadUsage, prm.useSparseOptimization, prm.takePumpSamples, prm.asynchronousUpdates,
                prm.whichMatrixFixed,
                prm.workerID,
                prm.runningDistributed,
                prm.dataIndicesSubset, prm.fixedPatterns);
            },
            [](py::tuple t) {
                if (t.size() != 32){
                    throw std::runtime_error("Invalid state!");
                }
                GapsParameters prm;
                prm.checkpointOutFile    = t[0].cast<std::string>();
                prm.checkpointFile    = t[1].cast<std::string>();
                prm.seed    = t[2].cast<uint32_t>();
                prm.nGenes    = t[3].cast<unsigned>();
                prm.nSamples    = t[4].cast<unsigned>();
                prm.nPatterns    = t[5].cast<unsigned>();
                prm.nIterations    = t[6].cast<unsigned>();
                prm.maxThreads    = t[7].cast<unsigned>();
                prm.outputFrequency    = t[8].cast<unsigned>();
                prm.checkpointInterval    = t[9].cast<unsigned>();
                prm.snapshotFrequency    = t[10].cast<unsigned>();
                prm.alphaA    = t[11].cast<float>();
                prm.alphaP    = t[12].cast<float>();
                prm.maxGibbsMassA    = t[13].cast<float>();
                prm.maxGibbsMassP    = t[14].cast<float>();
                prm.pumpThreshold    = t[15].cast<PumpThreshold>();
                prm.snapshotPhase    = t[16].cast<GapsAlgorithmPhase>();
                prm.useFixedPatterns    = t[17].cast<bool>();
                prm.subsetData    = t[18].cast<bool>();
                prm.useCheckPoint    = t[19].cast<bool>();
                prm.transposeData    = t[20].cast<bool>();
                prm.printMessages    = t[21].cast<bool>();
                prm.subsetGenes    = t[22].cast<bool>();
                prm.printThreadUsage    = t[23].cast<bool>();
                prm.useSparseOptimization    = t[24].cast<bool>();
                prm.takePumpSamples    = t[25].cast<bool>();
                prm.asynchronousUpdates    = t[26].cast<bool>();
                prm.whichMatrixFixed    = t[27].cast<char>();
                prm.workerID    = t[28].cast<unsigned>();
                prm.runningDistributed    = t[29].cast<bool>();
                prm.dataIndicesSubset    = t[30].cast<std::vector<unsigned>>();
                prm.fixedPatterns    = t[31].cast<Matrix>();
                return prm;
            }
        ));
    m.def("getBuildReport", &getBuildReport, "Return build report.");
    m.def("isCheckpointsEnabled", &isCheckpointsEnabled, "Return whether checkpoints enabled.");
    m.def("isCompiledWithOpenMPSupport", &isCompiledWithOpenMPSupport, "Return whether compiled with Open MP Support.");
    m.def("getFileInfo", &getFileInfo, "Get info of inputted file.");

    py::class_<GapsResult>(m, "GapsResult")
        .def(py::init<>())
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
        .def_readwrite("averageQueueLengthP", &GapsResult::averageQueueLengthP)
        .def(py::pickle(
        [](const GapsResult &r) { // __getstate__
                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(r.Amean, r.Asd, r.Pmean, r.Psd, r.pumpMatrix,
                r.meanPatternAssignment, r.equilibrationSnapshotsA, r.equilibrationSnapshotsP,
                r.samplingSnapshotsA, r.samplingSnapshotsP, r.chisqHistory, r.atomHistoryA,
                r.atomHistoryP, r.totalUpdates, r.seed, r.totalRunningTime, r.meanChiSq,
                r.averageQueueLengthA, r.averageQueueLengthP);
            },
            [](py::tuple t) { // __setstate__
//                if (t.size() != 2)
//                    throw std::runtime_error("Invalid state!");

                /* Create a new C++ instance */
                GapsResult r;
                r.Amean = t[0].cast<Matrix>();
                r.Asd = t[1].cast<Matrix>();
                r.Pmean = t[2].cast<Matrix>();
                r.Psd = t[3].cast<Matrix>();
                r.pumpMatrix = t[4].cast<Matrix>();
                r.meanPatternAssignment = t[5].cast<Matrix>();
                r.equilibrationSnapshotsA = t[6].cast<std::vector<Matrix>>();
                r.equilibrationSnapshotsP = t[7].cast<std::vector<Matrix>>();
                r.samplingSnapshotsA = t[8].cast<std::vector<Matrix>>();
                r.samplingSnapshotsP = t[9].cast<std::vector<Matrix>>();
                r.chisqHistory = t[10].cast<std::vector<float>>();
                r.atomHistoryA = t[11].cast<std::vector<unsigned>>();
                r.atomHistoryP = t[12].cast<std::vector<unsigned>>();
                r.totalUpdates = t[13].cast<uint64_t>();
                r.seed = t[14].cast<uint32_t>();
                r.totalRunningTime = t[15].cast<unsigned>();
                r.meanChiSq = t[16].cast<float>();
                r.averageQueueLengthA = t[17].cast<float>();
                r.averageQueueLengthP = t[18].cast<float>();
                return r;
            }
        ));
  
    py::class_<Vector>(m, "Vector")
        .def(py::init<unsigned &>())
        .def("size", &Vector::size);


    py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
        .def(py::init<>())
        .def(py::init<unsigned &, unsigned &>())
        .def(py::init<const Matrix &, bool &, bool &,
        std::vector<unsigned> &>())
        .def(py::init<const std::string &, bool &, bool &,
        std::vector<unsigned> &>())
        // Matrix constructed from numpy array
        .def(py::init([](py::array_t<float> b) {
            py::buffer_info info = b.request();

            if (info.ndim != 2)
            {
                throw std::runtime_error("Incompatible buffer dimension");
            }

            Matrix mat = Matrix(info.shape[0], info.shape[1]);
            float *ptr = static_cast<float *>(info.ptr);

            for(int i = 0; i < info.shape[0]; i++)
            {
                for (int j = 0; j < info.shape[1]; j++)
                {
                    mat.operator()(i,j) = ptr[i*info.shape[1] + j];
                }
            }

            return mat.getMatrix();
         }))

        .def("nCol", &Matrix::nCol)
        .def("nRow", &Matrix::nRow)
        .def("getCol", static_cast<Vector& (Matrix::*)(unsigned)>(&Matrix::getCol), "Get a column of the matrix")

        .def_buffer([](Matrix &m) -> py::buffer_info {
            return py::buffer_info(
                &(m.getMatrix().operator()(0,0)),
                sizeof(float),
                py::format_descriptor<float>::format(),
                2,
                {m.nRow(), m.nCol()},
                {sizeof(float) * m.nCol(), sizeof(float)}
            );
        })
        .def(py::pickle(
            [](const Matrix &m) { // __getstate__
                std::vector<std::vector<float>> a(m.nRow(), std::vector<float>(m.nCol()));
                for (int i = 0; i < m.nRow(); i++) {
                    for (int j = 0; j < m.nCol(); j++) {
                        a[i][j] = m.operator()(i,j);
                    }
                }
                return py::make_tuple(m.nCol(), m.nRow(), a);
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state!");

                /* Create a new C++ instance */
                unsigned ncol = t[0].cast<unsigned>();
                unsigned nrow = t[1].cast<unsigned>();
                Matrix m(nrow, ncol);
                std::vector<std::vector<float>> ptr = t[2].cast<std::vector<std::vector<float>>>();
                for(int i = 0; i < (int)nrow; i++)
                {
                    for (int j = 0; j < (int)ncol; j++)
                    {
                        m.operator()(i,j) = ptr[i][j];
                    }
                }
                return m;
            }
        ));
}
