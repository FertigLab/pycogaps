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

GapsResult runCogapsFromMatrix(Matrix mat, GapsParameters params)
{
    GapsRandomState randState(params.seed);
    GapsResult result(gaps::run(mat, params, Matrix(), &randState));
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
    // m.def("containsZeros", &containsZeros, "Check whether a Matrix contains zeros");
    // m.def("replaceZeros", &replaceZeros, "Replace a Matrix's zeros with small values");
    // m.def("divideMatrices", &divideMatrices, "Divide m1 / m2 element-wise; return result");
    // m.def("multiplyMatrices", &multiplyMatrices, "Multiply m1*m2, return result");
    // m.def("transposeMatrix", &transposeMatrix, "Transpose a matrix");
    py::enum_<GapsAlgorithmPhase>(m, "GapsAlgorithmPhase")
        .value("GAPS_EQUILIBRATION_PHASE", GAPS_EQUILIBRATION_PHASE)
        .value("GAPS_SAMPLING_PHASE", GAPS_SAMPLING_PHASE)
        .value("GAPS_ALL_PHASES", GAPS_ALL_PHASES)
        .export_values();
    py::class_<GapsParameters>(m, "GapsParameters")
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
            [](const GapsParameters &p) { // __getstate__
                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(p.fixedPatterns, p.dataIndicesSubset, p.maxThreads);
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state!");

                /* Create a new C++ instance */
                GapsParameters p("./data/GIST.csv");
                p.maxThreads = t[2].cast<int>();
//                /* Assign any additional state */
//                p.setExtra(t[1].cast<int>());

                return p;
            }
        ));
    m.def("getBuildReport", &getBuildReport, "Return build report.");
    m.def("isCheckpointsEnabled", &isCheckpointsEnabled, "Return whether checkpoints enabled.");
    m.def("isCompiledWithOpenMPSupport", &isCompiledWithOpenMPSupport, "Return whether compiled with Open MP Support.");
    m.def("getFileInfo", &getFileInfo, "Get info of inputted file.");

    py::class_<GapsResult>(m, "GapsResult")
        .def(py::init<const GapsStatistics &>())
        .def(py::init<>())
        .def(py::pickle([](const GapsResult &p) { // __getstate__
            /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(p.Amean, p.Asd, p.Pmean, p.Psd, p.pumpMatrix, p.meanPatternAssignment, p.equilibrationSnapshotsA,
                p.equilibrationSnapshotsP, p.samplingSnapshotsA, p.samplingSnapshotsP, p.chisqHistory, p.atomHistoryA, p.atomHistoryP,
                p.totalUpdates, p.seed, p.totalRunningTime, p.meanChiSq, p.averageQueueLengthA, p.averageQueueLengthP);
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state!");

            /* Create a new C++ instance amd reassign object state completely */
                GapsResult p;
//                p.Amean = py::cast<Matrix>(t[0]);
                p.seed = t[14].cast<std::uint64_t>();
                return p;
            }
         ))
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
  
    py::class_<Vector>(m, "Vector")
        .def(py::init<unsigned &>())
        .def("size", &Vector::size);


    py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
        .def(py::init<const Matrix &, bool &, bool &, std::vector<unsigned> &>())
        .def(py::init<>())
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
        .def(py::init<unsigned &, unsigned &>())
        .def(py::init<const std::string &, bool &, bool &,
        std::vector<unsigned> &>())

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
            [](const Matrix &mat) { // __getstate__
                /* Return a tuple that fully encodes the state of the object */
                std::cout << "pickling matrix...";
                int cols = mat.nCol();
                int rows = mat.nRow();
                std::cout << "cols: " << cols;
                std::cout << "rows: " << rows;
                std::vector<std::vector<int>> a(rows, std::vector<int>(cols));
                std::cout << "C array:\n";
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        a[i][j] = mat.operator()(i,j);
                        std::cout << a[i][j];
                    }
                    std::cout << std::endl;
                }
                return py::make_tuple(a, rows, cols);
            },
            [](py::tuple t) { // __setstate__
//                if (t.size() != 2)
//                    throw std::runtime_error("Invalid state!");
                std::cout << "trying to unpickle";
                int rows = t[1].cast<int>();
                int cols = t[2].cast<int>();
                std::cout << "cols: " << cols;
                std::cout << "rows: " << rows;

                Matrix mat = Matrix(rows, cols);
                std::vector<std::vector<int>> ptr = t[0].cast<std::vector<std::vector<int>>>();

                for(int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        std::cout << ptr[i][j] << " ";
                        mat.operator()(i,j) = ptr[i][j];
                    }
                    std::cout << std::endl;
                }
                return mat;
            }
        ));
}
