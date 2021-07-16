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

int containsZeros(Matrix m) {
    for(int i=0; i < m.nCol(); i++) {
        Vector vec = m.getCol(i);
        for(int j=0; j<vec.size(); j++){
            if (vec[i] == 0){
                return 1;
            }
        }
    }
    return 0;
}

Matrix replaceZeros(Matrix m) {
    for(int i=0; i < m.nCol(); i++) {
        Vector vec = m.getCol(i);
        for(int j=0; j<vec.size(); j++){
            if (vec[i] == 0){
                vec[i] = 1e-6;
            }
        }
    }
    return m;
}

Matrix divideMatrices(Matrix m1, Matrix m2){
    int nrow = m1.nRow();
    int ncol = m1.nCol();
    Matrix *mat = new Matrix(nrow, ncol);
    Matrix retmat = *mat;
    if (m1.nCol() != m2.nCol() || m1.nRow() != m2.nRow()){
        std::cout<<"Dimensions are not equal! Aborting.";
        return retmat;
    }
    for (int i = 0; i < nrow; i++) {
        for(int j = 0; j < ncol; j++) {
            retmat.getCol(j)[i] = m1.getCol(j)[i] / m2.getCol(j)[i];
        }
    }
    return retmat;
}

Matrix multiplyMatrices(Matrix m1, Matrix m2) {
    int m1rows = m1.nRow();
    int m1cols = m1.nCol();
    int m2cols = m2.nCol();
    int m2rows = m2.nRow();
    Matrix *mat = new Matrix(m1rows, m2cols);
    Matrix retmat = *mat;
    if (m1cols != m2rows) {
        std::cout<<"Matrices are not of proper dimensions. Please make sure m1 has the same number of rows as m2 has columns.";
        return retmat;
    }
    for (int i=0; i<m1rows; i++) {
        for (int j=0; j< m1cols; j++) {
            retmat(i,j) += m1(i,j) * m2(j,i);
        }
    }
    return retmat;
}

Matrix transposeMatrix(Matrix mat) {
    int rows = mat.nRow();
    int cols = mat.nCol();
    Matrix * pmat = new Matrix(cols, rows);
    Matrix retmat = *pmat;

    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            retmat(i,j) = mat(j,i);
        }
    }
    return retmat;
}



PYBIND11_MODULE(pycogaps, m)
{
    m.doc() = "CoGAPS Python Package";
    m.def("runCogaps", &runCogaps, "Run CoGAPS Algorithm");
    m.def("runCogapsFromMatrix", &runCogapsFromMatrix, "Run CoGAPS Algorithm");
    m.def("runCPPTests", &runCPPTests, "Run CoGAPS C++ Tests");
    m.def("getElement", &getElement, "Get an element of a Vector");
    m.def("containsZeros", &containsZeros, "Check whether a Matrix contains zeros");
    m.def("replaceZeros", &replaceZeros, "Replace a Matrix's zeros with small values");
    m.def("divideMatrices", &divideMatrices, "Divide m1 / m2 element-wise; return result");
    m.def("multiplyMatrices", &multiplyMatrices, "Multiply m1*m2, return result");
    m.def("transposeMatrix", &transposeMatrix, "Transpose a matrix");
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
  
        .def("nCol", &Matrix::nCol)
        .def("nRow", &Matrix::nRow)
        // .def("getRow", &getRow, "Get a row of the matrix");
        .def("getCol", static_cast<Vector& (Matrix::*)(unsigned)>(&Matrix::getCol), "Get a column of the matrix")

        .def_buffer([](Matrix &m) -> py::buffer_info {
            return py::buffer_info(
                &(m.getMatrix().operator()(0,0)),
                sizeof(float),
                py::format_descriptor<float>::format(),
                2,
                {m.nRow(), m.nCol()},
                {sizeof(float), sizeof(float)*m.nCol()}
            );
        })

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
        }));
}
