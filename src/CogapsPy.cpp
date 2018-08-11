#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <numpy/numpyconfig.h>
#include <numpy/ndarrayobject.h>
#include "Rpackage/src/GapsRunner.h"
#include "Rpackage/src/math/Random.h"
#include <unistd.h>
#include <iostream>

namespace p = boost::python;
namespace np = boost::python::numpy;

np::ndarray convertToPyMat(const RowMatrix &rmat)
{
    p::tuple dims = p::make_tuple(rmat.nRow(), rmat.nCol());
    np::ndarray pymat = np::empty(dims, np::dtype::get_builtin<float>());
    
    for (unsigned i = 0; i < rmat.nRow(); ++i)
    {
        for (unsigned j = 0; j < rmat.nCol(); ++j)
        {
            pymat[i][j] = rmat(i,j);
        }
    }

    return pymat;
}

np::ndarray convertToPyMat(const ColMatrix &cmat)
{
    p::tuple dims = p::make_tuple(cmat.nRow(), cmat.nCol());
    np::ndarray pymat = np::empty(dims, np::dtype::get_builtin<float>());
    
    for (unsigned i = 0; i < cmat.nRow(); ++i)
    {
        for (unsigned j = 0; j < cmat.nCol(); ++j)
        {
            pymat[i][j] = cmat(i,j);
        }
    }

    return pymat;
}

p::tuple CoGAPS(std::string dataPath, unsigned numPatterns, unsigned maxIterations)
{
  GapsRunner runner(dataPath, false, numPatterns, false,
        std::vector<unsigned>(1));

  gaps::random::setSeed(123);
  runner.recordSeed(123);
  runner.setMaxIterations(maxIterations);
  runner.setSparsity(0.01f, 0.01f, false);
  runner.setMaxGibbsMass(100.f, 100.f);
  runner.setMaxThreads(1);
  runner.setPrintMessages(true);
  runner.setOutputFrequency(250);
  runner.setCheckpointOutFile("gaps_checkpoint.out");
  runner.setCheckpointInterval(0);
  
  GapsResult result = runner.run();
  
  return p::make_tuple(convertToPyMat(result.Amean), convertToPyMat(result.Asd), convertToPyMat(result.Pmean), convertToPyMat(result.Psd));
}

BOOST_PYTHON_MODULE(CogapsPy)
{
    Py_Initialize();
    np::initialize();
    p::def("CoGAPS", CoGAPS, "");
}
