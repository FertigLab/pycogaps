#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <numpy/numpyconfig.h>
#include <numpy/ndarrayobject.h>
#include "./Rpackage/src/GapsDispatcher.h"
#include "./Rpackage/src/data_structures/Matrix.h"

namespace p = boost::python;
namespace np = boost::python::numpy;

RowMatrix convertPyMatrix(const np::ndarray pymat)
{
    RowMatrix mat(pymat.shape(0), pymat.shape(1));
    for (int i = 0; i < pymat.shape(0); ++i)
    {
        for (int j = 0; j < pymat.shape(1); ++j)
        {
            mat(i,j) = *pymat[i][j];
        }
    }
    return mat;
}

np::ndarray convertRowMatrix(const RowMatrix &rmat)
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

np::ndarray convertColMatrix(const ColMatrix &cmat)
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

p::tuple run(np::ndarray D, unsigned numPatterns, unsigned seed, unsigned nCores)
{
    Py_Initialize();
    np::initialize();

    RowMatrix convDMatrix = convertPyMatrix(D);

    GapsDispatcher dispatcher(seed);
    dispatcher.setNumPatterns(numPatterns);
    dispatcher.setNumCoresPerSet(nCores);
    dispatcher.loadData(convDMatrix);
    GapsResult result(dispatcher.run());

    np::ndarray AMatrix = convertColMatrix(result.Amean);
    np::ndarray PMatrix = convertRowMatrix(result.Pmean);

    return p::make_tuple(AMatrix, PMatrix);
}

BOOST_PYTHON_MODULE(Cogaps)
{
  p::def("run", run, "run Cogaps");
}
