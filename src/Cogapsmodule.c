#include <Python.h>
//#include <numpy/numpyconfig.h>
//#include <numpy/ndarrayobject.h>
//#include "Rpackage/src/GapsDispatcher.h"
//#include "Rpackage/src/data_structures/Matrix.h"
//#include "Cogapssupport.hpp"
/*
RowMatrix convertPyMatrix(const PyArrayObject *pymat)
{
    RowMatrix mat(pymat->dimensions[0], pymat->dimensions[1]);
    for (unsigned i = 0; i < pymat->dimensions[0]; ++i)
    {
        for (unsigned j = 0; j < pymat->dimensions[1]; ++j)
        {
            mat(i,j) = *(float *)(pymat->data + i*pymat->strides[0] + j*pymat->strides[1]);
        }
    }
    return mat;
}

PyArrayObject *convertRowMatrix(const RowMatrix &rmat)
{
    int dims[2];
    dims[0] = rmat.nRow();
    dims[1] = rmat.nCol();
    PyArrayObject *pymat = (PyArrayObject*) PyArray_FromDims(2, dims, PyArray_FLOAT);
    for (unsigned i = 0; i < rmat.nRow(); ++i)
    {
        for (unsigned j = 0; j < rmat.nCol(); ++j)
        {
            *(pymat->data + i*pymat->strides[0] + j*pymat->strides[1]) = rmat.operator()(i,j);
        }
    }

    return pymat;
}

PyArrayObject *convertColMatrix(const ColMatrix &cmat)
{
    int dims[2];
    dims[0] = cmat.nRow();
    dims[1] = cmat.nCol();
    PyArrayObject *pymat = (PyArrayObject*) PyArray_FromDims(2, dims, PyArray_FLOAT);
    for (unsigned i = 0; i < cmat.nRow(); ++i)
    {
        for (unsigned j = 0; j < cmat.nCol(); ++j)
        {
            *(pymat->data + i*pymat->strides[0] + j*pymat->strides[1]) = cmat.operator()(i,j);
        }
    }

    return pymat;
}

static PyObject *test_matrices(PyObject *self, PyObject *args)
{
    PyArrayObject *DMatrix;
    PyArrayObject *SMatrix;
    if (!PyArg_ParseTuple(args, "OO", &DMatrix, &SMatrix))
        return NULL;
    return PyTuple_Pack(2, DMatrix, SMatrix);
}

static PyObject *run(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyArrayObject *DMatrix;
    unsigned numPatterns;
    unsigned seed;
    unsigned nCores;

    static char *keywords[] = {"D", "numPatterns", "seed", "nCores", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|bbb", keywords, &DMatrix, &numPatterns, &seed, &nCores))
        return NULL;

    RowMatrix convDMatrix = convertPyMatrix(DMatrix);

    GapsDispatcher dispatcher(seed);
    dispatcher.setNumPatterns(numPatterns);
    dispatcher.setNumCoresPerSet(nCores);
    dispatcher.loadData(convDMatrix);
    GapsResult result(dispatcher.run());

    PyArrayObject *AMatrix = convertColMatrix(result.Amean);
    PyArrayObject *PMatrix = convertRowMatrix(result.Pmean);

    return PyTuple_Pack(2, AMatrix, PMatrix);
}
*/

static PyObject *test_do_nothing(PyObject *self, PyObject *args)
{
    return Py_None;
}

static PyMethodDef Cogaps_methods[] =
{
    //{"test_matrices", test_matrices, METH_VARARGS, NULL},
    {"test_do_nothing", test_do_nothing, METH_NOARGS, NULL},
    //{"run", (PyCFunction)run, METH_VARARGS|METH_KEYWORDS, NULL},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef Cogapsmodule =
{
    PyModuleDef_HEAD_INIT, "Cogaps", NULL, -1, Cogaps_methods
};

PyMODINIT_FUNC PyInit_Cogaps(void)
{
    import_array();
    return PyModule_Create(&Cogapsmodule);
}
