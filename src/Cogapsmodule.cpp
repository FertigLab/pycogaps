#include <Python.h>
#include <numpy/numpyconfig.h>
#include <numpy/ndarrayobject.h>
#include <iostream>

static PyObject *test_do_nothing(PyObject *self, PyObject *args)
{
    return Py_None;
}

static PyObject *test_matrices(PyObject *self, PyObject *args)
{
    PyArrayObject *DMatrix;
    PyArrayObject *SMatrix;
    if (!PyArg_ParseTuple(args, "OO", &DMatrix, &SMatrix))
        return NULL;
    return PyTuple_Pack(2, DMatrix, SMatrix);
}

static PyMethodDef Cogaps_methods[] =
  {
    {"test_matrices", test_matrices, METH_VARARGS, NULL},
    {"test_do_nothing", test_do_nothing, METH_NOARGS, NULL},
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
