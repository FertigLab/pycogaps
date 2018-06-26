#include "Cogapssupport.hpp"

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
            //*PyArray_GETPTR2(mat, i, j) = cmat.operator()(i,j);
        }
    }

    return pymat;
}
