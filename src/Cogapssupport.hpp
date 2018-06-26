#include <numpy/numpyconfig.h>
#include <numpy/ndarrayobject.h>
#include "Rpackage/src/data_structures/Matrix.h"

RowMatrix convertPyMatrix(const PyArrayObject *pymat);
PyArrayObject *convertRowMatrix(const RowMatrix &rmat);
PyArrayObject *convertColMatrix(const ColMatrix &cmat);
