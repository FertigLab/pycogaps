import anndata
import pycogaps
import scipy.io
import numpy as np

# input is anndata object

path = './data/GIST.mtx'
prm = pycogaps.GapsParameters(path)

adata = anndata.read_mtx(path)
adata_X = adata.X

# mat = scipy.io.mmread('./data/GIST.mtx')

print("\n######### Matrix from np array example ############")
arr = np.array([[10,2],[3,4]], dtype=float)
# arr = np.ones((2,2), dtype=float)*10
print(arr)
sample_mat = pycogaps.Matrix(arr)
np_ret = np.array(sample_mat, copy = False)
print(np_ret)


# print(sample_mat.nRow())
# pycogaps.runCogapsFromMatrix(sample_mat, prm)

# print(mat.shape)
# pycogaps.runCogapsFromMatrix(mat, prm)

print("\n######### Matrix from path example ############")
matrix = pycogaps.Matrix(path, False, False, [])
print(matrix.nRow())
matrix_np = np.array(matrix, copy=False)
print(matrix_np.shape)

print("\n######### Random Matrix example ############")
mat = pycogaps.Matrix(3,3)
print(np.array(mat, copy=False))