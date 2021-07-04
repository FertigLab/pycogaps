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

# matrix can be returned as a np array
# arr = np.array([[1,2],[3,4]], dtype=float)
arr = np.ones((5,5), dtype=float)
print(arr.dtype)
# returns Matrix object
sample_mat = pycogaps.Matrix(arr)
print(sample_mat.nRow())
pycogaps.runCogapsFromMatrix(sample_mat, prm)

# print(mat.shape)
# pycogaps.runCogapsFromMatrix(mat, prm)