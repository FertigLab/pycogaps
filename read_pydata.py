import anndata
import pycogaps
import scipy.io
import numpy as np

# input is anndata object

path = './data/GIST.csv'
prm = pycogaps.GapsParameters(path)

adata = anndata.read_csv(path)
adata_X = adata.X

# mat = scipy.io.mmread('./data/GIST.mtx')

print("\n######### Matrix from np array example ############")
sample_mat = pycogaps.Matrix(adata_X)
print("finished calling constructor!\n")
np_ret = np.array(sample_mat, copy = False)
print("array to Matrix: \n")
print(np_ret)

# print(sample_mat.nRow())
# pycogaps.runCogapsFromMatrix(sample_mat, prm)

# print(mat.shape)
# pycogaps.runCogapsFromMatrix(mat, prm)

# print("\n######### Matrix from path example ############")
# matrix = pycogaps.Matrix(path, False, False, [])
# matrix_np = np.array(matrix, copy=False)
# print(matrix_np)

# print("\n######### Random Matrix example ############")
# mat = pycogaps.Matrix(3,3)
# print(np.array(mat, copy=False))