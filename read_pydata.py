import anndata
import pycogaps
import scipy.io
import scipy.sparse
import numpy as np
import PyCoGAPS

# placeholder until we have anndata samples
path = './data/GIST.csv'
prm = pycogaps.GapsParameters(path)

adata = anndata.read_csv(path)
adata = adata.X

if scipy.sparse.issparse(adata):
    adata = adata.toarray()
matrix = pycogaps.Matrix(np.array(adata))

print("\n~~~~~~~~~~~~ pycogaps run from path: ~~~~~~~~~~~~~~~\n")
path_result = pycogaps.runCogaps(path, prm)

# converting Matrix object (GapsResult Amean) to numpy array
path_Amean = np.array(path_result.Amean, copy = False)


print("\n~~~~~~~~~~~~ pycogaps run from Matrix: ~~~~~~~~~~~~~~~~~\n")
mat_result = pycogaps.runCogapsFromMatrix(matrix, prm)

mat_Amean = np.array(mat_result.Amean, copy = False)


