import anndata
import pycogaps
import scipy.io
import scipy.sparse
import numpy as np

# placeholder for when we have anndata samples
path = './data/GIST.mtx'
prm = pycogaps.GapsParameters(path)

adata = anndata.read_mtx(path)
adata = adata.X

if scipy.sparse.issparse(adata):
    adata_X = adata.toarray()

sample_mat = pycogaps.Matrix(adata)
np_ret = np.array(sample_mat, copy = False)

print("~~~~~`pycogaps run from path: ~~~~~~~~~~\n")
pycogaps.runCogaps(path, prm)

print("~~~~~~pycogaps run from Matrix: ~~~~~~~~~~\n")
pycogaps.runCogapsFromMatrix(sample_mat, prm)