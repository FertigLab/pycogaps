import anndata
import pycogaps
import scipy.io
import scipy.sparse
import numpy as np
import PyCoGAPS
import pandas as pd

# placeholder until we have anndata samples
path = './data/GIST.csv'
prm = pycogaps.GapsParameters(path)

adata = anndata.read_csv(path)
adataX = adata.X

print(adata.obs_names)
print("var: ", adata.var_names)

if scipy.sparse.issparse(adataX):
    adataX = adataX.toarray()

matrix = pycogaps.Matrix(adataX)

print("\n~~~~~~~~~~~~ pycogaps run from path: ~~~~~~~~~~~~~~~\n")
path_result = pycogaps.runCogaps(path, prm)

# converting Matrix object (GapsResult Amean) to numpy array
path_Amean = np.array(path_result.Amean, copy = False)

print("\n~~~~~~~~~~~~ pycogaps run from Matrix: ~~~~~~~~~~~~~~~~~\n")
mat_result = pycogaps.runCogapsFromMatrix(matrix, prm)

mat_Amean = np.array(mat_result.Amean, copy = False)

mat_Pmean = np.array(mat_result.Pmean, copy = False)

patterns = ["Pattern" + str(i) for i in range(1, prm.nPatterns+1)]

# TODO: figure out what column labels should be 
A_mat = pd.DataFrame(data=mat_Amean, index=adata.obs_names, columns=patterns)
adata.obsm = A_mat

P_mat = pd.DataFrame(data=mat_Pmean, index=adata.var_names, columns=patterns)
adata.varm = P_mat
