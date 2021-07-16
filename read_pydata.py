import anndata
import pycogaps
import scipy.io
import scipy.sparse
import numpy as np
from PyCoGAPS import *
import pandas as pd
import scanpy as sc

# placeholder until we have anndata samples
path = './data/GIST.csv'
prm = pycogaps.GapsParameters(path)

adata = anndata.read_csv(path)
adataX = adata.X

if scipy.sparse.issparse(adataX):
    adataX = adataX.toarray()

matrix = pycogaps.Matrix(adataX)

all_vecdata = np.empty((matrix.nRow(), matrix.nCol()))
for i in range(matrix.nCol()):
    vector = matrix.getCol(i)
    vecdata = []
    for j in range(vector.size()):
        vecdata.append(getElement(vector, j))
    all_vecdata[:,i] = vecdata

print(all_vecdata)


print("\n~~~~~~~~~~~~ pycogaps run from path: ~~~~~~~~~~~~~~~\n")
path_result = pycogaps.runCogaps(path, prm)

all_vecdata = np.empty((path_result.Amean.nRow(), path_result.Amean.nCol()))
for i in range(path_result.Amean.nCol()):
    vector = path_result.Amean.getCol(i)
    vecdata = []
    for j in range(vector.size()):
        vecdata.append(getElement(vector, j))
    all_vecdata[:,i] = vecdata

print(all_vecdata)


all_vecdata = np.empty((path_result.Pmean.nRow(), path_result.Pmean.nCol()))
for i in range(path_result.Pmean.nCol()):
    vector = path_result.Pmean.getCol(i)
    vecdata = []
    for j in range(vector.size()):
        vecdata.append(getElement(vector, j))
    all_vecdata[:,i] = vecdata

print(all_vecdata)

# # converting Matrix object (GapsResult Amean) to numpy array
# path_Amean = np.array(path_result.Amean, copy = False)


print("\n~~~~~~~~~~~~ pycogaps run from Matrix: ~~~~~~~~~~~~~~~~~\n")
mat_result = pycogaps.runCogapsFromMatrix(matrix, prm)

all_vecdata = np.empty((mat_result.Amean.nRow(), mat_result.Amean.nCol()))
for i in range(mat_result.Amean.nCol()):
    vector = mat_result.Amean.getCol(i)
    vecdata = []
    for j in range(vector.size()):
        vecdata.append(getElement(vector, j))
    all_vecdata[:,i] = vecdata

print(all_vecdata)


all_vecdata = np.empty((mat_result.Pmean.nRow(), mat_result.Pmean.nCol()))
for i in range(mat_result.Pmean.nCol()):
    vector = mat_result.Pmean.getCol(i)
    vecdata = []
    for j in range(vector.size()):
        vecdata.append(getElement(vector, j))
    all_vecdata[:,i] = vecdata

print(all_vecdata)

# patterns = ["Pattern" + str(i) for i in range(1, prm.nPatterns+1)]

# # TODO: figure out what column labels should be 
# A_mat = pd.DataFrame(data=mat_Amean, index=adata.obs_names, columns=patterns)
# adata.obs = A_mat

# P_mat = pd.DataFrame(data=mat_Pmean, index=adata.var_names, columns=patterns)
# adata.var = P_mat

# print("~~~~ testing CoGAPS Pattern Markers ~~~~~~")
# # print(patternMarkers(adata))

# print(mat_Amean)
# print(mat_Pmean)

# mat_result.writeToFile('./data/')
