import sys
sys.path.append(".") # Adds higher directory to python modules path.

import anndata
import pycogaps
import scipy.io
import scipy.sparse
import numpy as np
from PyCoGAPS import *

# placeholder until we have anndata samples
# maybe also read files into an anndata object?
path = './data/GIST.csv'
prm = pycogaps.GapsParameters(path)

adata = anndata.read_csv(path)
adataX = adata.X

if scipy.sparse.issparse(adataX):
    adataX = adataX.toarray()

# create Matrix object from anndata X
matrix = pycogaps.Matrix(adataX)

result = pycogaps.runCogapsFromMatrix(matrix, prm)

# convert Amean and Pmean results to numpy arrays
Amean = toNumpy(result.Amean)
Pmean = toNumpy(result.Pmean)

# anndata labels
print('obs names: ', adata.obs_names)
print('var names: ', adata.var_names)
pattern_labels = ["Pattern" + str(i) for i in range(1, prm.nPatterns+1)]

# load adata obs and var from Amean and Pmean results
A_mat = pd.DataFrame(data=Amean, index=adata.obs_names, columns=pattern_labels)
adata.obs = A_mat

P_mat = pd.DataFrame(data=Pmean, index=adata.var_names, columns=pattern_labels)
adata.var = P_mat

print("~~~~~~~~~~~~~ testing CoGAPS Pattern Markers ~~~~~~~~~~~~~~")
print(patternMarkers(adata))