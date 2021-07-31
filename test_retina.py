from PyCoGAPS import *

path = './src/CoGAPS/inst/extdata/retina_subset_1.h5'

# adata = anndata.read_hdf(path, "counts")

# print(adata.X.shape)

result = CoGAPS(path)

# matrix = pycogaps.Matrix(adata.X)

# prm = pycogaps.GapsParameters(matrix)


# result = pycogaps.runCogapsFromMatrix(matrix, prm)