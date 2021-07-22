'''
check if equivalent outputs attained with same seed value
'''

from PyCoGAPS import *
import anndata

path = './data/GIST.csv'
params = GapsParameters(path)
setParam(params, "seed", 42) 

adata = anndata.read_csv(path)
adataX = adata.X

if scipy.sparse.issparse(adataX):
    adataX = adataX.toarray()

res1 = CoGAPS(path, params=params)
res2 = CoGAPS(path, params=params)

res1_Amean = toNumpy(res1.Amean)
res2_Amean = toNumpy(res2.Amean)

res1_Pmean = toNumpy(res1.Pmean)
res2_Pmean = toNumpy(res2.Pmean)

print((res1_Amean == res2_Amean).all())
print((res1_Pmean == res2_Pmean).all())
