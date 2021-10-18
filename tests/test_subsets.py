import sys
sys.path.append(".") # Adds higher directory to python modules path.

from PyCoGAPS import *
import subset_data

path = "./data/GIST.csv" 
params = CoParams(path)

adata = toAnndata(path)

print('\n### Test Explicit Subsets ###')
setParam(params, "explicitSets", ['IM00', 'IM02'])
params.setDistributedParams(nSets=2)
sets = subset_data.createSets(adata, params)
print(sets)

print('\n### Test Sample Uniformly Subsets ###')
setParam(params, "explicitSets", None)
sets = subset_data.createSets(adata, params)
print(sets)

print('\n### Test Annotation Weights Subsets ###')
names = ['IM00', 'IM02']
wt = [2, 0.5]
weight = dict(zip(names, wt))
annotation = ['IM00', 'IM02', 'IM00']
params.setAnnotationWeights(annotation, weight)
sets = subset_data.createSets(adata, params)
print(sets)