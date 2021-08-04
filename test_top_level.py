from PyCoGAPS import *

def no_na_in_result(result):
    return not (np.isnan(toNumpy(result['GapsResult'].Amean)).any()+
    np.isnan(toNumpy(result['GapsResult'].Asd)).any()+
    np.isnan(toNumpy(result['GapsResult'].Pmean)).any()+
    np.isnan(toNumpy(result['GapsResult'].Psd)).any())

csv_path = "./data/GIST.csv" 
mtx_path = "./data/GIST.mtx" 
tsv_path = "./data/GIST.tsv" 

# standard running
print('\n ### Testing Standard Run ###\n')

csv_params = CoParams(csv_path)
mtx_params = CoParams(mtx_path)
tsv_params = CoParams(tsv_path)

adata = anndata.read_csv(csv_path)


setParams(csv_params, {"nIterations": 100,
                    "nPatterns": 7,
                    'hdfKey': 'counts'})
setParams(mtx_params, {"nIterations": 100})
setParams(tsv_params, {"nIterations": 100})

csv_params.print_all()

res = [None] * 3
res[0] = CoGAPS(csv_path, csv_params, outputFrequency=50, messages=False)
res[1] = CoGAPS(mtx_path, mtx_params, outputFrequency=50, messages=False)
res[2] = CoGAPS(tsv_path, tsv_params, outputFrequency=50, messages=False)

for r in res:
    assert(no_na_in_result(r) == True)

assert(toNumpy(res[0]['GapsResult'].Amean).shape[0] == 1363)
assert(toNumpy(res[0]['GapsResult'].Amean).shape[1] == 7)
assert(toNumpy(res[0]['GapsResult'].Pmean).shape[0] == 9)
assert(toNumpy(res[0]['GapsResult'].Pmean).shape[1] == 7)

# transposing data
print('\n ### Testing Transpose Run ###\n')
res = [None] * 3
res[0] = CoGAPS(csv_path, transposeData=True, outputFrequency=50, messages=False)
res[1] = CoGAPS(mtx_path, transposeData=True, outputFrequency=50, messages=False)
res[2] = CoGAPS(tsv_path, transposeData=True, outputFrequency=50, messages=False)

for r in res:
    assert(no_na_in_result(r) == True)

assert(toNumpy(res[0]['GapsResult'].Amean).shape[0] == 9)
assert(toNumpy(res[0]['GapsResult'].Amean).shape[1] == 3)
assert(toNumpy(res[0]['GapsResult'].Pmean).shape[0] == 1363)
assert(toNumpy(res[0]['GapsResult'].Pmean).shape[1] == 3)

# multiple threads
print('\n ### Testing Multiple Threads Run ###\n')
res = [None] * 3
res[0] = CoGAPS(csv_path, outputFrequency=50, messages=False, nThreads=2)
res[1] = CoGAPS(csv_path, outputFrequency=50, messages=False, nThreads=6)
res[2] = CoGAPS(csv_path, outputFrequency=50, messages=False, nThreads=12)
for r in res:
    assert(no_na_in_result(r) == True)



# test running with fixed matrix
print('\n ### Testing Fixed Matrix Run ###\n')
nPat = 3
fixedA = np.random.uniform(1, 10, (adata.X.shape[0], nPat))
fixedP = np.random.uniform(1, 10, (adata.X.shape[1], nPat))
params = CoParams(csv_path)
params.setFixedPatterns(fixedA, "A")
setParams(params, {'nIterations': 100,
                    'seed': 42,
                    'nPatterns': nPat})
res = CoGAPS(csv_path, params, outputFrequency=100, messages=False)

assert(toNumpy(res['GapsResult'].Amean).shape == fixedA.shape)
for i in range(fixedA.shape[1]):
    fixedA[:,i] = fixedA[:,i] * (toNumpy(res['GapsResult'].Amean)[0,i] / fixedA[0,i])

assert(np.allclose(fixedA, toNumpy(res['GapsResult'].Amean), rtol=1e-3))


params = CoParams(csv_path)
params.setFixedPatterns(fixedP, "P")
setParams(params, {'nIterations': 100,
                    'seed': 42,
                    'nPatterns': nPat})
res = CoGAPS(csv_path, params, outputFrequency=100, messages=False)

assert(toNumpy(res['GapsResult'].Pmean).shape == fixedP.shape)
for i in range(fixedP.shape[1]):
    fixedP[:,i] = fixedP[:,i] * (toNumpy(res['GapsResult'].Pmean)[0,i] / fixedP[0,i])

assert(np.allclose(fixedP, toNumpy(res['GapsResult'].Pmean), rtol=1e-3))


# testing that None gets converted to NULL for distributed
print('\n ### Testing Distributed is None Run ###\n')
params = CoParams(csv_path)
setParams(params, {'nIterations': 100,
                    'seed': 42,
                    'nPatterns': 3,
                    'distributed': None})
res = CoGAPS(csv_path, params, outputFrequency=100, messages=False)

'''
# test using saved parameters in file
# use pickle ??
# pickle isn't compatible with c++ objects (GapsParameters)
# if needed can add support for this later
matP = getSampleFactors(res['anndata'])
params = CoParams(csv_path)
setParams(params, {'nPatterns': matP.shape[1],
                    'nIterations': 175,
                    'seed': 42,
                    'useSparseOptimization': True,
                    'distributed': "genome-wide",
                    'explicitSets': [np.arange(0,201), np.arange(201, 401), np.arange(401, 601), np.arange(601,801), np.arange(801,1001)]})
params.setDistributedParams(nSets=5, cut=matP.shape[1] + 1)
params.setFixedPatterns(matP, "P")


with open('temp_params.pkl', 'wb') as outp:
    pickle.dump(params, outp)

with open('temp_params.pkl', 'rb') as inp:
    temp_params = pickle.load(inp)

res1 = CoGAPS(csv_path, params)
res2 = CoGAPS(csv_path, temp_params)
'''