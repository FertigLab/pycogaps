from PyCoGAPS import *

def results_equal(res1, res2):
    res1 = res1['GapsResult']
    res2 = res2['GapsResult']
    assert(np.allclose(toNumpy(res1.Amean), toNumpy(res2.Amean), rtol=0.1))
    assert(np.allclose(toNumpy(res1.Asd), toNumpy(res2.Asd), rtol=0.1))
    assert(np.allclose(toNumpy(res1.Pmean), toNumpy(res2.Pmean), rtol=0.1))
    assert(np.allclose(toNumpy(res1.Psd), toNumpy(res2.Psd), rtol=0.1))
    assert(res1.atomHistoryA == res2.atomHistoryA)
    assert(res1.atomHistoryP == res2.atomHistoryP)

mtx_path = "./data/GIST.mtx" 

# standard cogaps
params = CoParams(mtx_path)
setParams(params, {"nIterations": 100,
                    "seed": 42})
res1 = CoGAPS(mtx_path, params, outputFrequency=10, messages=False)
res2 = CoGAPS(mtx_path, params, outputFrequency=10, messages=False)
results_equal(res1, res2)

'''
# TODO: run once implement distributed cogaps
# distributed cogaps
params = CoParams(mtx_path)
setParams(params, {"nIterations": 100,
                    "seed": 42,
                    "distributed": "genome-wide"})
res1 = CoGAPS(mtx_path, params, outputFrequency=10, messages=False)
res2 = CoGAPS(mtx_path, params, outputFrequency=10, messages=False)
results_equal(res1, res2)
'''

# multiple threads, dense sampler
params = CoParams(mtx_path)
setParams(params, {"nIterations": 100,
                    "seed": 42,
                    "useSparseOptimization": False})
res1 = CoGAPS(mtx_path, params, outputFrequency=10, messages=False, nThreads=1)
res2 = CoGAPS(mtx_path, params, outputFrequency=10, messages=False, nThreads =3)
res3 = CoGAPS(mtx_path, params, outputFrequency=10, messages=False, nThreads =6)
results_equal(res1, res2)
results_equal(res1, res3)
results_equal(res2, res3)

# multiple threads, sparse sampler
params = CoParams(mtx_path)
setParams(params, {"nIterations": 100,
                    "seed": 42,
                    "useSparseOptimization": True})
res1 = CoGAPS(mtx_path, params, outputFrequency=10, messages=False, nThreads=1)
res2 = CoGAPS(mtx_path, params, outputFrequency=10, messages=False, nThreads =3)
res3 = CoGAPS(mtx_path, params, outputFrequency=10, messages=False, nThreads =6)
results_equal(res1, res2)
results_equal(res1, res3)
results_equal(res2, res3)