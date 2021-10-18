import sys
sys.path.append(".") # Adds higher directory to python modules path.

from PyCoGAPS import *

mtx_path = "./data/GIST.mtx" 

if isCheckpointsEnabled():
    params = CoParams(mtx_path)
    setParams(params, {"nIterations": 100,
                        "seed": 42})
    res1 = CoGAPS(mtx_path, params, checkpointInterval=100, checkpointOutFile="test.out", messages=False)
    res2 = CoGAPS(mtx_path, params, checkpointInFile="test.out", messages=False)

    res1 = res1['GapsResult']
    res2 = res2['GapsResult']

    assert(np.allclose(toNumpy(res1.Amean), toNumpy(res2.Amean), rtol=0.1))
    assert(np.allclose(toNumpy(res1.Pmean), toNumpy(res2.Pmean), rtol=0.1))