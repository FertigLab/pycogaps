from PyCoGAPS import *
import time

# replace with the path to your data, or use this provided example
path = "src/CoGAPS/inst/extdata/retina_subset_1.h5" 

params = CoParams(path, hdfKey="counts", transposeData=True)

setParams(params, {"seed": 0,
                    "nIterations": 10000,
                    "nPatterns": 3,
                    "useSparseOptimization": True
                    })


start = time.time()
CoGAPS(path, params, transposeData=True)
end = time.time()
print(end - start)
