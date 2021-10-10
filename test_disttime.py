from PyCoGAPS import *
import time
import numpy as np
import h5py
# #
path = "data/GIST.csv"
# np.savetxt("src/CoGAPS/inst/extdata/retina_subset_1.csv", h5py.File(path)['counts'], '%g', ',')
# replace with the path to your data, or use this provided example
# path = "./data/GIST.csv"
# df = pd.read_hdf(path)
# df.to_csv("src/CoGAPS/inst/extdata/retina_subset_1.csv", index=False)
# csvpath = "src/CoGAPS/inst/extdata/retina_subset_1.csv"
params = CoParams(path)

setParams(params, {"seed": 0,
                    "nIterations": 10000,
                    "nPatterns": 10,
                    "useSparseOptimization": True,
                    "hdfKey": "counts",
                   "hdfColKey": "geneNames",
                   "hdfRowKey": "cellNames"
                    })


start = time.time()
if __name__ == '__main__':
    params.setDistributedParams()
    result = distributedCoGAPS(path, params, None)
end = time.time()
print("TIME:", end - start)
    # plotResiduals(result)
    # binaryA(result, threshold=3)
    # plotPatternMarkers(result, colDendrogram=True)
    # plotPatternMarkers(result, rowDendrogram=True)
    # print(result)
# CoGAPS(path, params, transposeData=True)
# Subset data
