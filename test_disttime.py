from PyCoGAPS import *
import time
import numpy as np
import h5py
# #
# path = "src/CoGAPS/inst/extdata/retina_subset_1.h5"
# np.savetxt("src/CoGAPS/inst/extdata/retina_subset_1.csv", h5py.File(path)['counts'], '%g', ',')
# replace with the path to your data, or use this provided example
path = "./data/GIST.csv"
# df = pd.read_hdf(path)
# df.to_csv("src/CoGAPS/inst/extdata/retina_subset_1.csv", index=False)
# csvpath = "src/CoGAPS/inst/extdata/retina_subset_1.csv"
params = CoParams(path, hdfKey="counts")

setParams(params, {"seed": 0,
                    "nIterations": 100,
                    "nPatterns": 3,
                    "useSparseOptimization": False,
                    "hdfKey": "counts",
                    })


start = time.time()
if __name__ == '__main__':
    params.setDistributedParams(nSets=3, minNS=1, cut=3)
    result = distributedCoGAPS(path, params, None)
    print("RUN FINISHED")
    plotResiduals(result)
    binaryA(result, threshold=3)
    plotPatternMarkers(result, colDendrogram=True)
    plotPatternMarkers(result, rowDendrogram=True)
    # print(result)
# CoGAPS(path, params, transposeData=True)
end = time.time()
print(end - start)