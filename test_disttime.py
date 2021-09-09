from PyCoGAPS import *
import time

# replace with the path to your data, or use this provided example
path = "src/CoGAPS/inst/extdata/retina_subset_1.h5"
# path = "./data/GIST.csv"

params = CoParams(path, hdfKey="counts", transposeData=False)

setParams(params, {"seed": 0,
                    "nIterations": 100,
                    "nPatterns": 3,
                    "useSparseOptimization": False,
                   "transposeData": False,
                    "hdfKey": "counts"
                    })


start = time.time()
if __name__ == '__main__':
    params.setDistributedParams(nSets=3, minNS=1, cut=3)
    result = distributedCoGAPS(path, params, None)
    print("RUN FINISHED")
    # print(result)
# CoGAPS(path, params, transposeData=True)
end = time.time()
print(end - start)