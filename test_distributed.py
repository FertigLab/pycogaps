import distributed_cogaps
import main
from PyCoGAPS import *


print("testing, name is", __name__)
if __name__ == "__main__":
    # replace with the path to your data, or use this provided example
    path = "./data/GIST.csv"

    # create a CoParams object and set desired parameters
    params = CoParams(path=path)
    # set distributed params, annotation weights, fixed patterns by calling specific methods
    params.setDistributedParams(nSets=2)
    params.setDistributedParams(sampleNames=['IM00', 'IM02', 'IM04', 'IM06', 'IM09', 'IM12', 'IM18', 'IM24', 'IM48'])
    params.setDistributedParams(subsetIndices=1)
    params.setFixedPatterns(None, None)
    result = distributed_cogaps.distributedCoGAPS(path, params)
    print("returned for real, attempting to get result")
    trueres = result.get()
    print(trueres)
