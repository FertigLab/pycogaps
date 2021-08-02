from distributed_cogaps import *

# replace with the path to your data, or use this provided example
path = "./data/GIST.csv"

# create a CoParams object and set desired parameters
params = CoParams(path=path)

callInternalCoGAPS(path, params, None, "IM00", 4)
# set distributed params, annotation weights, fixed patterns by calling specific methods
params.setDistributedParams(nSets=5)

result = distributedCoGAPS(path, params, None)