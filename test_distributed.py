import distcaller
from PyCoGAPS import *
# replace with the path to your data, or use this provided example
path = "./data/GIST.csv"

# create a CoParams object and set desired parameters
params = CoParams(path=path)
# set distributed params, annotation weights, fixed patterns by calling specific methods
params.setDistributedParams(nSets=2)
params.setFixedPatterns(None, None)

print("testing, name is", __name__)
distcaller.calldistributed(path, params)
#
# # dr = runDistributed(path, params)
# if __name__ == "__main__":
#     distresult = runDistributedCoGAPS(path, params)