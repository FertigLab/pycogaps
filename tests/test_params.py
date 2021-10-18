import sys
sys.path.append(".") # Adds higher directory to python modules path.

from PyCoGAPS import *


# replace with the path to your data, or use this provided example
path = "./data/GIST.csv" 

# run CoGAPS on your dataset
result = CoGAPS(path)

# create a CoParams object and set desired parameters
params = CoParams(path)
setParam(params, "maxThreads", 4) 
# and/or:
setParams(params, {
            'printMessages': True,
            'maxThreads': 4
        })

# set distributed params, annotation weights, fixed patterns by calling specific methods
params.setDistributedParams(nSets=5)

result = CoGAPS(path, params)


