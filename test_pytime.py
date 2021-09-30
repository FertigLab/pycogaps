from PyCoGAPS import *
import time

# replace with the path to your data, or use this provided example
path = "data/GIST.csv"

params = CoParams(path)

setParams(params, {"seed": 0,
                    "nIterations": 10000,
                    "nPatterns": 3,
                    "useSparseOptimization": True
                    })


start = time.time()
result = CoGAPS(path, params)
end = time.time()
print(end - start)
# plotPatternMarkers(result)
