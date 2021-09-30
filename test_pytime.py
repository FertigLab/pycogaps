from PyCoGAPS import *
import time

# replace with the path to your data, or use this provided example
path = "data/GSE98638_HCC.TCell.S5063.count.txt"

params = CoParams(path)

setParams(params, {"seed": 0,
                    "nIterations": 10000,
                    "nPatterns": 10,
                    "useSparseOptimization": True
                    })


start = time.time()
result = CoGAPS(path, params)
end = time.time()
print(end - start)
# plotPatternMarkers(result)
