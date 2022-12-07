from PyCoGAPS.parameters import *
from PyCoGAPS.pycogaps_main import CoGAPS
import scanpy as sc

modsimpath = "/Users/jeanette/fertiglab/pycogaps/data/ModSimData.txt"
modsimbasespath = "/Users/jeanette/fertiglab/pycogaps/data/ModSimBases.txt"

modsim = sc.read_text(modsimpath)
modsimbases = sc.read_text(modsimbasespath)

params = CoParams(path=modsimpath)
params.printParams()

setParams(params, {
    'nIterations': 50000,
    'seed': 42,
    'nPatterns': 3
})

# many people find it helpful to time cogaps runs
start = time.time()
# command that calls CoGAPS
# TIMING: on ModSim data, this only should take about 3 sec
modsimresult = CoGAPS(modsim, params)
end = time.time()
print("TIME:", end - start)

# always write cogaps result to disc before doing anything else!
modsimresult.write("/Users/jeanette/fertiglab/pycogaps/data/ModSimPyCoGAPSResult.h5ad")





