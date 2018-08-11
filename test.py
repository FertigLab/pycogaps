import sys
sys.path.insert(0,"./src")
from Cogaps import Cogaps

result = Cogaps(load=False, dataPath="data/GIST.tsv", numPatterns=5, maxIterations=50000)
result.graphPmean()
result.graphPatternMarkerStats(0)
