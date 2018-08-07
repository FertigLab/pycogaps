import sys
sys.path.insert(0,"./src")
from Cogaps import Cogaps

result = Cogaps("data/GIST.tsv")
result.graphPmean()
