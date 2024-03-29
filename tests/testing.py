'''
a script for testing the python interface to cogaps
jeanette johnson 6/14/21
before running:
navigate to your pycogaps direcotry and run 'pip install .'
'''
import sys
sys.path.append(".") # Adds higher directory to python modules path.

from PyCoGAPS import *  # gonna try to only use things from this module

path = "./data/GIST.csv"

print("-- Passing params object into runCogaps function --\n")
prm = CoParams(path)
prm.printAllParams()  # c++ method to display all parameter values
CoGAPS(path, prm)
cogapsrunresult = CoGAPS(path)
result = cogapsrunresult['GapsResult']
anndata = cogapsrunresult['anndata']
plotResiduals(cogapsrunresult)
plotPatternMarkers(cogapsrunresult, legend_pos=None)
plotPatternMarkers(cogapsrunresult, samplePalette=["green", "teal", "red", "violet", "crimson", "antiquewhite", "lightblue", "hotpink", "orange"], patternPalette=["pink", "teal", "gold"],
                   legend_pos=None)
binaryA(cogapsrunresult, threshold=3)
binaryA(cogapsrunresult, threshold=3, cluster=True)
print("AMEAN: ", result.Amean)
print("chisqHistory: ", result.chisqHistory)
print(getBuildReport())
print(isCheckpointsEnabled())
print(isCompiledWithOpenMPSupport())
print(getFileInfo(path))

print("--Testing CogapsResult Object\n")
print("calling show(result)\n")
show(anndata)
print("calling plot(result)\n")
plot(cogapsrunresult)

calcZP = calcZ(anndata, "sampleFactors")
print(calcZP)
calcZA = calcZ(anndata, "featureLoadings")
print(calcZA)
getVersion()

print("~~~~~~~~~~~~~ testing CoGAPS Stat Functions ~~~~~~~~~~~~~~")

dict = calcCoGAPSStat(cogapsrunresult, sets=['Hs.101174', 'Hs.1012'])
print(dict)

outStats = calcGeneGSStat(cogapsrunresult, GStoGenes=['Hs.101174', 'Hs.1012'], numPerm=1000)
print(outStats)

finalStats = computeGeneGSProb(cogapsrunresult, GStoGenes=['Hs.101174', 'Hs.1012'])
print(finalStats)