'''
a script for testing the python interface to cogaps
jeanette johnson 6/14/21
before running:
navigate to your pycogaps direcotry and run 'pip install .'
'''
from PyCoGAPS import *  # gonna try to only use things from this module

path = "./data/GIST.csv"

print("-- Passing params object into runCogaps function --\n")
# prm = GapsParameters(path)
# prm.print()  # c++ method to display all parameter values
# CoGAPS(path, prm)
cogapsrunresult = CoGAPS(path)
result = cogapsrunresult['GapsResult']
anndata = cogapsrunresult['anndata']
# binaryA(result, threshold=3)
originaldata = genfromtxt(path, delimiter=",")[1:, 1:]  # need to get rid of the first row and column
# plotResiduals(result, originaldata, None)
print("AMEAN: ", result.Amean)
print("chisqHistory: ", result.chisqHistory)
print(getBuildReport())
print(isCheckpointsEnabled())
print(isCompiledWithOpenMPSupport())
print(getFileInfo(path))

print("--Testing CogapsResult Object\n")
print("calling show(result)\n")
show(result)
print("calling plot(result)\n")
plot(result)

calcZP = calcZ(result, "sampleFactors")
print(calcZP)
calcZA = calcZ(result, "featureLoadings")
print(calcZA)
getVersion()
