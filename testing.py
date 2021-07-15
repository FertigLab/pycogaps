'''
a script for testing the python interface to cogaps
jeanette johnson 6/14/21
before running:
navigate to your pycogaps direcotry and run 'pip install .'
'''
from PyCoGAPS import * # gonna try to only use things from this module

path = "./data/GIST.csv"

print("-- Passing params object into runCogaps function --\n")
prm = GapsParameters(path)
prm.print()  # c++ method to display all parameter values
result = CoGAPS(path, prm)
print("AMEAN: ", result.Amean)
print("chisqHistory: ", result.chisqHistory)
print(getBuildReport())
print(isCheckpointsEnabled())
print(isCompiledWithOpenMPSupport())
print(getFileInfo(path))

print("--Testing CogapsResult Object\n")
show(result)
plot(result)

calcZP = calcZ(result, "sampleFactors")
print(calcZP)
calcZA = calcZ(result, "featureLoadings")
print(calcZA)

