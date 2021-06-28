'''
a script for testing the python interface to cogaps
jeanette johnson 6/14/21
before running:
navigate to your pycogaps direcotry and run 'pip install .'
'''
from PyCoGAPS import * # gonna try to only use things from this module

path = "./data/GIST.csv"

# print("-- Testing python bindings --\n")
# pycogaps.runCPPTests()

print("-- Passing params object into runCogaps function --\n")
prm = GapsParameters(path)
prm.print()  # c++ method to display all parameter values
CoGAPS(path, prm)
CoGAPS(path, prm)
result = CoGAPS(path)
print("AMEAN: ", result.Amean)
print("chisqHistory: ", result.chisqHistory)
print(getBuildReport())
print(isCheckpointsEnabled())
print(isCompiledWithOpenMPSupport())
print(getFileInfo(path))
