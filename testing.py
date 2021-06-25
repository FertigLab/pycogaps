'''
a script for testing the python interface to cogaps
jeanette johnson 6/14/21
before running:
navigate to your pycogaps direcotry and run 'pip install .'
'''
import pycogaps

path = "./data/GIST.csv"

# print("-- Testing python bindings --\n")
# pycogaps.runCPPTests()

print("-- Passing params object into runCogaps function --\n")
prm = pycogaps.GapsParameters(path)
prm.print()  # c++ method to display all parameter values
pycogaps.runCogaps(path, prm)
result = pycogaps.runCogaps("./data/GIST.csv")
print("AMEAN: ", result.Amean)
print("chisqHistory: ", result.chisqHistory)
print(pycogaps.getBuildReport())
print(pycogaps.isCheckpointsEnabled())
print(pycogaps.isCompiledWithOpenMPSupport())
print(pycogaps.getFileInfo("./data/GIST.csv"))
