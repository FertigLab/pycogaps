'''
a script for testing the python interface to cogaps
jeanette johnson 6/14/21
before running:
navigate to your pycogaps direcotry and run 'pip install .'
'''
import pycogaps

result = pycogaps.runCogaps("./data/GIST.csv")
print("AMEAN: ", result.Amean)
print("chisqHistory: ", result.chisqHistory)
print(pycogaps.getBuildReport())
print(pycogaps.isCheckpointsEnabled())
print(pycogaps.isCompiledWithOpenMPSupport())
print(pycogaps.getFileInfo("./data/GIST.csv"))
