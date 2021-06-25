'''
a script for testing the python interface to cogaps
jeanette johnson 6/14/21
before running:
navigate to your pycogaps direcotry and run 'pip install .'
'''
import pycogaps
import CoGAPSParameters
path = "./data/GIST.csv"

print("-- Testing python bindings --\n")
# pycogaps.runCogaps("./data/GIST.csv")
pycogaps.runCPPTests()

# print("-- Testing CogapsParams object --\n")
# params = CoGAPSParameters.CogapsParams()
# params.show()

print("-- Passing params object into runCogaps function --\n")
prm = pycogaps.GapsParameters(path)
pycogaps.runCogaps(path, prm)



# print("-- Making pho --\n")
# largepho = pypho.pypho()
# pycogaps.testStructCasting(largepho)
