'''
a script for testing the python interface to cogaps
jeanette johnson 6/14/21
before running:
navigate to your pycogaps direcotry and run 'pip install .'
'''
import pycogaps

pycogaps.runCogaps("./data/GIST.csv")
pycogaps.runCPPTests() 
