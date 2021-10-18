import sys
sys.path.append(".") # Adds higher directory to python modules path.

from PyCoGAPS import *

path = './src/CoGAPS/inst/extdata/retina_subset_1.h5'

# if input is an hdf file, then need to create params object, and set key name when initializing
params = CoParams(path, hdfKey='counts')

result = CoGAPS(path, params=params)
