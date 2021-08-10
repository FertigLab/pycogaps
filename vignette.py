from PyCoGAPS import *

path = "/Users/jeanettejohnson/Desktop/FertigLab/CensusImmune-CordBlood-10x_cell_type_2020-03-12.csv"

prm = GapsParameters(path)
prm.print()  # c++ method to display all parameter values
CoGAPS(path, prm)
cogapsrunresult = CoGAPS(path)