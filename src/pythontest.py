import numpy as np
import Cogaps

my_data = np.genfromtxt('./Rpackage/inst/extdata/GIST.csv', delimiter=',')
my_data = np.delete(my_data, (0), axis=0)
my_data = np.delete(my_data, (0), axis=1)
print(my_data)
