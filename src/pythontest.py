import numpy as np
import Cogaps

ones = np.ones((4,2))
twos = np.full((3,5), 2)
Cogaps.test_do_nothing()
newtuple = Cogaps.test_matrices(ones, twos)
print(newtuple[0])
print(newtuple[1])
