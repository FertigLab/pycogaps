import unittest
from PyCoGAPS.parameters import *
from PyCoGAPS.pycogaps_main import CoGAPS
import scanpy as sc

import logging as log



# Tests that run pycogaps with known conditions and compare run outcomes


class ModSimDefaultCase(unittest.TestCase):
    #runs on ModSim - the toy dataset
    PATH = 'data/ModSimData.txt'

    def setUp(self):
        modsim = sc.read_text(self.PATH)
        params = CoParams(path=self.PATH)
        setParams(params, {'seed': 0})
        self.res = CoGAPS(modsim, params)

    def test(self):
        #both recorded based on 4e7702f
        self.assertEqual(4.033180764186739, np.mean(self.res.obs))
        self.assertEqual(0.2936586886479366, np.mean(self.res.var))


if __name__=='__main__':
    unittest.main()
        
        