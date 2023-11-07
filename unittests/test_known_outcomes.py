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
        self.res = CoGAPS(modsim, params)

    def test(self):
        #both recorded based on 4e7702f
        self.assertEqual(4.033180764186739, np.mean(self.res.obs))
        self.assertEqual(0.2936586886479366, np.mean(self.res.var))


class ModSimDistributedRuns(unittest.TestCase):
    #runs on ModSim - the toy dataset
    PATH = 'data/ModSimData.txt'

    def setUp(self):
        modsim = sc.read_text(self.PATH)
        params = CoParams(path=self.PATH)
        setParams(params, {'distributed': 'genome-wide'})
        params.setDistributedParams(nSets=7)
        self.res_dist = CoGAPS(modsim, params)

    def test(self):
        #dummy test to check that distributed has ran
        return True



if __name__=='__main__':
    unittest.main()
        
        