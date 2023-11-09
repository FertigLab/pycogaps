import unittest
from PyCoGAPS.parameters import *
from PyCoGAPS.pycogaps_main import CoGAPS
import scanpy as sc

import logging as log

#tests that run different run types with default params

class ModSimDefault(unittest.TestCase):
    #runs on ModSim - the toy dataset
    PATH = 'data/ModSimData.txt'

    def setUp(self):
        modsim = sc.read_text(self.PATH)
        params = CoParams(path=self.PATH)
        self.res_dist = CoGAPS(modsim, params)

    def test(self):
        #dummy test to check that test has ran
        return True

class ModSimDistributedDefault(unittest.TestCase):
    #runs on ModSim - the toy dataset
    PATH = 'data/ModSimData.txt'

    def setUp(self):
        modsim = sc.read_text(self.PATH)
        params = CoParams(path=self.PATH)
        setParams(params, {'distributed': 'genome-wide'})
        self.res_dist = CoGAPS(modsim, params)

    def test(self):
        #dummy test to check that test has ran
        return True
    

if __name__=='__main__':
    unittest.main()