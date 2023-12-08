import unittest
from PyCoGAPS.parameters import *
from PyCoGAPS.pycogaps_main import CoGAPS
import scanpy as sc
import sys



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
        platform = sys.platform
        #recorded based on 4e7702f:
        mean_res_obs = {'darwin': 4.998108364973596,
                        'linux': 4.9981058453431855,
                        'win': 4.9981058453431855} #change once known

        self.assertAlmostEqual(mean_res_obs[platform], np.mean(self.res.obs))


if __name__=='__main__':
    unittest.main()
        
        