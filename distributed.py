import PyCoGAPS
from PyCoGAPS import *
import subset_data
import multiprocessing
import helper_functions
import numpy as np


def distributedCoGAPS(path, params, uncertainty=None):
    data = helper_functions.toAnndata(path)
    sets = subset_data.createSets(data, params)
    with multiprocessing.get_context("spawn").Pool(processes=len(sets)) as pool:
        m = pycogaps.Matrix(4, 4)
        result = pool.apply_async(callInternalCoGAPS, args=[path, params])
        pool.close()
        print("closed the pool")
        pool.join()
        print("joined the pool")
    return result


def callInternalCoGAPS(path, params, uncertainty=None, subsetIndices=None, workerID=None):
    gapsresult = PyCoGAPS.CoGAPS(path, params)
    return gapsresult


def callback(mat, params):
    return "returned from callback"