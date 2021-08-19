from PyCoGAPS import *
import subset_data
import multiprocessing
import helper_functions


def distributedCoGAPS(path, params, uncertainty=None):
    data = helper_functions.toAnndata(path)
    sets = subset_data.createSets(data, params)
    with multiprocessing.get_context("spawn").Pool(processes=len(sets)) as pool:
        m = pycogaps.Matrix(4, 4)
        result = pool.apply_async(callback, args=[m])
        pool.close()
        print("closed the pool")
        pool.join()
        print("joined the pool")
    return result


def callInternalCoGAPS(path, params, uncertainty, subsetIndices, workerID):
    return "okokok"


def callback(mat):
    return "returned from callback"