import PyCoGAPS
import subset_data
import multiprocessing
import helper_functions
import numpy as np


def distributedCoGAPS(path, params, uncertainty=None):
    data = helper_functions.toAnndata(path)
    sets = subset_data.createSets(data, params)
    PyCoGAPS.setParams(params, {'checkpointOutFile':""})
    with multiprocessing.get_context("spawn").Pool(processes=len(sets)) as pool:
        m = PyCoGAPS.pycogaps.Matrix(4, 4)
        result = pool.apply_async(callInternalCoGAPS, args=[path, params, 1, sets, None])
        pool.close()
        print("closed the pool")
        pool.join()
        print("joined the pool")
    return result


def callInternalCoGAPS(path, params, workerID, subsetIndices=None, uncertainty=None):
    if subsetIndices is None:
        print("No subset indices provided; generating random sets...")
        adata = helper_functions.toAnndata(path)
        subsetIndices = subset_data.createSets(adata, params)
    if params.coparams['distributed'] == "genome-wide":
        genes = np.array(params.coparams['geneNames'])
        params.coparams['geneNames'] = np.take(genes, subsetIndices)
        params.coparams['subsetDim'] = 1
    else:
        samples = np.array(params.coparams['sampleNames'])
        params.coparams['sampleNames'] = np.take(samples, subsetIndices)
        params.coparams['subsetDim'] = 2

    params.coparams['subsetIndices'] = subsetIndices
    params.gaps.workerID = workerID
    params.gaps.asynchronousUpdates = False
    params.gaps.maxThreads = 1
    gapsresult = PyCoGAPS.CoGAPS(path, params, uncertainty)

    return gapsresult


def callback(mat, params):
    return "returned from callback"
