import warnings

import PyCoGAPS
import subset_data
import multiprocessing
import helper_functions
import numpy as np


def distributedCoGAPS(path, params, uncertainty=None):
    data = helper_functions.toAnndata(path)
    sets = subset_data.createSets(data, params)
    if min(map(len, sets)) < params.gaps.nPatterns:
        warnings.warn("Data subset dimension less than nPatterns. Aborting.")
        return 1

    PyCoGAPS.setParams(params, {'checkpointOutFile':""})

    if params.coparams["fixedPatterns"] is None:
        print("Running Across Subsets...\n\n")
        with multiprocessing.get_context("spawn").Pool(processes=len(sets)) as pool:
            # make a list of parameters for each function call so they can easily be mapped to processes
            paramlst = []
            for i in range(len(sets)):
                paramlst.append([path, params, i, sets[i], None])

            result = pool.map(callInternalCoGAPS, paramlst)
            pool.close()
            print("closed the pool")
            pool.join()
            print("joined the pool")
            # if params.coparams['distributed'] == "genome-wide":
            #     unmatched = result.get().Pmean
            # else:
            #     unmatched = result.get().Amean

    else:
        print("not implemented")
    return result


def callInternalCoGAPS(paramlst):
    # take out parameters passed as a list to the worker process
    path = paramlst[0]
    params = paramlst[1]
    workerID = paramlst[2]
    subsetIndices = paramlst[3]
    uncertainty = paramlst[4]
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
