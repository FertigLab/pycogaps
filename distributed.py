import warnings

from pandas.io.clipboard import paste

import PyCoGAPS
import subset_data
import multiprocessing
import helper_functions
import numpy as np
import pandas as pd
import itertools

def distributedCoGAPS(path, params, uncertainty=None):
    data = helper_functions.toAnndata(path)
    sets = subset_data.createSets(data, params)
    if min(map(len, sets)) < params.gaps.nPatterns:
        warnings.warn("Data subset dimension less than nPatterns. Aborting.")
        return 1

    PyCoGAPS.setParams(params, {'checkpointOutFile': ""})

    if params.coparams["fixedPatterns"] is None:
        print("Running Across Subsets...\n\n")
        with multiprocessing.get_context("spawn").Pool(processes=len(sets)) as pool:
            # make a list of parameters for each function call so they can easily be mapped to processes
            paramlst = []
            for i in range(len(sets)):
                paramlst.append([path, params, i, sets[i], None])

            result = pool.map(callInternalCoGAPS, paramlst)
            pool.close()
            pool.join()

            if params.coparams['distributed'] == "genome-wide":
                # unmatched = list(map(lambda x: np.array(x["GapsResult"].Pmean), result))
                unmatched = []
                print("length of result", len(result))
                for r in result:
                    unmatched.append(np.array(r["GapsResult"].Pmean))
            else:
                unmatched = map(lambda x: np.array(x["GapsResult"].Amean), result)
            print("Matching patterns across subsets...\n")
            matched = findConsensusMatrix(unmatched, params)
    else:
        matched = params.gaps.fixedPatterns

    return matched

    # params.gaps.nPatterns = matched.consensus.shape[1]
    # params.gaps.fixedPatterns = matched.consensus
    # if params.coparams["distributed"] == "genome-wide":
    #     params.gaps.whichMatrixFixed = "P"
    # else:
    #     params.gaps.whichMatrixFixed = "A"
    #
    # print("Running final stage...")
    # with multiprocessing.get_context("spawn").Pool(processes=len(sets)) as pool:
    #     paramlst = []
    #     for i in range(len(sets)):
    #         paramlst.append([path, params, i, sets[i], None])
    #     finalresult = pool.map(callInternalCoGAPS, paramlst)
    #     pool.close()
    #     pool.join()
    #
    # fullresult = stitchTogether(finalresult, params, sets)
    #
    # # add diagnostics...
    # fullresult.diagnostics.firstPass = result
    #
    # return fullresult


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


def findConsensusMatrix(unmatched, params):
    print("not yet implemented")
    allpatterns = pd.DataFrame(np.hstack(unmatched))
    comb = expandgrid(range(params.coparams["nSets"]), range(params.gaps.nPatterns))
    comb = list(comb.values())
    comb = pd.DataFrame(comb)
    comb = pd.DataFrame.transpose(comb)
    comb = comb.to_numpy()
    names = []
    for i in range(comb.shape[0]):
        names.append(str(comb[i, 0]+1) + "." + str(comb[i, 1]+1))
    allpatterns.columns = names
    print(allpatterns)
    return patternMatch()
    # allPatterns < - do.call(cbind, unmatchedPatterns)
    # comb < - expand.grid(1: gapsParams @ nSets, 1: gapsParams @ nPatterns)
    # colnames(allPatterns) < - paste(comb[, 1], comb[, 2], sep = ".")
    # return (patternMatch(allPatterns, gapsParams))
    return 1


def expandgrid(*itrs):
   product = list(itertools.product(*itrs))
   return {'Var{}'.format(i+1):[x[i] for x in product] for i in range(len(itrs))}


def patternMatch():
    print("not yet implemented")


def stitchTogether(result, params, sets):
    print("not yet implemented")