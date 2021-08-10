from PyCoGAPS import *
from subset_data import *
import warnings
import pathos
import dill
import multiprocessing
import main
from multiprocessing import Pool, TimeoutError
import time
import os


def callInternalCoGAPS(data, allParams, uncertainty=None, subsetIndices=None, workerID=None):
    print("in callinternalcogaps function, name is", __name__)
    print("assigning genome wide")
    subsetIndices = subsetIndices[workerID]
    genomewide = allParams.coparams["distributed"] == "genome-wide"

    if genomewide and allParams.coparams["geneNames"] is not None:
        allParams.coparams["geneNames"] = allParams.coparams["geneNames"][subsetIndices]
    else:
        if allParams.coparams["sampleNames"] is not None:
            allParams.coparams["sampleNames"] = allParams.coparams["sampleNames"][subsetIndices]
    print("callinternal")
    allParams.coparams["subsetIndices"] = subsetIndices
    if genomewide:
        allParams.coparams["subsetDim"] = 1
    else:
        allParams.coparams["subsetDim"] = 2

    allParams.gaps.workerID = workerID
    allParams.gaps.asynchronousUpdates = None
    allParams.gaps.nThreads = 1

    return CoGAPS(data, allParams, uncertainty)


def callback():
    print("in callback")


def handleResult():
    print("handler was called")
    return "this is the hidden treasure baby"


def goofin(mat):
    print("in goofin")
    return "flomillishit"


def distributedCoGAPS(data, allParams, uncertainty=None):
    print("in distributed cogaps function, name is", __name__)
    # step 0: if the user passed in a path, fetch data and convert to anndata object
    if isinstance(data, str):
        path = data
        data = toAnndata(data)
    # step 1: randomly break up the data
    sets = createSets(data, allParams)
    print("number of sets is", len(sets))
    if min(len(s) for s in sets) < allParams.gaps.nPatterns:
        warnings.warn("data subset dimension less than nPatterns--terminating execution")
        return

    if allParams.coparams["fixedPatterns"] is None:
        print("beginning multithreading, name is...", __name__)
        # start 4 worker processes
        i=0
        with multiprocessing.get_context("spawn").Pool(processes=len(sets)) as pool:
            # for i in range(len(sets)):
            #     print("i=",i)
            # pool.apply_async(main.callinternal, args=(path, allParams, uncertainty, sets, 1), callback=callback)
            print("started the pool")
            # result = pool.apply_async(callInternalCoGAPS, args=(path, allParams, uncertainty, sets[1], 1))
            # result = pool.apply_async(CoGAPS, args=(path, allParams, uncertainty))
            matrix = pycogaps.Matrix(5,5)
            print("made a matrix")
            result = pool.apply_async(goofin, args=[matrix])
            # print("async call returned. now trying to get()")
            # result = result.get()
            # print("get() call returned")
            pool.close()
            print("closed the pool")
            pool.join()
            print("joined the pool")
    print("out of the loop baby")
    return result


def findConsensusMatrix(unmatchedPatterns, gapsParams):
    print("Not yet implemented")
    return


def patternMatch(allPatterns, gapsParams):
    print("Not yet implemented")
    return


def corrToMeanPattern(cluster):
    print("Not yet implemented")
    return


def corcut(allPatterns, cut, minNS):
    print("Not yet implemented")
    return


def stitchTogether(result, allParams, sets):
    print("Not yet implemented")
    return

def runDistributed(data, allParams, uncertainty=None):
    if __name__ == "__main__":
        return distributedCoGAPS(data, allParams, uncertainty=None)