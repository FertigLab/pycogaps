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
    genomewide = allParams.coparams["distributed"] == "genome-wide"

    # if genomewide:
    #     allParams.coparams["geneNames"] = allParams.coparams["geneNames"][subsetIndices]
    # else:
    #     allParams.coparams["sampleNames"] = np.array(allParams.coparams["sampleNames"])[subsetIndices]

    allParams.coparams["subsetIndices"] = subsetIndices
    if genomewide:
        allParams.coparams["subsetDim"] = 1
    else:
        allParams.coparams["subsetDim"] = 2

    allParams.gaps.workerID = workerID
    allParams.gaps.asynchronousUpdates = None
    # allParams.gaps.nThreads = 1

    return CoGAPS(data, allParams, uncertainty)


def add(a,b):
    print("in add fxn")
    return a+b


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
        with Pool(processes=len(sets)) as pool:
            for i in range(len(sets)):
            # print(pool.apply(main.callinternal, args=(data, allParams, uncertainty, sets[1], 1)))
                print("run", i, "result is", pool.apply(main.callinternal, args=(path, allParams, uncertainty, sets[1], 1)))
                # if __name__ == '__main__':
                # print("in multithreading loop: run", i)
                # p = multiprocessing.Process(target=distcaller.callinternal, args=(data, allParams, uncertainty, sets[i], i))
                # processes.append(p)
                # print(p.name, "created")
                # print("starting", p.name, "with id", p.pid)
                #
                # print("returned from start() function for", p.name)
            # for process in processes:
            #     print("joining", p.name)
            #     process.join()
            #     print("joined", p.name)
    return


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