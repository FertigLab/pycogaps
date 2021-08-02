from PyCoGAPS import *
import warnings
import multiprocessing


def callInternalCoGAPS(data, allParams, uncertainty=None, subsetIndices=None, workerID=None):
    genomewide = allParams.coparams["distributed"] == "genome-wide"

    if genomewide:
        allParams.coparams["geneNames"] = allParams.coparams["geneNames"][subsetIndices]
    else:
        allParams.coparams["sampleNames"] = allParams.coparams["sampleNames"][subsetIndices]

    allParams.coparams["subsetIndices"] = subsetIndices
    if genomewide:
        allParams.coparams["subsetDim"] = 1
    else:
        allParams.coparams["subsetDim"] = 2

    allParams.gaps.workerID = workerID
    allParams.gaps.asynchronousUpdates = None
    allParams.gaps.nThreads = 1

    return CoGAPS(data, allParams, uncertainty)


def distributedCoGAPS (data, allParams, uncertainty):
    print("Not yet implemented")

    # step 1: randomly break up the data
    sets = createSets(data, allParams)
    if min(len(s) for s in sets) < allParams.nPatterns:
        warnings.warn("data subset dimension less than nPatterns--terminating execution")
        return

    processes = []
    if allParams.gaps.fixedPatterns is None:
        for i in range(len(sets)):
            p = multiprocessing.Process(target=callInternalCoGAPS, args=(data, allParams, uncertainty, sets[[i]], i))
            processes.append(p)
            p.start()
        for process in processes:
            process.join()
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
