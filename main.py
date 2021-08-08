from distributed_cogaps import *


def calldistributed(path, params):
    print("name is ", __name__)
    if __name__ == "main":
        return distributedCoGAPS(path, params)
    else:
        return


def callinternal(data, allParams, uncertainty=None, subsetIndices=None, workerID=None):
    print("name is ", __name__)
    if __name__ == "main":
        return callInternalCoGAPS(data, allParams, uncertainty, subsetIndices, workerID)
        # return callInternalCoGAPS(data, allParams)
    else:
        return