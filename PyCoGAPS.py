import time
import math
import pycogaps


# contains wrappers for every C++ bound function
# TODO: getters, setters, input sanitation

def CoGAPS(path, params=None):
    if params is None:
        return pycogaps.runCogaps(path)
    else:
        return pycogaps.runCogapsWithParams(path, params)


def GapsParameters(path):
    return pycogaps.GapsParameters(path)


def getBuildReport():
    return pycogaps.getBuildReport()


def isCheckpointsEnabled():
    return pycogaps.isCheckpointsEnabled()


def isCompiledWithOpenMPSupport():
    return pycogaps.isCompiledWithOpenMPSupport()


def getFileInfo(path):
    return pycogaps.getFileInfo(path)


def current_milli_time():
    return round(time.time() * 1000)


# class CogapsParams
# constructor has default values for each parameter
def setParam(paramobj, whichParam, value):
    if whichParam == "nPatterns":
        paramobj.nPatterns = value
    elif whichParam == "nIterations":
        paramobj.nIterations = value
    elif whichParam == "alpha":
        paramobj.alphaA = value
        paramobj.alphaP = value
    elif whichParam == "maxGibbsMass":
        paramobj.maxGibbsMassA = value
        paramobj.maxGibbsMassP = value
    elif whichParam in ("nSets", "cut", "minNS", "maxNS"):
        print("please set \'", whichParam, "\' with setDistributedParams")
        return
    elif whichParam in ("samplingAnnotation", "samplingWeight"):
        print("please set \'", whichParam, "\' with setAnnotationWeights")
        return
    elif whichParam in ("fixedPatterns", "whichMatrixFixed"):
        print("please set \'", whichParam, "\' with setFixedPatterns")
        return
    # elif whichParam == "nPatterns":
    #     paramobj.nPatterns = value
    #     paramobj.cut = min(paramobj.cut, paramobj.nPatterns)
    # elif whichParam == "distributed":
    #     if value == "none":
    #         paramobj.distributed = None
    #     else:
    #         paramobj.distributed = value
    elif whichParam in "singleCell":
        print(whichParam, " has been deprecated, this parameter will be ignored")
        return
    else:
        setattr(paramobj, whichParam, value)
    return paramobj


def show(paramobj):
    paramobj.print()
