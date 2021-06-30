import time
import math
import pycogaps
import helper_functions


def CoGAPS(path, params=None, nThreads=1, messages=True,
           outputFrequency=1000, uncertainty=None, checkpointOutFile="gaps_checkpoint.out",
           checkpointInterval=0, checkpointInFile="", transposeData=None,
           BPPARAM=None, workerID=1, asynchronousUpdates=None, nSnapshots=0,
           snapshotPhase='sampling'):
    """
    Python wrapper to run CoGAPS via bindings
    @param path: path to data
    @param params: GapsParameters object
    @param nThreads: number of threads to use
    @param messages: whether to print messages
    @param outputFrequency:
    @param uncertainty:
    @param checkpointOutFile: path to where checkpoint info should be written
    @param checkpointInterval: how often to make a checkpoint
    @param checkpointInFile:
    @param transposeData:
    @param BPPARAM:
    @param workerID:
    @param asynchronousUpdates:
    @param nSnapshots:
    @param snapshotPhase: one of "sampling", "equilibration", "all"
    @return: a CogapsResult object
    """
    

    # check OpenMP support
    if isCompiledWithOpenMPSupport() is False:
        if asynchronousUpdates is not None and nThreads > 1:
            print("requesting multi-threaded version of CoGAPS but compiler did not support OpenMP")
        asynchronousUpdates = False
        nThreads = 1
    # convert sampling phase to enum
    if snapshotPhase == "sampling":
        snapshotPhase = pycogaps.GAPS_SAMPLING_PHASE
    elif snapshotPhase == "equilibration":
        snapshotPhase = pycogaps.GAPS_EQUILIBRATION_PHASE
    elif snapshotPhase == "all":
        snapshotPhase = pycogaps.GAPS_ALL_PHASES
    else:
        print("The snapshot phase you indicated is not recognized.")
        print("Please choose one of: sampling, equilibration, all")
        return

    if params is None:
        # construct a parameters object using whatever was passed in
        prm = GapsParameters(path)
        opts = {
            'maxThreads': nThreads,
            'printMessages': messages,
            'outputFrequency': outputFrequency,
            'checkpointOutFile': checkpointOutFile,
            'checkpointInterval': checkpointInterval,
            'checkpointFile': checkpointInFile,
            'transposeData': transposeData,
            'workerID': workerID,
            'asynchronousUpdates': asynchronousUpdates,
            'snapshotFrequency': nSnapshots,
            'snapshotPhase': snapshotPhase,
        }
        setParams(prm, opts)
        # check data input
        helper_functions.checkData(path, prm)
        return pycogaps.runCogaps(path, prm)
    else:
        # should we allow them to pass in params?
        # it's hard because we can't distinguish
        # between defaults and user-supplied params AFAIK
        # check data input
        helper_functions.checkData(path, params)
        return pycogaps.runCogaps(path, params)


# TODO: should we pass uncertainty into runCogaps?


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


def setParams(paramobj, list):
    """

    @param paramobj: a GapsParameters object
    @param list: key:value pairings for each parameter you wish to set
    """
    for (k, v) in list.items():
        setParam(paramobj, k, v)


# class CogapsParams
# constructor has default values for each parameter
def setParam(paramobj, whichParam, value):
    """

    @param paramobj: a GapsParameters object
    @param whichParam: the name of the parameter you wish to change
    @param value: the value to set whichParam as
    @return: nothing; paramobj will be modified
    """
    if whichParam == "alpha":
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
    elif whichParam in "singleCell":
        print(whichParam, " has been deprecated, this parameter will be ignored")
        return
    else:
        setattr(paramobj, whichParam, value)
    return paramobj


def getParam(paramobj, whichParam):
    return getattr(paramobj, whichParam)


def show(paramobj):
    paramobj.print()
