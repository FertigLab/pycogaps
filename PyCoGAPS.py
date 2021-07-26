import time
import math

import pandas as pd
import pycogaps
import anndata
from helper_functions import *


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

    gapsresultobj = None
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
        adata = toAnndata(path)
        checkData(adata, prm, uncertainty)
        matrix = pycogaps.Matrix(adata.X)
        gapsresultobj = pycogaps.runCogapsFromMatrix(matrix, prm)
    else:
        # should we allow them to pass in params?
        # it's hard because we can't distinguish
        # between defaults and user-supplied params AFAIK
        # check data input
        adata = toAnndata(path)
        checkData(adata, params, uncertainty)
        matrix = pycogaps.Matrix(adata.X)
        gapsresultobj = pycogaps.runCogapsFromMatrix(matrix, params)
        prm = params
    result = {
        "GapsResult": gapsresultobj,
        "anndata": GapsResultToAnnData(gapsresultobj, adata, prm)
    }
    return result


# TODO: should we pass uncertainty into runCogaps?


def GapsResultToAnnData (gapsresult:GapsResult, adata, prm:GapsParameters):
    # convert Amean and Pmean results to numpy arrays
    Amean = toNumpy(gapsresult.Amean)
    Pmean = toNumpy(gapsresult.Pmean)
    Asd = toNumpy(gapsresult.Asd)
    Psd = toNumpy(gapsresult.Psd)
    pattern_labels = ["Pattern" + str(i) for i in range(1, prm.nPatterns + 1)]
    # load adata obs and var from Amean and Pmean results
    A_mat = pd.DataFrame(data=Amean, index=adata.obs_names, columns=pattern_labels)
    adata.obs = A_mat
    P_mat = pd.DataFrame(data=Pmean, index=adata.var_names, columns=pattern_labels)
    adata.var = P_mat
    adata.uns["asd"] = pd.DataFrame(data=Asd, index=adata.obs_names, columns=pattern_labels)
    adata.uns["psd"] = pd.DataFrame(data=Psd, index=adata.var_names, columns=pattern_labels)
    return adata

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

