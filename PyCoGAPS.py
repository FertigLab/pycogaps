import time
import math

import pandas as pd
import pycogaps
import anndata
from helper_functions import *

class CoParams:
    '''
    self.gaps : GapsParameters object
    self.cogaps : dictionary of additional parameters (not in GapsParameters)
    '''
    def __init__(self, path=None, matrix=None, params=None):
        if matrix is not None:
            self.gaps = GapsParameters(matrix)
        if path is not None:
            adata = toAnndata(path)
            matrix = pycogaps.Matrix(adata.X)
            self.gaps = GapsParameters(path)
        if params is not None:
            self.gaps = params

        self.coparams = {'cut': self.gaps.nPatterns,
                            'nSets': 4,
                            'minNS': None,
                            'maxNS': None,
                            'explicitSets': None,
                            'samplingAnnotation': None,
                            'samplingWeight': None,
                            'subsetIndices': None,
                            'subsetDim': 0,
                            'geneNames': None,
                            'sampleNames': None,
                            'fixedPatterns': None,
                        }
        self.coparams['minNS'] = math.ceil(self.coparams['cut'] / 2)
        self.coparams['maxNS'] = self.coparams['minNS'] + self.coparams['nSets']


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

    # convert data to anndata and matrix obj
    adata = toAnndata(path)
    matrix = pycogaps.Matrix(adata.X)

    gapsresultobj = None
    if params is None:
        # construct a parameters object using whatever was passed in
        prm = CoParams(matrix=matrix)
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
    else:
        # params passed in should probably be type CoParams, but just in case
        if isinstance(params, CoParams):
            prm = params
        else:
            prm = CoParams(params=params)

    # check data input
    checkData(adata, prm.gaps, uncertainty)  
    gapsresultobj = pycogaps.runCogapsFromMatrix(matrix, prm.gaps)

    result = {
        "GapsResult": gapsresultobj,
        "anndata": GapsResultToAnnData(gapsresultobj, adata, prm.gaps)
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


def setParams(paramobj: CoParams, list):
    """

    @param paramobj: a GapsParameters object
    @param list: key:value pairings for each parameter you wish to set
    """
    for (k, v) in list.items():
        setParam(paramobj, k, v)


# class CogapsParams
# constructor has default values for each parameter
def setParam(paramobj: CoParams, whichParam, value):
    """

    @param paramobj: a GapsParameters object
    @param whichParam: the name of the parameter you wish to change
    @param value: the value to set whichParam as
    @return: nothing; paramobj will be modified
    """

    coparams = paramobj.coparams

    if whichParam == "alpha":
        paramobj.gaps.alphaA = value
        paramobj.gaps.alphaP = value
    elif whichParam == "maxGibbsMass":
        paramobj.gaps.maxGibbsMassA = value
        paramobj.gaps.maxGibbsMassP = value
    elif whichParam in ("nSets", "cut", "minNS", "maxNS"):
        print("setting distributed parameters - call this again if you change nPatterns")
        coparams[whichParam] = value
        ### TODO: update the rest
    elif whichParam in ("samplingAnnotation", "samplingWeight"):
        coparams[whichParam] = value
    elif whichParam in ("fixedPatterns", "whichMatrixFixed"):
        coparams[whichParam] = value

    elif whichParam in "singleCell":
        print(whichParam, " has been deprecated, this parameter will be ignored")
        return
    else:
        setattr(paramobj, whichParam, value)
    return paramobj


def getParam(paramobj, whichParam):
    return getattr(paramobj, whichParam)

