import time
import math
from numpy.core.fromnumeric import transpose

import pandas as pd
import pycogaps
import anndata
from helper_functions import *

class CoParams:
    '''
    self.gaps : GapsParameters object
    self.cogaps : dictionary of additional parameters (not in GapsParameters)
    '''
    def __init__(self, path=None, matrix=None, transposeData=False, hdfKey=None):
        if matrix is not None:
            self.gaps = GapsParameters(matrix)
        elif path is not None:
            if not path.lower().endswith(".h5"):
                adata = toAnndata(path, transposeData=transposeData)
            else:
                adata = toAnndata(path, hdfKey, transposeData=transposeData)
            matrix = pycogaps.Matrix(adata.X)
            self.gaps = GapsParameters(matrix)
        else:
            raise Exception('initialize with path= or params=')

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
                            'distributed': None,
                            'hdfKey': hdfKey,
                            'useSparseOptimization': None,
                            'transposeData': transposeData,
                        }
        self.coparams['minNS'] = math.ceil(self.coparams['cut'] / 2)
        self.coparams['maxNS'] = self.coparams['minNS'] + self.coparams['nSets']

    def setDistributedParams(self, nSets=None, cut=None, minNS=None, maxNS=None):
        print("setting distributed parameters - call this again if you change nPatterns")
        if self.coparams['distributed'] != "genome-wide":
            print("if you wish to perform genome-wide distributed cogaps, please run setParams(params, "
                  "\"distributed\", ""\"genome-wide\")")
        if nSets is None:
            self.coparams['nSets'] = self.coparams['nSets']
        else:
            self.coparams['nSets'] = nSets
        if cut is None:
            self.coparams['cut'] = self.gaps.nPatterns
        else:
            self.coparams['cut'] = cut
        if minNS is None:
            self.coparams['minNS'] = math.ceil(self.coparams['cut'] / 2)
        else:
            self.coparams['minNS'] = minNS
        if maxNS is None:
            self.coparams['maxNS'] = self.coparams['minNS'] + self.coparams['nSets']
        else:
            self.coparams['maxNS'] = minNS

    # samplingWeight is a dictionary
    # can use: dict(zip(names, weights))
    def setAnnotationWeights(self, annotation, weight):
        self.coparams['samplingAnnotation'] = annotation
        self.coparams['samplingWeight'] = weight

    def setFixedPatterns(self, fixedPatterns, whichMatrixFixed):
        self.coparams['fixedPatterns'] = fixedPatterns
        self.coparams['whichMatrixFixed'] = whichMatrixFixed
        self.gaps.useFixedPatterns = True
        self.gaps.fixedPatterns = pycogaps.Matrix(fixedPatterns)
        self.gaps.whichMatrixFixed = whichMatrixFixed

    # print standard and sparsity parameters
    def printParams(self):
        print('\n-- Standard Parameters --')
        print('nPatterns: ', self.gaps.nPatterns)
        print('nIterations: ', self.gaps.nIterations)
        print('seed: ', self.gaps.seed)
        print('sparseOptimization: ', self.gaps.useSparseOptimization)
        print('\n')
        print('-- Sparsity Parameters --')
        print('alpha: {:0.2f}'.format(self.gaps.alphaA))
        print('maxGibbsMass: ', self.gaps.maxGibbsMassA)
        print('\n')
        if self.gaps.runningDistributed:
            print('-- Distributed Parameters --')
            print('cut: ', self.coparams['cut'])
            print('nSets: ', self.coparams['nSets'])
            print('minNS: ', self.coparams['minNS'])
            print('maxNS: ', self.coparams['maxNS'])
            print('\n')
        
    # print all GapsParameters 
    def printAllParams(self):
        self.gaps.print()


def CoGAPS(path, params=None, nThreads=1, messages=True,
           outputFrequency=1000, uncertainty=None, checkpointOutFile="gaps_checkpoint.out",
           checkpointInterval=0, checkpointInFile="", transposeData=False,
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

    # convert data to anndata and matrix obj
    if params is not None:
        adata = toAnndata(path, params.coparams['hdfKey'], transposeData=transposeData)
    else:
        adata = toAnndata(path, transposeData=transposeData)

    matrix = pycogaps.Matrix(adata.X)

    if params is None:
        prm = CoParams(matrix=matrix, transposeData=transposeData)
    else:
        prm = params

    opts = {
        'maxThreads': nThreads,
        'printMessages': messages,
        'outputFrequency': outputFrequency,
        'checkpointOutFile': checkpointOutFile,
        'checkpointInterval': checkpointInterval,
        'checkpointFile': checkpointInFile,
        # 'transposeData': transposeData,
        'workerID': workerID,
        'asynchronousUpdates': asynchronousUpdates,
        'snapshotFrequency': nSnapshots,
        'snapshotPhase': snapshotPhase,
    }
    setParams(prm, opts)

    if uncertainty is not None:
        unc = toAnndata(uncertainty)
        unc = pycogaps.Matrix(unc.X)
    else:
        unc = pycogaps.Matrix()


    prm = getDimNames(adata, prm)

    # check data input
    checkData(adata, prm.gaps, uncertainty) 
    checkInputs(uncertainty, prm)

    startupMessage(prm, path)
    gapsresultobj = pycogaps.runCogapsFromMatrix(matrix, prm.gaps, unc)
    prm.gaps.transposeData = transposeData

    if prm.gaps.transposeData != prm.coparams["transposeData"]:
        raise Exception("make sure to pass transposeData=True argument in both CoParams() and CoGAPS()")

    result = {
        "GapsResult": gapsresultobj,
        "anndata": GapsResultToAnnData(gapsresultobj, adata, prm.gaps)
    }

    show(result["anndata"])
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
    print(adata.X)
    print(Amean)
    print(Pmean)
    adata.obs = pd.DataFrame(data=Amean, index=adata.obs_names, columns=pattern_labels)
    adata.var = pd.DataFrame(data=Pmean, index=adata.var_names, columns=pattern_labels)
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
    @return: nothing paramobj will be modified
    """
    coparam_params = ['hdfKey', 'explicitSets', 'subsetDim', 'geneNames', 'sampleNames']
    if whichParam == "alpha":
        paramobj.gaps.alphaA = value
        paramobj.gaps.alphaP = value
    elif whichParam == "maxGibbsMass":
        paramobj.gaps.maxGibbsMassA = value
        paramobj.gaps.maxGibbsMassP = value
    elif whichParam in coparam_params:
        paramobj.coparams[whichParam] = value
    elif whichParam in ("distributed"):
        if value == "genome-wide":
            paramobj.gaps.runningDistributed = True
            paramobj.coparams['distributed'] = value
        elif (value is not None) and (value is not False):
            print("if you wish to perform genome-wide distributed cogaps, please run setParams(params, "
                  "\"distributed\", ""\"genome-wide\")")
    elif whichParam in ("nSets", "cut", "minNS", "maxNS"):
        paramobj.gaps.runningDistributed = True
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
        setattr(paramobj.gaps, whichParam, value)
    return paramobj


def getParam(paramobj, whichParam):
    return getattr(paramobj, whichParam)

