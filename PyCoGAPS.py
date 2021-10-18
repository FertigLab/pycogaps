import time
import math

import pandas as pd
import pycogaps
import anndata
from helper_functions import *
from subset_data import *
# from distributed import *
import sys

sys.setrecursionlimit(10000)

print("""\
______      _____       _____   ___  ______  _____ 
| ___ \    /  __ \     |  __ \ / _ \ | ___ \/  ___|
| |_/ /   _| /  \/ ___ | |  \// /_\ \| |_/ /\ `--. 
|  __/ | | | |    / _ \| | __ |  _  ||  __/  `--. |
| |  | |_| | \__/\ (_) | |_\ \| | | || |    /\__/ /
\_|   \__, |\____/\___/ \____/\_| |_/\_|    \____/ 
       __/ |                                       
      |___/                                        
                    """)


class CoParams:
    """ Encapsulates all parameters for PyCoGAPS.

    """

    def __init__(self, path=None, matrix=None, transposeData=False, hdfKey=None, hdfRowKey=None, hdfColKey=None):
        """ Initializes CoParams object. 
            self.gaps : GapsParameters object
            self.cogaps : dictionary of additional parameters (not in GapsParameters)

        Args:
            path (str, optional): Path to data. Defaults to None.
            matrix (anndata, optional): AnnData object containing supplied data matrix. Defaults to None.
            transposeData (bool, optional): Expects genes x samples. Defaults to False.
            hdfKey (str, optional): For reading .h5 files. Defaults to None.
            hdfRowKey (str, optional): For reading .h5 files. Defaults to None.
            hdfColKey (str, optional): For reading .h5 files. Defaults to None.

        Raises:
            Exception: If path or params not passed as an argument.
        """        

        if matrix is not None:
            self.gaps = GapsParameters(pycogaps.Matrix(matrix.X))
            adata = matrix
        elif path is not None:
            if path.lower().endswith(".h5"):
                adata = toAnndata(path, hdfKey, hdfRowKey, hdfColKey, transposeData=transposeData)
            elif path.lower().endswith(".txt"):
                table = pd.read_table(path)
                adata = anndata.AnnData(table.iloc[:, 2:])
                adata.obs_names = table["symbol"]
            else:
                adata = toAnndata(path, transposeData=transposeData)
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
                         'geneNames': adata.obs_names,
                         'sampleNames': adata.var_names,
                         'fixedPatterns': None,
                         'distributed': "genome-wide",
                         'hdfKey': hdfKey,
                         'hdfRowKey': hdfRowKey,
                         'hdfColKey': hdfColKey,
                         'useSparseOptimization': None,
                         'transposeData': transposeData,
                         }
        self.coparams['minNS'] = math.ceil(self.coparams['cut'] / 2)
        self.coparams['maxNS'] = self.coparams['minNS'] + self.coparams['nSets']

    def setDistributedParams(self, nSets=None, cut=None, minNS=None, maxNS=None):
        """ Sets parameters for running distributed CoGAPS.

        Args:
            nSets (int, optional): Number of sets to break data into. Defaults to None.
            cut (int, optional): Number of branches at which to cut dendrogram used in pattern matching. Defaults to None.
            minNS (int, optional): [description]. Minimum of individual set contributions a cluster must contain. Defaults to None.
            maxNS (int, optional): [description]. Maximum of individual set contributions a cluster can contain. Defaults to None.
        """        

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
            self.coparams['maxNS'] = maxNS

    # samplingWeight is a dictionary
    # can use: dict(zip(names, weights))
    def setAnnotationWeights(self, annotation, weight):
        """ Set annotation weights for distributed CoGAPS.

        Args:
            annotation (str list): Specify categories along the rows (cols) to use for weighted sampling.
            weight (int list): Weights associated with samplingAnnotation
        """    

        self.coparams['samplingAnnotation'] = annotation
        self.coparams['samplingWeight'] = weight

    def setFixedPatterns(self, fixedPatterns, whichMatrixFixed):
        """ Fix either 'A' or 'P' matrix to given values.

        Args:
            fixedPatterns (arr): Fix either 'A' or 'P' matrix to these values, 
            in the context of distributed CoGAPS, the first phase is skipped and 
            fixedPatterns is used for all sets - allowing manual pattern matching, 
            as well as fixed runs of standard CoGAPS.
            whichMatrixFixed (str): Either 'A' or 'P', indicating which matrix is fixed
        """ 

        self.coparams['fixedPatterns'] = fixedPatterns
        self.coparams['whichMatrixFixed'] = whichMatrixFixed
        self.gaps.useFixedPatterns = True
        self.gaps.fixedPatterns = pycogaps.Matrix(fixedPatterns)
        self.gaps.whichMatrixFixed = whichMatrixFixed

    def printParams(self):
        """ Print standard and sparsity parameters, and distributed if set.
        """        
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

    def printAllParams(self):
        """ Print all GapsParameters.
        """        
        self.gaps.print()


def CoGAPS(path, params=None, nThreads=1, messages=True,
           outputFrequency=1000, uncertainty=None, checkpointOutFile="",
           checkpointInterval=0, checkpointInFile="", transposeData=False,
           BPPARAM=None, workerID=1, asynchronousUpdates=None, nSnapshots=0,
           snapshotPhase='sampling'):
    """ Python wrapper to run CoGAPS via bindings

    Args:
        path (str): Path to data. 
        params (CoParams, optional): CoParams object of parameters. Defaults to None.
        nThreads (int, optional): Number of threads to use. Defaults to 1.
        messages (bool, optional): Whether to print messages. Defaults to True.
        outputFrequency (int, optional): How often to output messages. Defaults to 1000.
        uncertainty (arr, optional): Optional uncertainty matrix. Defaults to None.
        checkpointOutFile (str, optional): Path to where checkpoint info should be written. Defaults to "".
        checkpointInterval (int, optional): How often to make a checkpoint. Defaults to 0.
        checkpointInFile (str, optional): Path to existing checkpoint file to run CoGAPS from. Defaults to "".
        transposeData (bool, optional): Expects genes x samples. Defaults to False.
        BPPARAM ([type], optional): BiocParallel backend . Defaults to None.
        workerID (int, optional): If calling CoGAPS in parallel the worker ID can be specified,
        only worker 1 prints output and each worker outputs when it finishes, this
        is not neccesary when using the default parallel methods (i.e. distributed
        CoGAPS) but only when the user is manually calling CoGAPS in parallel. Defaults to 1.
        asynchronousUpdates (bool, optional): Enable asynchronous updating which allows for multi-threaded runs. Defaults to None.
        nSnapshots (int, optional): How many snapshots to take in each phase, setting this to 0 disables snapshots. Defaults to 0.
        snapshotPhase (str, optional): One of "sampling", "equilibration", "all". Defaults to 'sampling'.

    Raises:
        Exception: If transposeData=True is not passed as an argument to both CoParams and CoGAPS.

    Returns:
        CogapsResult: A CogapsResult object.
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
    if isinstance(path, str):
        if params is not None:
            adata = toAnndata(path, params.coparams['hdfKey'], params.coparams['hdfRowKey'],
                              params.coparams['hdfColKey'], transposeData=transposeData)
        else:
            adata = toAnndata(path, transposeData=transposeData)
    else:
        adata = path

    matrix = pycogaps.Matrix(adata.X)

    if params is None:
        prm = CoParams(matrix=adata, transposeData=transposeData)
    else:
        prm = params

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

    if uncertainty is not None:
        unc = toAnndata(uncertainty)
        unc = pycogaps.Matrix(unc.X)
    else:
        unc = pycogaps.Matrix()

    if prm.coparams["subsetIndices"] is None:
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
        "anndata": GapsResultToAnnData(gapsresultobj, adata, prm)
    }

    show(result["anndata"])
    return result


def GapsResultToAnnData(gapsresult, adata, prm: CoParams):
    """ Converts a CogapsResult object to anndata object.

    Args:
        gapsresult (CogapsResult): Dictionary result object.
        adata (anndata): Anndata object populated by CoGAPS.
        prm (CoParams): CoParams object.

    Returns:
        anndata: An anndata object.
    """    
    # need to subset matrices based on which dimension we're in...
    if prm.coparams['subsetDim'] == 1:
        Amean = toNumpy(gapsresult.Amean)[prm.coparams["subsetIndices"], :]
        Pmean = toNumpy(gapsresult.Pmean)
        Asd = toNumpy(gapsresult.Asd)[prm.coparams["subsetIndices"], :]
        Psd = toNumpy(gapsresult.Psd)
    else:
        Amean = toNumpy(gapsresult.Amean)
        Pmean = toNumpy(gapsresult.Pmean)[prm.coparams["subsetIndices"], :]
        Asd = toNumpy(gapsresult.Asd)
        Psd = toNumpy(gapsresult.Psd)[prm.coparams["subsetIndices"], :]
    pattern_labels = ["Pattern" + str(i) for i in range(1, prm.gaps.nPatterns + 1)]
    # load adata obs and var from Amean and Pmean results
    if len(Pmean.shape) > 2:
        Pmean = Pmean[0, :, :]
        Psd = Psd[0, :, :]

    adata.obs = pd.DataFrame(data=Amean, index=adata.obs_names, columns=pattern_labels)
    adata.var = pd.DataFrame(data=Pmean, index=adata.var_names, columns=pattern_labels)
    adata.uns["asd"] = pd.DataFrame(data=Asd, index=adata.obs_names, columns=pattern_labels)
    adata.uns["psd"] = pd.DataFrame(data=Psd, index=adata.var_names, columns=pattern_labels)
    return adata


def GapsParameters(path):
    """ Returns C++ GapsParameters object.

    Args:
        path (str): Path to data.

    Returns:
        GapsParameters: A GapsParameters object.
    """    
    return pycogaps.GapsParameters(path)


def getBuildReport():
    """ Returns information about how the package was compiled, i.e. which
    compiler/version was used, which compile time options were enabled, etc...

    Returns:
        str: String containing build report.
    """    
    return pycogaps.getBuildReport()


def isCheckpointsEnabled():
    """ Check if package was built with checkpoints enabled

    Returns:
        bool: true/false if checkpoints are enabled
    """    
    return pycogaps.isCheckpointsEnabled()


def isCompiledWithOpenMPSupport():
    """ Check if compiler supported OpenMP

    Returns:
        bool: true/false if OpenMP was supported
    """    
    return pycogaps.isCompiledWithOpenMPSupport()


def getFileInfo(path):
    """ Get info of inputted file.

    Args:
        path (str): Path to data.

    Returns:
        str: string of file info.
    """    
    return pycogaps.getFileInfo(path)


def current_milli_time():
    """ Return current time in milliseconds.

    Returns:
        int: Current time in milliseconds.
    """    
    return round(time.time() * 1000)


def setParams(paramobj: CoParams, list):
    """ Set CoParams from a list.

    Args:
        paramobj (CoParams): CoParams object.
        list (dict): Dictionary of parameter, value pairings for each parameter you wish to set.
    """    

    for (k, v) in list.items():
        setParam(paramobj, k, v)


# class CogapsParams
# constructor has default values for each parameter
def setParam(paramobj: CoParams, whichParam, value):
    """ Sets CoParams parameters.

    Args:
        paramobj (CoParams): a CoParams object
        whichParam ([type]): the name of the parameter you wish to change
        value ([type]): the value to set whichParam as

    Returns:
        CoParams: the modified CoParams object.
    """    

    coparam_params = ['hdfKey', 'hdfRowKey', 'hdfColKey', 'explicitSets', 'subsetDim', 'geneNames', 'sampleNames']
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
    """ Get parameter info and values.

    Args:
        paramobj (CoParams): a CoParams object.
        whichParam (str): which parameter to get the info of.

    Returns:
        [type]: the value of the parameter.
    """    
    return getattr(paramobj, whichParam)
