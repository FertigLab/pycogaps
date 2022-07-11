from PyCoGAPS.config import *
from PyCoGAPS.helper_functions import *
from PyCoGAPS.subset_data import *

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
                         'distributed': None,
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
        whichParam (str): the name of the parameter you wish to change
        value ([type]): the value to set whichParam as

    Returns:
        CoParams: the modified CoParams object.
    """    

    coparam_params = ['hdfKey', 'hdfRowKey', 'hdfColKey', 'explicitSets', 'subsetDim', 'geneNames', 'sampleNames', 'subsetIndices']
    if whichParam == "alpha":
        paramobj.gaps.alphaA = value
        paramobj.gaps.alphaP = value
    elif whichParam == "maxGibbsMass":
        paramobj.gaps.maxGibbsMassA = value
        paramobj.gaps.maxGibbsMassP = value
    elif whichParam in coparam_params:
        if value is not None:
            paramobj.coparams[whichParam] = value
    elif whichParam == "distributed":
        if value == "genome-wide":
            print("running genome-wide. if you wish to perform single-cell distributed cogaps, please run setParams(params, "
                  "\"distributed\", ""\"single-cell\")")
            paramobj.gaps.runningDistributed = True
            paramobj.coparams['distributed'] = value
        elif value == "single-cell":
            print("running single-cell. if you wish to perform genome-wide distributed cogaps, please run setParams(params, "
                  "\"distributed\", ""\"genome-wide\")")
            paramobj.coparams['distributed'] = 'single-cell'
            paramobj.gaps.runningDistributed = True
        else:
            paramobj.gaps.runningDistributed = False
            paramobj.coparams['distributed'] = None
    elif whichParam in ("nSets", "cut", "minNS", "maxNS"):
        paramobj.gaps.runningDistributed = True
        print("please set \'", whichParam, "\' with setDistributedParams")
        return
    elif whichParam in ("samplingAnnotation", "samplingWeight"):
        print("please set \'", whichParam, "\' with setAnnotationWeights")
        return
    elif whichParam in ("fixedPatterns", "whichMatrixFixed"):
        # print("please set \'", whichParam, "\' with setFixedPatterns")
        return
    elif whichParam == "singleCell":
        print(whichParam, " has been deprecated, this parameter will be ignored")
        return
    elif whichParam == "nThreads":
        whichParam = "maxThreads"
        setattr(paramobj.gaps, whichParam, value)
    elif whichParam == "messages":
        whichParam = "printMessages"
        setattr(paramobj.gaps, whichParam, value)
    elif whichParam == "nSnapshots":
        whichParam = "snapshotFrequency"
        setattr(paramobj.gaps, whichParam, value)
    elif whichParam == "checkpointInFile":
        whichParam = "checkpointFile"
        setattr(paramobj.gaps, whichParam, value)
    elif whichParam == "snapshotPhase":
        if value == "sampling":
            value = pycogaps.GAPS_SAMPLING_PHASE
        elif value == "equilibration":
            value = pycogaps.GAPS_EQUILIBRATION_PHASE
        elif value == "all":
            value = pycogaps.GAPS_ALL_PHASES
        # else:
        #     print("The snapshot phase you indicated is not recognized.")
        #     print("Please choose one of: sampling, equilibration, all")
        #     return
        setattr(paramobj.gaps, whichParam, value)
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