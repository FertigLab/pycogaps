import pandas as pd

from PyCoGAPS.config import *

import h5py
import scipy.io
import pkg_resources  # part of setuptools
from pycogaps import getElement

def supported(file):
    """ Checks whether file is supported type

    Args:
        file (str): path to data

    Returns:
        bool: true/false if file is supported
    """    
    return file.lower().endswith((".tsv", ".csv", ".mtx", ".h5ad", ".h5", ".gct", ".txt"))


def checkData(adata, params, uncertainty=None):
    """ Check validity of inputted data

    Args:
        adata (anndata): data as anndata object
        params (CoParams): CoParams object
        uncertainty (arr, optional): optional uncertainty matrix. Defaults to None.

    Raises:
        Exception: If NA values are present in data
        Exception: If data is not numeric
        Exception: If negative values are in data matrix
        Exception: If negative values in uncertainty matrix
        Warning: If small values in uncertainty matrix
        Exception: If nPatterns is greater than dimensions of data
    """    
    data = adata.X

    if np.isnan(data).any():
        raise Exception('NA values in data')
    if not np.issubdtype(data.dtype, np.number):
        raise Exception('data is not numeric')
    if np.any((data < 0)):
        raise Exception('negative values in data matrix')
    if uncertainty != None:
        if np.any((uncertainty < 0)):
            raise Exception('negative values in uncertainty matrix')
        if np.any(uncertainty < 1e-5):
            raise Warning('small values in uncertainty matrix detected')
    if data.shape[0] <= params.nPatterns | data.shape[1] <= params.nPatterns:
        raise Exception('nPatterns must be less than dimensions of data')


def toAnndata(file, hdf_counts_key=None, hdf_dim1_key=None, hdf_dim2_key=None, transposeData=False):
    """ Converts file to anndata object

    Args:
        file (str): path to data
        hdf_counts_key (str, optional): if .h5 data inputted. Defaults to None.
        hdf_dim1_key (str, optional): if .h5 data inputted. Defaults to None.
        hdf_dim2_key (str, optional): if .h5 data inputted. Defaults to None.
        transposeData (bool, optional): if data should be transposed. Defaults to False.

    Raises:
        Exception: if unsupported data type
        Exception: if dataset name from hdf file is not provided to CoParams

    Returns:
        anndata: anndata object
    """   

    '''
    TODO: CHANGE TO READ AS NUMPY ARRAY FIRST, THEN MAKE DATAFRAME, THEN MAKE ANNDATA
    '''
    if not supported(file):
        raise Exception("unsupported data type")
    if file.lower().endswith(".csv"):
        adata = anndata.read_csv(file)
    elif file.lower().endswith(".txt"):
        # table = pd.read_table(file)
        # adata = anndata.AnnData(table.iloc[:, 2:])
        # adata.obs_names = table["symbol"]
        pd_table = pd.read_table(file)
        table = pd.DataFrame(data=pd_table.values, index=pd_table.index, columns=pd_table.columns)
        adata = anndata.AnnData(table)
    elif file.lower().endswith(".tsv"):
        csv_table = pd.read_table(file,sep='\t')
        csv_table.to_csv('file.csv', index=False)
        adata = anndata.read_csv('{}.csv'.format(os.path.splitext(file)[0]))
    elif file.lower().endswith(".mtx"):
        adata = anndata.read_mtx(file)
    elif file.lower().endswith(".h5ad"):
        adata = anndata.read_h5ad(file)
    elif file.lower().endswith(".h5"):
        if hdf_counts_key is None:
            raise Exception("set dataset name from hdf file to use with params = CoParams(path=filename, hdfKey=key")
        adata = anndata.read_hdf(file, hdf_counts_key) # user supplied keydata
        if transposeData:
            if hdf_dim1_key is not None:
                adata.obs_names = h5py.File(file, 'r')[hdf_dim1_key]
            if hdf_dim2_key is not None:
                adata.var_names = h5py.File(file, 'r')[hdf_dim2_key]
        else:
            if hdf_dim1_key is not None:
                adata.var_names = h5py.File(file, 'r')[hdf_dim1_key]
            if hdf_dim2_key is not None:
                adata.obs_names = h5py.File(file, 'r')[hdf_dim2_key]
    elif file.lower().endswith(".gct"):
        csv_table = pd.read_csv(file, sep='\t', skiprows=2)
        csv_table.to_csv('file.csv', index=False)
        adata = anndata.read_csv('{}.csv'.format(os.path.splitext(file)[0]))

    if scipy.sparse.issparse(adata.X):
        adata.X = (adata.X).toarray()
    
    if transposeData:
        adata = adata.transpose()

    return adata


def checkInputs(uncertainty, allParams):
    """ Check validity of inputs to CoGAPS.

    Args:
        uncertainty (arr): uncertainty matrix
        allParams (CoParams): CoParams object

    Raises:
        Exception: If unsupported file extension provided for uncertainty matrix.
        Exception: If default uncertainty not used with useSparseOptimization=True
        Exception: If CoGAPS was built with checkpoints disabled
        Exception: If checkpoints not supported for distributed CoGAPS
    """    
    if uncertainty is not None and not supported(uncertainty):
            raise Exception("unsupported file extension for uncertainty")
    if uncertainty is not None and allParams.coparams["useSparseOptimization"] is True:
        raise Exception("must use default uncertainty when enabling useSparseOptimization")
    if allParams.gaps.checkpointFile is not None and not isCheckpointsEnabled():
        raise Exception("CoGAPS was built with checkpoints disabled")
    if allParams.gaps.snapshotFrequency > 0:
        warnings.warn("snapshots slow down computatioin and shouldo nly be used for testing")

    if allParams.coparams["distributed"] is not None:
        if allParams.gaps.maxThreads > 1:
            warnings.warn("can't run multi-threaded and distributed CoGAPS at the same time, ignoring nThreads")
        if allParams.gaps.checkpointFile != "":
            raise Exception("checkpoints not supported for distributed CoGAPS")


def nrowHelper(data):
    return data.shape[0]


def ncolHelper(data):
    return data.shape[1]


def getGeneNames(data, transpose):
    """ Return gene names

    Args:
        data (anndata): data as matrix
        transpose (bool): if data was transposed

    Returns:
        list: list of names
    """    
    if transpose:
        return getSampleNames(data, False)
    names = data.obs_names

    if names.all() == None or len(names) == 0:
        return ["Gene" + str(i) for i in range(1, nrowHelper(data))]
    return names


def getSampleNames(data, transpose):
    """ Return sample names

    Args:
        data (anndata): data as matrix
        transpose (bool): if data was transposed

    Returns:
        list: list of names
    """  
    if transpose:
        return getGeneNames(data, False)
    names = data.var_names

    if names.all() == None or len(names) == 0:
        return ["Sample" + str(i) for i in range(1, ncolHelper(data))]
    return names



def getDimNames(data, allParams):
    # support both path and anndata object as data input
    """ Get dimension names

    Args:
        data (arr): data as matrix
        allParams (CoParams): CoParams object

    Returns:
        CoParams: updated CoParams object
    """  
    if isinstance(data, str):
        if not supported(data):
            raise Exception("unsupported data type")
        else:
            data = toAnndata(data).X

    geneNames = getGeneNames(data, allParams.gaps.transposeData)
    sampleNames = getSampleNames(data, allParams.gaps.transposeData)

    if allParams.gaps.transposeData:
        nGenes = ncolHelper(data)
        nSamples = nrowHelper(data)
    else:
        nGenes = nrowHelper(data)
        nSamples = ncolHelper(data)

    if allParams.coparams['subsetDim'] == 1:
        nGenes = len(allParams.coparams['subsetIndices'])
        geneNames = np.take(geneNames, allParams.coparams['subsetIndices'])
    elif allParams.coparams['subsetDim'] == 2:
        nSamples = len(allParams.coparams['subsetIndices'])
        sampleNames = np.take(sampleNames, allParams.coparams['subsetIndices'])

    if len(geneNames) != nGenes:
        raise Exception(len(geneNames), " != ", nGenes, " incorrect number of gene names given")
    if len(sampleNames) != nSamples:
        raise Exception(len(sampleNames), " != ", nSamples, " incorrect number of sample names given")

    # store processed gene/sample names directly in allParams list
    # this is an important distinction - allParams@gaps contains the
    # gene/sample names originally passed by the user, allParams contains
    # the procseed gene/sample names to be used when labeling the result

    allParams.coparams['geneNames'] = geneNames
    allParams.coparams['sampleNames'] = sampleNames
    return (allParams)

def startupMessage(params, path):
    """ Message to display at startup

    Args:
        params (CoParams): CoParams object
        path (str): path to data
    """    
    print("\nThis is ", end='')
    getVersion()

    dist_message = "Standard"
    if params.coparams["distributed"] is not None and params.coparams["distributed"] is not False:
        dist_message = params.coparams["distributed"]
    if isinstance(path, str):
        data_name = os.path.basename(path)
    else:
        data_name = "provided data object"

    print("Running", dist_message, "CoGAPS on", data_name, "(", len(params.coparams['geneNames']), "genes and", len(params.coparams['sampleNames']),"samples)",
    "with parameters: ")
    params.printParams()


def show(obj: anndata):
    """ Concluding message after CoGAPS completes run

    Args:
        obj (anndata): anndata object
    """    
    nfeatures = obj.n_obs
    nsamples = obj.n_vars
    npatterns = len(obj.obs_keys())
    print("\nGapsResult result object with", nfeatures, "features and", nsamples, "samples")
    print(npatterns, "patterns were learned\n")
    return



def getFeatureLoadings(object: anndata):
    """ Get feature loadings matrix

    Args:
        object (anndata): anndata object

    Returns:
        arr: array of matrix
    """    
    return object.obs


def getAmplitudeMatrix(object):
    """ Get amplitude matrix

    Args:
        object (anndata): anndata object

    Returns:
        arr: array of matrix
    """  
    return object.obs


def getSampleFactors(object):
    """ Get sample factors matrix

    Args:
        object (anndata): anndata object

    Returns:
        arr: array of matrix
    """  
    return object.var


def getPatternMatrix(object):
    """ Get pattern matrix

    Args:
        object (anndata): anndata object

    Returns:
        arr: array of matrix
    """  
    return object.var

def getMeanChiSq(object):
    """ get mean chi sq

    Args:
        object (CogapsResult): CogapsResult object

    Returns:
        [type]: mean chi sq value
    """    
    object = object["GapsResult"]
    return object.meanChiSq


def getVersion():
    """ Prints version of PyCoGAPS package

    Returns:
        str: version number
    """    
    version = pkg_resources.require("pycogaps")[0].version
    print("pycogaps version ", version)
    return version

'''
# TODO need to access original params through anndata
def getOriginalParameters(object: GapsResult):
    print("Not yet implemented")
    return


def getUnmatchedPatterns(object):
    print("Not yet implemented")
    return


def getClusteredPatterns(object):
    print("Not yet implemented")
    return


def getCorrelationToMeanPattern(object):
    print("Not yet implemented")
    return


def getSubsets(object):
    print("Not yet implemented")
    return
'''


# convert matrix object to numpy array
def toNumpy(matrix):
    """ Convert matrix object to numpy array

    Args:
        matrix (Matrix): a Matrix object

    Returns:
        arr: a numpy array
    """    
    all_vecdata = np.empty((matrix.nRow(), matrix.nCol()))
    for i in range(matrix.nCol()):
        vector = matrix.getCol(i)
        vecdata = []
        for j in range(vector.size()):
            vecdata.append(getElement(vector, j))
        all_vecdata[:, i] = vecdata
    return all_vecdata

def toNumpyFromVector(vector):
    """ Convert vector object to numpy array

    Args:
        vector (Matrix): a vector<Matrix> object

    Returns:
        arr: a numpy array
    """    
    arr = np.empty(len(vector))
    for j in range(len(vector)):
        matrix = getElement(vector, j)
        arr[j] = toNumpy(matrix)
    return arr

def GapsResultToAnnData(gapsresult, adata, prm):
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
    # if prm.coparams["distributed"] == "genome-wide":
    adata.obs = pd.DataFrame(data=Amean, index=adata.obs_names, columns=pattern_labels)
    adata.var = pd.DataFrame(data=Pmean, index=adata.var_names, columns=pattern_labels)
    adata.uns["asd"] = pd.DataFrame(data=Asd, index=adata.obs_names, columns=pattern_labels)
    adata.uns["psd"] = pd.DataFrame(data=Psd, index=adata.var_names, columns=pattern_labels)
    # else:
    #     adata.obs = pd.DataFrame(data=Pmean, index=adata.obs_names, columns=pattern_labels)
    #     adata.var = pd.DataFrame(data=Amean, index=adata.var_names, columns=pattern_labels)
    #     adata.uns["asd"] = pd.DataFrame(data=Asd, index=adata.var_names, columns=pattern_labels)
    #     adata.uns["psd"] = pd.DataFrame(data=Psd, index=adata.obs_names, columns=pattern_labels)
    adata.uns["atomhistoryA"] = pd.Series(gapsresult.atomHistoryA)
    adata.uns["atomhistoryP"] = pd.Series(gapsresult.atomHistoryP)
    adata.uns["averageQueueLengthA"] = float(gapsresult.averageQueueLengthA)
    adata.uns["averageQueueLengthP"] = float(gapsresult.averageQueueLengthP)
    adata.uns["chisqHistory"] = pd.Series(gapsresult.chisqHistory)
    adata.uns["equilibrationSnapshotsA"] = toNumpyFromVector(gapsresult.equilibrationSnapshotsA)
    adata.uns["equilibrationSnapshotsP"] = toNumpyFromVector(gapsresult.equilibrationSnapshotsP)
    adata.uns["meanChiSq"] = float(gapsresult.meanChiSq)
    adata.uns["meanPatternAssignment"] = toNumpy(gapsresult.meanPatternAssignment)
    adata.uns["pumpMatrix"] = toNumpy(gapsresult.pumpMatrix)
    adata.uns["samplingSnapshotsA"] = toNumpyFromVector(gapsresult.samplingSnapshotsA)
    adata.uns["samplingSnapshotsP"] = toNumpyFromVector(gapsresult.samplingSnapshotsP)
    adata.uns["seed"] = int(gapsresult.seed)
    adata.uns["totalRunningTime"] = int(gapsresult.totalRunningTime)
    adata.uns["totalUpdates"] = int(gapsresult.totalUpdates)

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



