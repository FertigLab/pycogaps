import pycogaps
import numpy as np
import pandas as pd
import scipy.io

def supported(file):
    return file.lower().endswith((".tsv", ".csv", ".mtx", ".gct"))

# convert file types to pandas dataframe
def dataToDF(file):
    if file.lower().endswith(".csv"):
        df = pd.read_csv(file, index_col=0)
    elif file.lower().endswith(".tsv"):
        df = pd.read_csv(file, sep='\t', index_col=0)
    elif file.lower().endswith(".mtx"):
        data = scipy.io.mmread(file)
        if not isinstance(data, np.ndarray):
            df = pd.DataFrame.sparse.from_spmatrix(data)
        else:
            df = pd.DataFrame(data)
    elif file.lower().endswith(".gct"):
        df = pd.read_csv(file, sep='\t',skiprows=2, index_col=0)
    return df

def checkData(file, params, uncertainty=None):
    if supported(file):
        data = dataToDF(file).to_numpy()
    else:
        raise Exception("unsupported data type")

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

# TODO: add support for reading AnnData, HDF5, and R data (?)

# not implemented yet - reads HDF5 file
# we can use this for testing later 
def getRetinaSubset(n=1):
    if not (1 <= n <= 4):
        raise Exception("invalide number of subsets requested")

def nrowHelper(data):
    # if isinstance(data, str):
    #     return int(pycogaps.getFileInfo(data)["dimensions"][0])
    return data.shape[0] # assuming data is pandas dataframe

def ncolHelper(data):
    # if isinstance(data, str):
    #     return int(pycogaps.getFileInfo(data)["dimensions"][1])
    return data.shape[1] # assuming data is pandas dataframe

def getGeneNames(data, transpose):
    if transpose:
        return getSampleNames(data, False)
    # names = pycogaps.getFileInfo(data)["rowNames"]
    names = data.index.values
    
    if names.all() == None or len(names) == 0:
        return ["Gene" + str(i) for i in range(1, nrowHelper(data))]
    return names

def getSampleNames(data, transpose):
    if transpose:
        return getGeneNames(data, transpose)
    # names = pycogaps.getFileInfo(data)["colNames"]
    names = data.columns.values

    if names.all() == None or len(names) == 0:
        return ["Sample" + str(i) for i in range(1, ncolHelper(data))]
    return names

# allParams doesn't have geneNames & sampleNames yet
# currently doesn't support user-supplied param inputs
def getDimNames(data, allParams):
    # geneNames = allParams.gaps.geneNames
    # sampleNames = allParams.gaps.sampleNames

    # if allParams.gaps.geneNames == None:
    #     geneNames = getGeneNames(data, allParams.transposeData)
    # if allParams.gaps.sampleNames == None:
    #     sampleNames = getSampleNames(data, allParams.transposeData)

    if supported(data):
        data = dataToDF(data)

    geneNames = getGeneNames(data, allParams.transposeData)
    sampleNames = getSampleNames(data, allParams.transposeData)
    
    if allParams.transposeData:
        nGenes = ncolHelper(data)
        nSamples = nrowHelper(data)
    else:
        nGenes = nrowHelper(data)
        nSamples = ncolHelper(data)

    # if allParams.gaps.subsetDim == 1:
    #     nGenes = len(allParams.gaps.subsetIndices)
    #     geneNames = geneNames[allParams.subsetIndices]
    # elif allParams.gaps.subsetDim == 2:
    #     nSamples = len(allParams.subsetIndices)
    #     sampleNames = sampleNames[allParams.subsetIndices]

    if len(geneNames) != nGenes:
        raise Exception(len(geneNames), " != ", nGenes, " incorrect number of gene names given")
    if len(sampleNames) != nSamples:
        raise Exception(len(sampleNames), " != ", nSamples, " incorrect number of sample names given")

    # store processed gene/sample names directly in allParams list
    # this is an important distinction - allParams@gaps contains the
    # gene/sample names originally passed by the user, allParams contains
    # the procseed gene/sample names to be used when labeling the result
    allParams.geneNames = geneNames
    allParams.sampleNames = sampleNames
    return(allParams)