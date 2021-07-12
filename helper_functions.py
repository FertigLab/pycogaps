from pycogaps import *
import numpy as np
import pandas as pd
import scipy.io
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
import colorspacious
from colorspacious import cspace_converter
from collections import OrderedDict


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
        df = pd.read_csv(file, sep='\t', skiprows=2, index_col=0)
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
    return data.shape[0]  # assuming data is pandas dataframe


def ncolHelper(data):
    # if isinstance(data, str):
    #     return int(pycogaps.getFileInfo(data)["dimensions"][1])
    return data.shape[1]  # assuming data is pandas dataframe


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
    return (allParams)


# TODO: implement CogapsResults helper functions

def createCogapsResult(object: GapsResult, params: GapsParameters):
    print("Not yet implemented")
    return


def show(obj: GapsResult):
    nfeatures = obj.Amean.nRow()
    nsamples = obj.Pmean.nRow()
    npatterns = obj.Pmean.nCol()
    print("GapsResult result object with ", nfeatures, " features and ", nsamples, " samples")
    print(npatterns, " patterns were learned")
    return


def plot(obj: GapsResult):
    nsamples = obj.Pmean.nRow()
    nfactors = obj.Pmean.nCol()
    mtx = obj.Pmean
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(nfactors):
        vector = mtx.getCol(i)
        vecdata = []
        for j in range(vector.size()):
            vecdata.append(getElement(vector, j))
        ax.plot(np.array(range(1, nsamples + 1)), np.array(vecdata), label="Pattern " + str(i + 1))
    ax.legend()
    plt.xlabel("Samples")
    plt.ylabel("Relative Amplitude")
    plt.show()


def getFeatureLoadings(object: GapsResult):
    return object.Amean


def getAmplitudeMatrix(object):
    return object.Amean


def getSampleFactors(object):
    return object.Pmean


def getPatternMatrix(object):
    return object.Pmean


def getetMeanChiSq(object: GapsResult):
    return object.meanChiSq


def getVersion(object):
    print("Not yet implemented")
    return


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


def calcZ(object: GapsResult, whichMatrix):
    if whichMatrix in "sampleFactors":
        mean = object.Pmean
        stddev = object.Psd
    elif whichMatrix in "featureLoadings":
        mean = object.Amean
        stddev = object.Asd
    else:
        print('whichMatrix must be either \'featureLoadings\' or \'sampleFactors\'')
        return
    nonzero = containsZeros(stddev)
    print("nonzero count", nonzero)
    if nonzero > 0:
        print("zeros detected in the standard deviation matrix; they have been replaced by small values")
        stddev = replaceZeros(stddev)
    return divideMatrices(mean, stddev)


def reconstructGene(object:GapsResult, genes):
    D = multiplyMatrices(object.Amean, transposeMatrix(object.Pmean))
    if genes is not None:
        # TODO: subset genes... i'm confused as to what's supposed to be happening here
        return D

#TODO: figure out what this one actully does lol
def binaryA(object:GapsResult, threshold):
    if calcZ(object) > threshold:
        binA=1
    else:
        binA = 0;

    a = np.random.random((16, 16))
    plt.imshow(a, cmap='hot', interpolation='nearest')
    plt.show()
    return


def plotResiduals(object):
    print("Not yet implemented")
    return


def unitVector():
    print("Not yet implemented")
    return


def patternMarkers(object, threshold, lp, axis):
    print("Not yet implemented")
    return


def calcCoGAPSStat(object):
    print("Not yet implemented")
    return


def calcGeneGSStat(object, GStoGenes, numPerm, Pw, nullGenes):
    print("Not yet implemented")
    return


def computeGeneGSProb(object, GStoGenes, numPerm, Pw, PwNull):
    print("Not yet implemented")
    return


def plotPatternMarkers(object, data, patternPalette, sampleNames,
                       samplePalette=None, heatmapCol="bluered",
                       colDenogram=True, scale="row"):
    print("Not yet implemented")
    return
