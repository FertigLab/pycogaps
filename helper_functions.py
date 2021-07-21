from pycogaps import *
import numpy as np
import pandas as pd
import scipy.io
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
# import colorspacious
import pkg_resources  # part of setuptools
from numpy import genfromtxt
import anndata


def supported(file):
    return file.lower().endswith((".tsv", ".csv", ".mtx", ".h5ad", ".hdf"))  # currently gct not supported w/ anndata


def checkData(adata, params, uncertainty=None):
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


def toAnndata(file):
    if not supported(file):
        raise Exception("unsupported data type")
    if file.lower().endswith(".csv"):
        adata = anndata.read_csv(file)
    elif file.lower().endswith(".tsv"):
        adata = anndata.read_csv(file, delimeter='\t')
    elif file.lower().endswith(".mtx"):
        adata = anndata.read_mtx(file)
    elif file.lower().endswith(".h5ad"):
        adata = anndata.read_h5ad(file)
    elif file.lower().endswith(".hdf"):
        adata = anndata.read_hdf(file)
    # elif file.lower().endswith(".gct")

    if scipy.sparse.issparse(adata.X):
        adata.X = (adata.X).toarray()

    return adata


# not implemented yet - reads HDF5 file
# we can use this for testing later 
def getRetinaSubset(n=1):
    if not (1 <= n <= 4):
        raise Exception("invalide number of subsets requested")


def nrowHelper(data):
    return data.shape[0]  # assuming data is pandas dataframe


def ncolHelper(data):
    return data.shape[1]  # assuming data is pandas dataframe


def getGeneNames(data, transpose):
    if transpose:
        return getSampleNames(data, False)
    names = data.index.values

    if names.all() == None or len(names) == 0:
        return ["Gene" + str(i) for i in range(1, nrowHelper(data))]
    return names


def getSampleNames(data, transpose):
    if transpose:
        return getGeneNames(data, transpose)
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

    if not supported(data):
        raise Exception("unsupported data type")

    data = toAnndata(data).X

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
    nfeatures = np.shape(toNumpy(obj.Amean))[0]
    nsamples = np.shape(toNumpy(obj.Pmean))[0]
    npatterns = np.shape(toNumpy(obj.Pmean))[1]
    print("GapsResult result object with ", nfeatures, " features and ", nsamples, " samples")
    print(npatterns, " patterns were learned")
    return


def plot(obj: GapsResult):
    samples = toNumpy(obj.Pmean)
    nsamples = np.shape(samples)[0]
    nfactors = np.shape(samples)[1]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(nfactors):
        ax.plot(np.array(range(1, nsamples + 1)), samples[:, i], label="Pattern " + str(i + 1))
    ax.legend()
    plt.xlabel("Samples")
    plt.ylabel("Relative Amplitude")
    plt.xlim(0, nsamples + 1)
    plt.ylim(0, np.argmax(samples) * 1.1)
    plt.show()


def getFeatureLoadings(object: GapsResult):
    return object.Amean


def getAmplitudeMatrix(object):
    return object.Amean


def getSampleFactors(object):
    return object.Pmean


def getPatternMatrix(object):
    return object.Pmean


def getMeanChiSq(object: GapsResult):
    return object.meanChiSq


def getVersion():
    version = pkg_resources.require("pycogaps")[0].version
    print("pycogaps version ", version)
    return version


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
        mean = toNumpy(object.Pmean)
        stddev = toNumpy(object.Psd)
    elif whichMatrix in "featureLoadings":
        mean = toNumpy(object.Amean)
        stddev = toNumpy(object.Asd)
    else:
        print('whichMatrix must be either \'featureLoadings\' or \'sampleFactors\'')
        return
    if 0 in stddev:
        print("zeros detected in the standard deviation matrix; they have been replaced by small values")
        stddev[stddev == 0] = 1 ** -6
    return mean / stddev


def reconstructGene(object: GapsResult, genes=None):
    D = np.dot(toNumpy(object.Amean), np.transpose(toNumpy(object.Pmean)))
    if genes is not None:
        D = D[genes,]
    return D


def binaryA(object: GapsResult, threshold, nrows="all"):
    """
    plots a binary heatmap with each entry representing whether
    that position in the A matrix has a value greater than (black)
    or lesser than (white) the specified threshold * the standard
    deviation for that element
    @param object: GapsResult object
    @param threshold: threshold to compare to A/Asd
    @param nrows: how many rows should be plotted (for very long
    and skinny feature matrices)
    @return: matplotlib plot object
    """
    binA = calcZ(object, whichMatrix="featureLoadings")
    if nrows != "all":
        binA = binA[1:nrows, :]
    overthresh = binA > threshold
    underthresh = binA < threshold
    binA[overthresh] = 1
    binA[underthresh] = 0
    plt.matshow(binA, cmap='hot', interpolation='nearest')
    plt.show()
    return plt


def plotResiduals(object: GapsResult, data, uncertainty):
    """
    generate a residual plot
    @param object: GapsResult object
    @param data: original data matrix on which GAPS was run
    @param uncertainty: original SD matrix with which GAPS was run
    @return: matplotlib plot object
    """
    data = np.array(data)
    if uncertainty is None:
        uncertainty = np.where(data * 0.1 > 0.1, data * 0.1, 0.1)
    uncertainty = np.array(uncertainty)

    M = reconstructGene(object)
    residual = (data - M) / uncertainty
    plt.matshow(residual, cmap='hot', interpolation='nearest')
    plt.show()
    return plt


def unitVector(n, length):
    vec = np.repeat(0, length)
    vec[n] = 1
    return vec


def patternMarkers(adata, threshold='all', lp=None, axis=1):
    if threshold.lower() not in ["cut", "all"]:
        raise Exception("threshold must be either 'cut' or 'all'")
    if lp is not None and (np.size(lp) != adata.obs.shape[1]):
        raise Exception("lp length must equal the number of patterns")
    if axis not in [1, 2]:
        raise Exception("axis must be either 1 or 2")

    if axis == 1:
        resultMatrix = adata.obs
    else:
        resultMatrix = adata.var

    row_max = np.nanmax(resultMatrix.values, axis=1, keepdims=True)
    row_max = np.where(row_max == 0, 1, row_max)

    normedMatrix = resultMatrix / row_max

    if lp is not None:
        markerScores = pd.DataFrame(np.sqrt(np.sum((normedMatrix.values - lp) ** 2, axis=1)), index=normedMatrix.index)
        markersByPattern = markerScores.sort_values(0).index.values
        dict = {"PatternMarkers": markersByPattern, "PatternMarkerRanks": np.argsort(markerScores, axis=0),
                "PatternMarkerScores": markerScores}
        return dict

    markerScores_arr = np.empty_like(normedMatrix)
    for i in range(normedMatrix.shape[1]):
        lp = unitVector(i, normedMatrix.shape[1])
        markerScores_arr[:, i] = np.sqrt(np.sum((normedMatrix.values - lp) ** 2, axis=1))

    markerScores = pd.DataFrame(markerScores_arr, index=normedMatrix.index, columns=normedMatrix.columns)

    markerRanks = pd.DataFrame(np.argsort(markerScores.values, axis=0), index=markerScores.index,
                               columns=markerScores.columns)

    rankCutoff = np.empty(markerRanks.shape[1])
    markersByPattern = {}
    if threshold == "cut":
        for i in range(markerRanks.shape[1]):
            patternRank = markerRanks.values[:, i]
            rankCutoff[i] = np.max(patternRank[patternRank == np.amin(markerRanks, axis=1)])
            markersByPattern['Pattern' + str(i + 1)] = (
                markerRanks[markerRanks.values[:, i] <= rankCutoff[i]]).index.values

    elif threshold == "all":
        patternsByMarker = markerScores.columns[np.argmin(markerScores.values, axis=1)]
        for i in range(markerScores.shape[1]):
            markersByPattern['Pattern' + str(i + 1)] = markerScores[
                markerScores.columns[i] == patternsByMarker].index.values

    dict = {"PatternMarkers": markersByPattern, "PatternMarkerRanks": np.argsort(markerScores, axis=0),
            "PatternMarkerScores": markerScores}
    return dict


def calcCoGAPSStat(object):
    print("Not yet implemented")
    return


def calcGeneGSStat(object, GStoGenes, numPerm, Pw, nullGenes):
    print("Not yet implemented")
    return


def computeGeneGSProb(object, GStoGenes, numPerm, Pw, PwNull):
    print("Not yet implemented")
    return


def plotPatternMarkers(object: GapsResult, data, patternmarkers, patternPalette, sampleNames,
                       samplePalette=None, heatmapCol="bluered",
                       colDenogram=True, scale="row"):
    """

    @param object: GapsResult object
    @param data: original data in matrix format
    @param patternPalette: vector indicating which color should be used for each pattern
    @param sampleNames: names with which samples should be labelled
    @param samplePalette: vector indicating which color should be used for each sample
    @param heatmapCol: pallelet giving color scheme for heatmap
    @param colDenogram: logical indicating whether to display sample dendrogram
    @param scale: character indicating if the values should be centered and scaled in
    the row direction, the column direction, or none. the default is "row"
    @return:
    """
    if samplePalette is None:
        samplePalette = np.repeat("black", np.shape(data)[1])
    # this is the equivalent of length(unlist(patternmarkers[...])) in the R code
    nummarkers = len(np.concatenate(list(patternmarkers["PatternMarkers"].values())))
    patternCols = []
    i = 0
    for pattern in patternmarkers["PatternMarkers"]:
        patternlength = len(patternmarkers["PatternMarkers"][pattern])
        tmp = np.repeat(patternPalette[i], patternlength)
        patternCols.append(tmp)
        i+=1

    
    return


# convert matrix object to numpy array
def toNumpy(matrix):
    all_vecdata = np.empty((matrix.nRow(), matrix.nCol()))
    for i in range(matrix.nCol()):
        vector = matrix.getCol(i)
        vecdata = []
        for j in range(vector.size()):
            vecdata.append(getElement(vector, j))
        all_vecdata[:, i] = vecdata
    return all_vecdata
