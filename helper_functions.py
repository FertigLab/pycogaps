from pycogaps import *
import numpy as np
import pandas as pd
import scipy.io
import matplotlib as mpl
import matplotlib.pyplot as plt
import pkg_resources  # part of setuptools
import anndata
import seaborn as sns
from scipy.stats import zscore
import warnings
from PyCoGAPS import *


def supported(file):
    return file.lower().endswith((".tsv", ".csv", ".mtx", ".h5ad", ".h5"))  # currently gct not supported w/ anndata


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
    elif file.lower().endswith(".h5"):
        adata = anndata.read_hdf(file, "counts") # change to user supplied key
    # elif file.lower().endswith(".gct")

    if scipy.sparse.issparse(adata.X):
        adata.X = (adata.X).toarray()

    return adata


# not implemented yet - reads HDF5 file
# we can use this for testing later 
def getRetinaSubset(n=1):
    if not (1 <= n <= 4):
        raise Exception("invalid number of subsets requested")


def nrowHelper(data):
    return data.shape[0]


def ncolHelper(data):
    return data.shape[1]


def getGeneNames(data, transpose):
    if transpose:
        return getSampleNames(data, False)
    names = data.obs_names

    if names.all() == None or len(names) == 0:
        return ["Gene" + str(i) for i in range(1, nrowHelper(data))]
    return names


def getSampleNames(data, transpose):
    if transpose:
        return getGeneNames(data, False)
    names = data.var_names

    if names.all() == None or len(names) == 0:
        return ["Sample" + str(i) for i in range(1, ncolHelper(data))]
    return names



def getDimNames(data, allParams):
    # support both path and anndata object as data input
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
        geneNames = geneNames[allParams.coparams['subsetIndices']]
    elif allParams.coparams['subsetDim'] == 2:
        nSamples = len(allParams.coparams['subsetIndices'])
        sampleNames = sampleNames[allParams.coparams['subsetIndices']]

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


def show(obj: anndata):
    nfeatures = obj.n_obs
    nsamples = obj.n_var
    npatterns = len(obj.obs_keys())
    print("GapsResult result object with ", nfeatures, " features and ", nsamples, " samples")
    print(npatterns, " patterns were learned")
    return


def plot(obj: anndata):
    samples = obj.var
    nsamples = np.shape(samples)[0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for factor in list(samples):
        ax.plot(np.array(range(1, nsamples + 1)), samples[factor], label=factor)
    ax.legend()
    plt.xlabel("Samples")
    plt.ylabel("Relative Amplitude")
    plt.xlim(0, nsamples + 1)
    plt.ylim(0, np.argmax(samples) * 1.1)
    plt.show()


def getFeatureLoadings(object: anndata):
    return object.obs


def getAmplitudeMatrix(object):
    return object.obs


def getSampleFactors(object):
    return object.var


def getPatternMatrix(object):
    return object.var

# TODO need to access chisq through anndata
def getMeanChiSq(object: GapsResult):
    return object.meanChiSq


def getVersion():
    version = pkg_resources.require("pycogaps")[0].version
    print("pycogaps version ", version)
    return version

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


def calcZ(object: anndata, whichMatrix):
    if whichMatrix in "sampleFactors":
        mean = object.var
        stddev = object.uns["psd"]
    elif whichMatrix in "featureLoadings":
        mean = object.obs
        stddev = object.uns["asd"]
    else:
        print('whichMatrix must be either \'featureLoadings\' or \'sampleFactors\'')
        return
    if 0 in stddev:
        print("zeros detected in the standard deviation matrix; they have been replaced by small values")
        stddev[stddev == 0] = 1 ** -6
    return mean / stddev


def reconstructGene(object: anndata, genes=None):
    D = np.dot(object.obs, np.transpose(object.var))
    if genes is not None:
        D = D[genes, ]
    return D


def binaryA(object: anndata, threshold, nrows="all", cluster=False):
    """
    plots a binary heatmap with each entry representing whether
    that position in the A matrix has a value greater than (black)
    or lesser than (white) the specified threshold * the standard
    deviation for that element
    @param cluster: True or False, whether rows should be clustered
    (results in huge black and white blocks)
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
    if cluster:
        hm = sns.clustermap(binA, cbar_pos=None)
    else:
        hm = sns.heatmap(binA, cbar=False)
    plt.show()
    return hm


def plotResiduals(object: anndata, uncertainty=None, legend=False):
    """
    generate a residual plot
    @param object: AnnData object
    @param uncertainty: original SD matrix with which GAPS was run
    @return: matplotlib plot object
    """
    rawdata = object.X
    if uncertainty is None:
        uncertainty = np.where(rawdata * 0.1 > 0.1, rawdata * 0.1, 0.1)
    uncertainty = np.array(uncertainty)

    markerlabels = object.obs_names
    samplelabels = object.var_names
    M = reconstructGene(object)
    residual = (rawdata - M) / uncertainty
    residual = pd.DataFrame(residual, columns=samplelabels, index=markerlabels)
    hm = sns.heatmap(residual, cmap="Spectral", cbar=legend)
    plt.show()
    return hm


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


def calcCoGAPSStat(object, sets, whichMatrix='featureLoadings', numPerm=1000):

    if not isinstance(sets, list):
        raise Exception("Sets must be a list of either measurements of samples")

    zMatrix = calcZ(object['GapsResult'], whichMatrix)

    pattern_labels = (object['anndata'].obs).columns

    zMatrix = pd.DataFrame(zMatrix, index=object['anndata'].obs_names, columns=pattern_labels)
    pvalUpReg = []

    lessThanCount = np.zeros(zMatrix.shape[1])
    actualZScore = np.mean(zMatrix.loc[sets,:].values, axis=0)
    for n in range(numPerm):
        permutedIndices = np.random.choice(np.arange(1, zMatrix.shape[0]), size=len(sets), replace=False)
        permutedZScore = np.mean(zMatrix.iloc[permutedIndices,:].values, axis=0)
        lessThanCount = lessThanCount + (actualZScore < permutedZScore)
    pvalUpReg.append(lessThanCount / numPerm)

    pvalUpReg = np.array(pvalUpReg)
    pvalDownReg = 1 - pvalUpReg
    activityEstimate = 1 - 2 * pvalUpReg
    
    dict = {'twoSidedPValue': pd.DataFrame((np.maximum(np.minimum(pvalDownReg, pvalUpReg), 1 / numPerm)).T, index=pattern_labels),
        'GSUpreg': pd.DataFrame(pvalUpReg.T, index=pattern_labels),
        'GSDownreg': pd.DataFrame(pvalDownReg.T, index=pattern_labels),
        'GSActEst': pd.DataFrame(activityEstimate.T, index=pattern_labels)}
    
    return dict


def calcGeneGSStat(object, GStoGenes, numPerm, Pw=None, nullGenes=False):
    featureLoadings = toNumpy(object['GapsResult'].Amean)
    
    adata = object['anndata']

    if Pw is None:
        Pw = np.ones(featureLoadings.shape[1])
    gsStat = calcCoGAPSStat(object, GStoGenes, numPerm=numPerm)
    gsStat =  gsStat['GSUpreg'].values.T
    gsStat = -np.log(gsStat)

    if not np.isnan(Pw).all():
        if np.size(Pw) != gsStat.shape[1]:
            raise Exception('Invalid weighting')
        gsStat = gsStat*Pw
    
    stddev = toNumpy(object['GapsResult'].Asd)
    if 0 in stddev:
        print("zeros detected in the standard deviation matrix; they have been replaced by small values")
        stddev[stddev == 0] = 1 ** -6
    stddev = pd.DataFrame(stddev, index=adata.obs_names, columns=(adata.obs).columns)

    featureLoadings = pd.DataFrame(featureLoadings, index=adata.obs_names, columns=(adata.obs).columns)

    if nullGenes:
        ZD = featureLoadings.loc[(featureLoadings.index).difference(GStoGenes),:].values / stddev.loc[(featureLoadings.index).difference(GStoGenes),:].values
    else:
        ZD = featureLoadings.loc[GStoGenes,:].values / stddev.loc[GStoGenes,:].values

    ZD_apply = np.multiply(ZD, gsStat)
    ZD_apply = np.sum(ZD_apply, axis=1)

    outStats = ZD_apply / np.sum(gsStat)
    outStats = outStats / np.sum(ZD, axis=1)
    outStats[np.argwhere(np.sum(ZD, axis=1) < 1e-6)] = 0

    if np.sum(gsStat) < 1e-6:
        return 0

    if nullGenes:
        outStats = pd.DataFrame(outStats, index=(featureLoadings.index).difference(GStoGenes))
    else:
        outStats = pd.DataFrame(outStats, index=GStoGenes)

    return outStats


def computeGeneGSProb(object, GStoGenes, numPerm=500, Pw=None, PwNull=False):

    featureLoadings = toNumpy(getFeatureLoadings(object['GapsResult']))
    adata = object['anndata']
    
    if Pw is None:
        Pw = np.ones(featureLoadings.shape[1])

    geneGSStat = calcGeneGSStat(object, Pw=Pw, GStoGenes=GStoGenes, numPerm=numPerm).values

    if PwNull:
        permGSStat = calcGeneGSStat(object, GStoGenes=GStoGenes, numPerm=numPerm, Pw=Pw, nullGenes=True).values
    else:
        permGSStat = calcGeneGSStat(object, GStoGenes=GStoGenes, numPerm=numPerm, nullGenes=True).values

    finalStats = np.empty(len(GStoGenes))
    for i in range(len(GStoGenes)):
        finalStats[i] = np.size(np.argwhere(permGSStat > geneGSStat[i])) / np.size(permGSStat)

    finalStats = pd.DataFrame(finalStats, index=GStoGenes)

    return finalStats



def plotPatternMarkers(data: anndata, patternmarkers=None, patternPalette=None,
                       samplePalette=None, colorscheme="coolwarm",
                       colDendrogram=True, rowDendrogram=False, scale="row", legend_pos=None):
    """
    @param data: an anndata object, which should be your original data annotated with CoGAPS results
    @param patternmarkers: list of markers for each pattern, as determined by the "patternMarkers(data)" function
    @param patternPalette: a list of colors to be used for each pattern. if None, colors will be set automatically
    @param samplePalette: a list of colors to be used for each sample. if None, colors will be set automatically
    @param colorscheme: string indicating which color scheme should be used within the heatmap. more options at https://seaborn.pydata.org/tutorial/color_palettes.html
    @param colDendrogram: Whether or not to draw a column dendrogram, default true
    @param rowDendrogram: Whether or not to draw a row dendrogram, default false
    @param scale: whether you want data to be scaled by row, column, or none. default is row
    @param legend_pos: string indicating legend position, or none (no legend). default is none
)    @return: a clustergrid instance
    """
    if patternmarkers is None:
        patternmarkers=patternMarkers(data)
    if samplePalette is None:
        samplePalette=sns.color_palette("Spectral", np.shape(data)[1])
    if patternPalette is None:
        thiscmap = sns.color_palette("Spectral")
        palette = []
        patternkeys = list(patternmarkers["PatternMarkers"].keys())
        for i in range(len(patternkeys)):
            palette = np.concatenate((palette, np.repeat(mpl.colors.to_hex(thiscmap[i]), len(patternmarkers["PatternMarkers"][patternkeys[i]]))))
        patternPalette = palette
    elif patternPalette is not None:
        palette = []
        patternkeys = list(patternmarkers["PatternMarkers"].keys())
        for i in range(len(patternkeys)):
            palette = np.concatenate((palette, np.repeat(patternPalette[i], len(patternmarkers["PatternMarkers"][patternkeys[i]]))))
        patternPalette = palette

    markers = np.concatenate(list(patternmarkers["PatternMarkers"].values()))
    plotdata = data[markers].X
    markerlabels = data[markers].obs_names
    samplelabels = data[markers].var_names

    if scale not in ["row", "column", "none"]:
        warnings.warn("warning: scale must be one of \"row\", \"column\", or \"none\". data will not be scaled in "
                      "this plot")
    if scale == "row":
        t = np.transpose(pd.DataFrame(plotdata))
        z = zscore(t)
        plotdata_z = np.transpose(z)
    elif scale == "column":
        plotdata_z = zscore(pd.DataFrame(plotdata))
    else:
        plotdata_z = pd.DataFrame(plotdata)

    plotdata_z = pd.DataFrame(plotdata_z, columns=samplelabels, index=markerlabels)

    hm = sns.clustermap(plotdata_z, cmap=colorscheme, row_cluster=rowDendrogram, col_cluster=colDendrogram, row_colors=patternPalette, col_colors=samplePalette, cbar_pos=legend_pos)
    plt.show()
    return hm


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
