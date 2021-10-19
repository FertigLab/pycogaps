import h5py
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
import os



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

    if not supported(file):
        raise Exception("unsupported data type")
    if file.lower().endswith(".csv"):
        adata = anndata.read_csv(file)
    elif file.lower().endswith(".txt"):
        table = pd.read_table(file)
        adata = anndata.AnnData(table.iloc[:, 2:])
        adata.obs_names = table["symbol"]
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


def plot(obj, groups=None, title=None):
    """ Plots how patterns vary across samples

    Args:
        obj (CogapsResult): CogapsResult object
        groups (str list, optional): list of groups. Defaults to None.
        title (str, optional): title of plot. Defaults to None.

    Returns:
        fig: figure of plot
    """    
    
    obj = obj["anndata"]
    if groups is not None:
        if len(groups) == len(obj.var_names):
            obj.var_names = groups
        else:
            warnings.warn("length of groups does not match number of samples, aborting...")
            return
        samples = obj.var
        samplenames = list(set(obj.var_names))
        patterns = list(obj.var.columns)
        fig = plt.figure()
        ax = fig.add_subplot(111)

        for pattern in patterns:
            groupavgs = []
            for name in samplenames:
                groupavgs.append(samples.loc[name][pattern].mean())
            ax.plot(np.array(range(1, len(samplenames) + 1)), groupavgs, label=pattern)
        ax.legend()
        plt.xlabel("Groups")
        plt.ylabel("Relative Amplitude")
        plt.xticks(np.arange(1, len(samplenames) + 1), samplenames, rotation=45, ha="right")
        plt.subplots_adjust(bottom=0.15)
        if title is not None:
            ax.set_title(title)
        plt.show()
        return fig
    else:
        samples = obj.var
        nsamples = np.shape(samples)[0]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for factor in list(samples):
            ax.plot(np.array(range(1, nsamples + 1)), samples[factor], label=factor)
        ax.legend()
        plt.xlabel("Samples")
        plt.ylabel("Relative Amplitude")
        if title is not None:
            ax.set_title(title)
        plt.show()
    return fig


def patternBoxPlot(obj, groups):
    """ generate a boxplot where each subplot displays amplitudes for each group for each pattern

    Args:
        obj (CogapsResult): CogapsResult object
        groups (str list, optional): list of groups. Defaults to None.

    Returns:
        fig: figure of plot
    """    
    
    obj = obj['anndata']
    if len(groups) == len(obj.var_names):
        obj.var_names = groups
    else:
        warnings.warn("length of groups does not match number of samples, aborting...")
        return
    samples = obj.var
    samplenames = list(set(obj.var_names))
    patterns = list(obj.var.columns)
    for i in np.arange(0,4):
        thispattern = samples[patterns[i]]
        data = []
        for name in samplenames:
            data.append(thispattern.loc[name].values)
        df = pd.DataFrame(data)
        df = df.transpose()
        df.columns = samplenames
        ax = plt.subplot(2,2,i+1)
        ax.set_title(patterns[i])
        ax.set_xlabel("Groups")
        ax.set_ylabel("Amplitude")
        plt.tight_layout()
        df.boxplot(ax=ax, rot=20, fontsize=6)
    return df


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

def calcZ(object: anndata, whichMatrix):
    """ Calculates the Z-score for each element based on input mean and standard deviation matrices

    Args:
        object (anndata): Anndata result object
        whichMatrix (str): either "featureLoadings" or "sampleFactors" indicating which matrix to calculate the z-score for

    Returns:
        arr: matrix of z scores
    """    
    if whichMatrix in "sampleFactors":
        mean = object.var
        stddev = object.uns["asd"]
    elif whichMatrix in "featureLoadings":
        mean = object.obs
        stddev = object.uns["psd"]
    else:
        print('whichMatrix must be either \'featureLoadings\' or \'sampleFactors\'')
        return
    if 0 in stddev:
        print("zeros detected in the standard deviation matrix; they have been replaced by small values")
        stddev[stddev == 0] = 1 ** -6
    return mean / stddev


def reconstructGene(object: anndata, genes=None):
    """[summary]

    Args:
        object (anndata): Anndata result object
        genes (int, optional): an index of the gene or genes of interest. Defaults to None.

    Returns:
        arr: the D' estimate of a gene or set of genes
    """    
    D = np.dot(object.obs, np.transpose(object.var))
    if genes is not None:
        D = D[genes, ]
    return D


def binaryA(object, threshold, nrows="all", cluster=False):
    """ plots a binary heatmap with each entry representing whether
    that position in the A matrix has a value greater than (black)
    or lesser than (white) the specified threshold * the standard
    deviation for that element

    Args:
        object (CogapsResult): A CogapsResult object
        threshold (float): threshold to compare A/Asd
        nrows (str, optional): how many rows should be plotted (for very long
        and skinny feature matrices). Defaults to "all".
        cluster (bool, optional): True or False, whether rows should be clustered
        (results in huge black and white blocks). Defaults to False.

    Returns:
        fig: a matplotlib plot object
    """    

    object = object["anndata"]
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


def plotResiduals(object, uncertainty=None, legend=False, groups=None, ids=None):
    """ Generate a residual plot

    Args:
        object (CogapsResult): A CogapsResult object
        uncertainty (arr, optional): original SD matrix with which GAPS was run. Defaults to None.
        legend (bool, optional): Add legend to plot. Defaults to False.
        groups (list, optional): group genes for plotting. Defaults to None.
        ids (list, optional): [description]. Defaults to None.

    Returns:
        fig: matplotlib figure
    """   

    object = object["anndata"]
    # if groups is not None:
    #
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
    """ Return unit vector of length with value 1 at pos n

    Args:
        n (int): pos of value 1
        length (int): length of unit vector

    Returns:
        arr: returns numpy array
    """    
    vec = np.repeat(0, length)
    vec[n] = 1
    return vec


def patternMarkers(adata, threshold='all', lp=None, axis=1):
    """ calculate the most associated pattern for each gene

    Args:
        adata (anndata): anndata result object
        threshold (str, optional): the type of threshold to be used. The default "all" will
        distribute genes into pattern with the lowest ranking. The "cut" thresholds
        by the first gene to have a lower ranking, i.e. better fit to, a pattern.. Defaults to 'all'.
        lp (arr, optional): a vector of weights for each pattern to be used for finding
        markers. If NA markers for each pattern of the A matrix will be used.. Defaults to None.
        axis (int, optional): either 0 or 1, specifying if pattern markers should be calculated using
        the rows of the data (1) or the columns of the data (2). Defaults to 1.

    Raises:
        Exception: If threshold is not 'cut' or 'all'
        Exception: If lp length is not equal to number of patterns
        Exception: If axis is not either 0 or 1

    Returns:
        dict: A dictionary of PatternMarkers, PatternMarkerRanks, PatternMarkerScores
    """    
    if threshold.lower() not in ["cut", "all"]:
        raise Exception("threshold must be either 'cut' or 'all'")
    if lp is not None and (np.size(lp) != adata.obs.shape[1]):
        raise Exception("lp length must equal the number of patterns")
    if axis not in [1, 2]:
        raise Exception("axis must be either 0 or 1")

    if axis == 1:
        resultMatrix = adata.obs
    else:
        resultMatrix = adata.var

    # Replacing infinite with 0
    resultMatrix.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
    resultMatrix.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
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
                markerRanks[markerRanks.values[:, i] > rankCutoff[i]]).index.values

    elif threshold == "all":
        patternsByMarker = markerScores.columns[np.argmin(markerScores.values, axis=1)]
        for i in range(markerScores.shape[1]):
            markersByPattern['Pattern' + str(i + 1)] = markerScores[
                markerScores.columns[i] == patternsByMarker].index.values

    dict = {"PatternMarkers": markersByPattern, "PatternMarkerRanks": np.argsort(markerScores, axis=0),
            "PatternMarkerScores": markerScores}
    return dict


def calcCoGAPSStat(object, sets, whichMatrix='featureLoadings', numPerm=1000):
    """ calculates a statistic to determine if a pattern is enriched in a
    a particular set of measurements or samples.

    Args:
        object (CogapsResult): a CogapsResult object
        sets (list): list of sets of measurements/samples
        whichMatrix (str, optional): either "featureLoadings" or "sampleFactors" indicating which matrix
        to calculate the statistics. Defaults to 'featureLoadings'.
        numPerm (int, optional): number of permutations to use when calculatin p-value. Defaults to 1000.

    Raises:
        Exception: If sets are not a list of measurements or samples

    Returns:
        dict: dict of gene set statistics for each column of A
    """    

    if not isinstance(sets, list):
        raise Exception("Sets must be a list of either measurements or samples")

    zMatrix = calcZ(object['anndata'], whichMatrix)

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
    """ calculates the probability that a gene
    listed in a gene set behaves like other genes in the set within
    the given data set

    Args:
        object (CogapsResult): a CogapsResult object
        GStoGenes (list): list with gene sets
        numPerm (int): number of permutations for null
        Pw (arr, optional): weight on genes. Defaults to None.
        nullGenes (bool, optional): logical indicating gene adjustment. Defaults to False.

    Raises:
        Exception: If weighting is invalid

    Returns:
        dataframe: gene similiarity statistic
    """    
    featureLoadings = toNumpy(object['GapsResult'].Amean)
    
    adata = object['anndata']

    if Pw is None:
        Pw = np.ones(featureLoadings.shape[1])
    gsStat = calcCoGAPSStat(object, GStoGenes, numPerm=numPerm)
    gsStat = gsStat['GSUpreg'].values.T
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
    """ Computes the p-value for gene set membership using the CoGAPS-based
    statistics developed in Fertig et al. (2012).  This statistic refines set
    membership for each candidate gene in a set specified in \code{GSGenes} by
    comparing the inferred activity of that gene to the average activity of the
    set.

    Args:
        object (CogapsResult): a CogapsResult object
        GStoGenes (list): list with gene sets
        numPerm (int, optional): number of permutations for null. Defaults to 500.
        Pw ([type], optional): weight on genes. Defaults to None.
        PwNull (bool, optional): logical indicating gene adjustment. Defaults to False.

    Returns:
        arr: A vector of length GSGenes containing the p-values of set membership
        for each gene containined in the set specified in GSGenes.
    """    

    featureLoadings = toNumpy(object['GapsResult'].Amean)
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



def plotPatternMarkers(data, patternmarkers=None, groups = None, patternPalette=None,
                       samplePalette=None, colorscheme="coolwarm",
                       colDendrogram=True, rowDendrogram=False, scale="row", legend_pos=None):
    """ Plots pattern markers of most associated pattern for each gene.

    Args:
        data (anndata):  an anndata object, which should be your original data annotated with CoGAPS results
        patternmarkers (list, optional): list of markers for each pattern, as determined by the "patternMarkers(data)" function. Defaults to None.
        groups (list, optional): list of genes to group. Defaults to None.
        patternPalette (list, optional): a list of colors to be used for each pattern. 
        if None, colors will be set automatically. Defaults to None.
        samplePalette (list, optional): a list of colors to be used for each sample. 
        if None, colors will be set automatically. Defaults to None.
        colorscheme (str, optional): string indicating which color scheme should be used within the heatmap. 
        more options at https://seaborn.pydata.org/tutorial/color_palettes.html. Defaults to "coolwarm".
        colDendrogram (bool, optional):  Whether or not to draw a column dendrogram, default true. Defaults to True.
        rowDendrogram (bool, optional): Whether or not to draw a row dendrogram, default false. Defaults to False.
        scale (str, optional): whether you want data to be scaled by row, column, or none. Defaults to "row".
        legend_pos (str, optional): string indicating legend position, or none (no legend). Defaults to None.
    
    Returns:
        fig: a clustergrid instance
    """    

    data = data["anndata"]
    if patternmarkers is None:
        patternmarkers=patternMarkers(data)
    if samplePalette is None:
        if groups is None:
            # color for each sample
            samplePalette=sns.color_palette("Spectral", np.shape(data)[1])
        else:
            # color for each group
            samplePalette = sns.color_palette("Spectral", len(set(groups)))
            palette = []
            groupkeys = list(set(groups))
            grplst = list(groups)
            for i in range(len(groupkeys)):
                palette = np.concatenate((palette, np.repeat(mpl.colors.to_hex(samplePalette[i]), grplst.count(groupkeys[i]))))
            samplePalette = palette
    if patternPalette is None:
        palette = []
        patternkeys = list(patternmarkers["PatternMarkers"].keys())
        thiscmap = sns.color_palette("Spectral", len(patternkeys))
        for i in range(len(patternkeys)):
            palette = np.concatenate((palette, np.repeat(mpl.colors.to_hex(thiscmap[i]), len(patternmarkers["PatternMarkers"][patternkeys[i]]))))
        patternPalette = palette
    elif patternPalette is not None:
        palette = []
        patternkeys = list(patternmarkers["PatternMarkers"].keys())
        for i in range(len(patternkeys)):
            palette = np.concatenate((palette, np.repeat(patternPalette[i], len(patternmarkers["PatternMarkers"][patternkeys[i]]))))
        patternPalette = palette
    if groups is not None:
        top = []
        markers = patternmarkers["PatternMarkers"]
        keys=markers.keys()
        for key in keys:
            top.append(markers[key][1:15])
        top=np.transpose(top)
        # top.columns = patterns[1:10]
        markers = [item for sublist in top for item in sublist]
        markers = [x for x in markers if str(x) != 'nan']
        if len(groups) == len(data.var_names):
            data.var_names = groups
        else:
            warnings.warn("length of groups does not match number of samples, aborting...")
            return
        samples = data.var
        markermatrix = []
        for group in groups:
            grplst = []
            for marker in markers:
                print(group, marker)
                grplst.append(np.average(data[marker, group].X).tolist())
            markermatrix.append(grplst)


        samplenames = list(set(data.var_names))
        patterns = list(data.var.columns)
    else:
        markers = np.concatenate(list(patternmarkers["PatternMarkers"].values()))
        plotinfo = data[data.obs_names.isin(markers)]
        plotdata = plotinfo.X
        markerlabels = plotinfo.obs_names
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
    plotdata_z.columns = samplelabels
    plotdata_z.index = markerlabels
    plotdata_z.replace([np.inf, -np.inf, np.nan], 0, inplace=True)

    hm = sns.clustermap(plotdata_z, cmap=colorscheme, row_cluster=rowDendrogram, col_cluster=colDendrogram,
                        row_colors=patternPalette, col_colors=samplePalette, cbar_pos=legend_pos)
    plt.show()
    return hm


def plotUMAP(result, genes_in_rows=True):
    """ Create a UMAP plot

    Args:
        result (anndata or CogapsResult): An anndata object of result or CogapsResult object
        genes_in_rows (bool, optional): Scanpy needs genes in columns, cells in rows. Defaults to True.
    """    
    print("not implemented")
    if not isinstance(result, anndata):
        result=result["anndata"]
    if genes_in_rows:
        # scanpy needs genes in columns, cells in rows
        result = result.transpose()
    import scanpy as sc
    # set up environment
    patterns = list(result.obs.columns)
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')
    # result.var_names_make_unique()
    # filter genes and cells
    sc.pl.highest_expr_genes(result, n_top=20, )
    sc.pp.filter_cells(result, min_genes=200)
    sc.pp.filter_genes(result, min_cells=3)
    sc.pp.log1p(result)
    sc.pp.highly_variable_genes(result, min_mean=0.0125, max_mean=3, min_disp=0.5)
    result = result[:, result.var.highly_variable]
    sc.pp.scale(result, max_value=10)
    sc.tl.pca(result, svd_solver='arpack')
    sc.pp.neighbors(result)
    sc.tl.umap(result)
    sc.pl.umap(result, color=patterns)





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
