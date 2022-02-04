import pandas as pd
import numpy as np
import anndata
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore

def plot(obj, groups=None, title=None, fn=""):
    """ Plots how patterns vary across samples

    Args:
        obj (CogapsResult): CogapsResult object
        groups (str list, optional): list of groups. Defaults to None.
        title (str, optional): title of plot. Defaults to None.

    Returns:
        fig: figure of plot
    """    
    
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
        else:
            ax.set_title('Patterns over Samples')
        
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
        else:
            ax.set_title('Patterns over Samples')
    plt.savefig("{}_plot.png".format(fn))
    plt.show()
    return fig


def patternBoxPlot(obj, groups, fn=""):
    """ generate a boxplot where each subplot displays amplitudes for each group for each pattern

    Args:
        obj (CogapsResult): CogapsResult object
        groups (str list): list of groups. 

    Returns:
        fig: figure of plot
    """    
    
    # obj = obj['anndata']
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
    plt.savefig("{}_patternBoxPlot.png".format(fn))
    return df


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


def binaryA(object, threshold, nrows="all", cluster=False, fn=""):
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

    # object = object["anndata"]
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
    plt.title('Binary Heatmap')
    plt.savefig("{}_binaryA.png".format(fn))
    plt.show()
    return hm


def plotResiduals(object, uncertainty=None, legend=False, groups=None, ids=None, fn=""):
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

    # object = object["anndata"]
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
    plt.title('Residuals Plot')
    plt.savefig("{}_residualsPlot.png".format(fn))
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
        otherMatrix = adata.var
    else:
        resultMatrix = adata.var
        otherMatrix = adata.obs

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
        def simplicityGENES(As, Ps):
            # Using MinMaxScaler
            As.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
            Ps.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
            from sklearn.preprocessing import MinMaxScaler
            scaler = MinMaxScaler(feature_range=(0,1))
            # Stack everything into a single column to scale by the global min / max
            tmp = Ps.to_numpy().reshape(-1, 1)
            scaled2 = scaler.fit_transform(tmp).reshape(Ps.shape)
            pscale = pd.DataFrame(scaled2, index=Ps.index, columns=Ps.columns)
            # TODO: figure out how to translate this: As < - sweep(As, 2, pscale, FUN="*")
            # pmax is the most significant pattern for each gene
            # Arowmax is the A matrix with every element divided by the max in the row
            pmax = np.nanmax(As.values, axis=1, keepdims=True)
            Arowmax = As / pmax

            ssl = pd.DataFrame().reindex_like(As.drop_duplicates())
            import math
            for i in np.arange(As.shape[1]):
                lp = np.repeat(0, As.shape[1])
                lp[i] = 1
                def stat(x):
                    return (math.sqrt(np.matmul(np.transpose(x-lp), (x-lp))))

                ssl.stat = Arowmax.drop_duplicates().apply(func=stat, axis=1)
                order = np.argsort(ssl.stat)
                ssl["Pattern"+str(i+1)] = order.values

            return ssl[ssl >= 0]

        simGenes = simplicityGENES(As=resultMatrix, Ps=otherMatrix)
        nP = simGenes.shape[1]


        for i in np.arange(nP):
            # order gene names by significance for this pattern
            pname = "Pattern"+str(i+1)
            sortSim = simGenes[pname].sort_values().index
            sortedGenes = simGenes.loc[sortSim, :]
            globalmins = sortedGenes.min(axis=1)
            thispattern = simGenes.loc[sortSim, pname]

            geneThresh = int(thispattern[thispattern > globalmins].min())

            markerGenes = sortSim[1:geneThresh]
            markersByPattern[pname] = markerGenes.values

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

    zMatrix = calcZ(object, whichMatrix)

    pattern_labels = (object.obs).columns

    zMatrix = pd.DataFrame(zMatrix, index=object.obs_names, columns=pattern_labels)
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
    featureLoadings = object.obs
    
    # adata = object

    if Pw is None:
        Pw = np.ones(featureLoadings.shape[1])
    gsStat = calcCoGAPSStat(object, GStoGenes, numPerm=numPerm)
    gsStat = gsStat['GSUpreg'].values.T
    gsStat = -np.log(gsStat)

    if not np.isnan(Pw).all():
        if np.size(Pw) != gsStat.shape[1]:
            raise Exception('Invalid weighting')
        gsStat = gsStat*Pw
    
    stddev = object.uns['asd']
    if 0 in stddev:
        print("zeros detected in the standard deviation matrix; they have been replaced by small values")
        stddev[stddev == 0] = 1 ** -6
    stddev = pd.DataFrame(stddev, index=object.obs_names, columns=(object.obs).columns)

    featureLoadings = pd.DataFrame(featureLoadings, index=object.obs_names, columns=(object.obs).columns)

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

    featureLoadings = object.obs
    # adata = object['anndata']
    
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
                       colDendrogram=True, rowDendrogram=False, scale="row", legend_pos=None, fn=""):
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

    # data = data["anndata"]
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
        plotdata_z = pd.DataFrame(np.transpose(z))
    elif scale == "column":
        plotdata_z = pd.DataFrame(zscore(pd.DataFrame(plotdata)))
    else:
        plotdata_z = pd.DataFrame(plotdata)
    plotdata_z.columns = samplelabels
    plotdata_z.index = markerlabels
    plotdata_z.replace([np.inf, -np.inf, np.nan], 0, inplace=True)

    hm = sns.clustermap(plotdata_z, cmap=colorscheme, row_cluster=rowDendrogram, col_cluster=colDendrogram,
                        row_colors=patternPalette, col_colors=samplePalette, cbar_pos=legend_pos)
    plt.title('Pattern Markers Plot')
    plt.savefig("{}_patternMarkers.png".format(fn))
    plt.show()
    return hm


def plotUMAP(result, genes_in_rows=True, fn=""):
    """ Create a UMAP plot

    Args:
        result (anndata or CogapsResult): An anndata object of result or CogapsResult object
        genes_in_rows (bool, optional): Scanpy needs genes in columns, cells in rows. Defaults to True.
    """    
    # result=result["anndata"]
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
    sc.pl.highest_expr_genes(result, n_top=20, save="{}_highestExpressedGenes.png".format(fn))
    sc.pp.filter_cells(result, min_genes=200)
    sc.pp.filter_genes(result, min_cells=3)
    sc.pp.log1p(result)
    sc.pp.highly_variable_genes(result, min_mean=0.0125, max_mean=3, min_disp=0.5)
    result = result[:, result.var.highly_variable]
    sc.pp.scale(result, max_value=10)
    sc.tl.pca(result, svd_solver='arpack')
    sc.pp.neighbors(result)
    sc.tl.umap(result)
    sc.pl.umap(result, color=patterns, save="{}_UMAP.png".format(fn))

if __name__ == '__main__':
    import pickle
    import sys
    import os
    import matplotlib as mpl
    mpl.use('tkagg')

    # path to your result file, from command line
    pkl_path = sys.argv[1]
    
    # this unpickles the result object for use
    result = pickle.load(open(pkl_path, "rb"))

    # get filename, to name the saved plots
    filename = os.path.basename(pkl_path).split('.')[0]

    # call some of the plotting functions and save
    plot(result, fn=filename)
    binaryA(result, threshold=2, fn=filename)
    plotPatternMarkers(result, fn=filename)
    plotResiduals(result, fn=filename)
    plotUMAP(result, fn=filename)
