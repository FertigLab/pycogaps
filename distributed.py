import warnings
from scipy.stats.stats import pearsonr
import pycogaps
from pandas.io.clipboard import paste

import PyCoGAPS
import subset_data
import multiprocessing
import helper_functions
import numpy as np
import pandas as pd
import itertools
from sklearn.cluster import AgglomerativeClustering

def distributedCoGAPS(path, params, uncertainty=None):
    if isinstance(path, str):
        data = helper_functions.toAnndata(path, hdf_counts_key=params.coparams["hdfKey"],
                                          hdf_dim1_key=params.coparams["hdfRowKey"],
                                          hdf_dim2_key=params.coparams["hdfColKey"],
                                          transposeData=params.coparams["transposeData"])
    else:
        data = path
    sets = subset_data.createSets(data, params)
    if min(map(len, sets)) < params.gaps.nPatterns:
        warnings.warn("Data subset dimension less than nPatterns. Aborting.")
        return 1

    PyCoGAPS.setParams(params, {'checkpointOutFile': ""})

    if params.coparams["fixedPatterns"] is None:
        print("Running Across Subsets...\n\n")
        with multiprocessing.get_context("spawn").Pool(processes=len(sets)) as pool:

            # make a list of parameters for each function call so they can easily be mapped to processes
            paramlst = []
            for i in range(len(sets)):
                paramlst.append([data, params, i, sets[i], uncertainty])

            result = pool.imap(callInternalCoGAPS, paramlst)
            pool.close()
            pool.join()
            result = list(result)

            if params.coparams['distributed'] == "genome-wide":
                unmatched = map(lambda x: np.array(x["GapsResult"].Pmean), result)
            else:
                unmatched = map(lambda x: np.array(x["GapsResult"].Amean), result)
            print("Matching patterns across subsets...\n")
            matched = findConsensusMatrix(unmatched, params)
    else:
        matched = params.gaps.fixedPatterns

    params.gaps.nPatterns = matched["consensus"].shape[1]
    params.gaps.fixedPatterns = pycogaps.Matrix(matched["consensus"])
    print("FIXED PATTERNS\n", matched["consensus"])
    if params.coparams["distributed"] == "genome-wide":
        params.gaps.whichMatrixFixed = "P"
    else:
        params.gaps.whichMatrixFixed = "A"

    print("Running final stage...")
    with multiprocessing.get_context("spawn").Pool(processes=len(sets)) as pool:
        paramlst = []
        for i in range(len(sets)):
            paramlst.append([data, params, i, sets[i], uncertainty])
        finalresult = pool.imap(callInternalCoGAPS, paramlst)
        pool.close()
        pool.join()
        finalresult = list(finalresult)

    stitched = stitchTogether(finalresult, result, params, sets)
    gapsresult = pycogaps.GapsResult
    if params.coparams["distributed"] == "genome-wide":
        gapsresult.Amean = pycogaps.Matrix(np.array(stitched["Amean"]))
        gapsresult.Asd = pycogaps.Matrix(np.array(stitched["Asd"]))
        gapsresult.Pmean = pycogaps.Matrix(np.array(stitched["Pmean"]))
        gapsresult.Psd = pycogaps.Matrix(np.array(stitched["Psd"]))

    else:
        gapsresult.Amean = pycogaps.Matrix(np.array(stitched["Amean"]))
        gapsresult.Asd = pycogaps.Matrix(np.array(stitched["Asd"]))
        gapsresult.Pmean = pycogaps.Matrix(np.array(stitched["Pmean"]))
        gapsresult.Psd = pycogaps.Matrix(np.array(stitched["Psd"]))

    return {
        "GapsResult": gapsresult,
        "anndata": PyCoGAPS.GapsResultToAnnData(gapsresult, data, params)
    }


def callInternalCoGAPS(paramlst):
    # take out parameters passed as a list to the worker process
    path = paramlst[0]
    params = paramlst[1]
    workerID = paramlst[2]
    subsetIndices = paramlst[3]
    uncertainty = paramlst[4]
    if isinstance(path, str):
        adata = helper_functions.toAnndata(path)
    else:
        adata = path
    if subsetIndices is None:
        print("No subset indices provided; generating random sets...")
        subsetIndices = subset_data.createSets(adata, params)
    if params.coparams['distributed'] == "genome-wide":
        genes = np.array(params.coparams['geneNames'])
        genesubset = np.take(genes, subsetIndices)
        params.coparams['geneNames'] = set(genesubset)
        adata = adata[subsetIndices, :]
        params.coparams['subsetDim'] = 1
    else:
        samples = np.array(params.coparams['sampleNames'])
        samplesubset = np.take(samples, subsetIndices)
        params.coparams['sampleNames'] = samplesubset
        adata = adata[samplesubset, :]
        params.coparams['subsetDim'] = 2

    params.coparams['subsetIndices'] = subsetIndices
    params.gaps.workerID = workerID
    params.gaps.asynchronousUpdates = False
    params.gaps.maxThreads = 1
    print("Calling internal CoGAPS...\n")
    gapsresult = PyCoGAPS.CoGAPS(adata, params, uncertainty, transposeData=params.coparams["transposeData"])

    return gapsresult


def findConsensusMatrix(unmatched, params):
    allpatterns = pd.DataFrame(np.hstack(unmatched))
    comb = expandgrid(range(params.coparams["nSets"]), range(params.gaps.nPatterns))
    comb = list(comb.values())
    comb = pd.DataFrame(comb)
    comb = pd.DataFrame.transpose(comb)
    comb = comb.to_numpy()
    names = []
    # print("COMB", comb)
    for i in range(comb.shape[0]):
        names.append(str(comb[i, 0] + 1) + "." + str(comb[i, 1] + 1))
    allpatterns.columns = names
    # print("NAMES", names)
    return patternMatch(allpatterns, params)


def expandgrid(*itrs):
    product = list(itertools.product(*itrs))
    return {'Var{}'.format(i + 1): [x[i] for x in product] for i in range(len(itrs))}


def patternMatch(allpatterns, params):
    clusters = corcut(allpatterns, params.coparams["cut"], params.coparams["minNS"])
    maxNS = params.coparams["maxNS"]
    # print("MAXNS", maxNS)

    def splitcluster(allpatterns, index, minNS):
        # print("LIST", allpatterns)
        # print("INDEX", index)
        for i in np.arange(len(allpatterns)):
            if len(allpatterns[i].columns.intersection(index.columns)) == allpatterns[i].shape[1]:
                    subindex= i
        split = corcut(allpatterns[subindex], 2, minNS)
        # print("Length of split", len(split))
        allpatterns[subindex] = split[0]
        if len(split) > 1:
            allpatterns.append(split[1])
        return allpatterns

    def toolarge(x):
        if x is None:
            return False
        return x.shape[1] > maxNS

    indx = [c for c in clusters if toolarge(c)]

    while len(indx) > 0:
        clusters = splitcluster(clusters, indx[0], params.coparams["minNS"])
        indx = [c for c in clusters if toolarge(c)]
        # print("SHAPE OF INDX:", len(indx))

    # print("AFTER SPlITTING--CLUSTERS\n", clusters)
    # create matrix of mean patterns - weighted by correlation to mean pattern
    meanpatterns = []
    for i in np.arange(len(clusters)):
        cluster = clusters[i]
        if cluster is not None:
            cr = corrToMeanPattern(cluster)
            meanpat = []
            for row in np.arange(cluster.shape[0]):
                avg = np.average(cluster.iloc[row], weights=cr)
                print(avg)
                meanpat.append(avg)
            meanpatterns.append(meanpat)
    meanpatterns = pd.DataFrame(meanpatterns)
    meanpatterns = meanpatterns.transpose()
    # print("MEAN PATTERNS\n", meanpatterns)
    # returned patterns after scaling max to 1
    result = {
        "clustered": clusters,
        "consensus": meanpatterns.apply(lambda x: x / x.max())
    }
    return result


def corrToMeanPattern(cluster):

    # print("cluster:", cluster)
    cluster = cluster.dropna()
    meanpat = cluster.mean(axis=1)
    corrmat = []
    for column in cluster:
        corrmat.append(pearsonr(cluster[column], meanpat)[0])
    return corrmat



def stitchTogether(finalresult, result, params, sets):
    """
    concatenate final results across subsets
    @param result: list of CogapsResult objects
    @param params: CoParams object (params used to generate these results)
    @param sets: sets used to break apart data
    @return final GapsResult object
    """

    print("Stitching results together...")
    if params.coparams["distributed"] == "genome-wide":

        Amean = pd.DataFrame()
        Asd = pd.DataFrame()
        for r in finalresult:
            df1 = pd.DataFrame(np.array(r["anndata"].obs))
            df2 = pd.DataFrame(np.array(r["anndata"].uns['asd']))
            Amean = pd.concat([df1, Amean], axis=0)
            Asd = pd.concat([df2, Asd], axis=0)
        Pmean = np.array(finalresult[0]["GapsResult"].Pmean)
        Psd = np.array(finalresult[0]["GapsResult"].Psd)


    else:
        Pmean = pd.DataFrame()
        Psd = pd.DataFrame()
        for r in finalresult:
            df1 = pd.DataFrame(np.array(r["GapsResult"].Pmean))
            df2 = pd.DataFrame(np.array(r["GapsResult"].Psd))
            Pmean = pd.concat([df1, Pmean], axis=1)
            Psd = pd.concat([df2, Psd], axis=1)
        Amean = np.array(finalresult[0]["GapsResult"].Amean)
        Asd = np.array(finalresult[0]["GapsResult"].Asd)
        
    reslist = {
        "Amean": Amean,
        "Asd": Asd,
        "Pmean": Pmean,
        "Psd": Psd
    }
    return reslist


def corcut(allpatterns, cut, minNS):
    dist = allpatterns.corr()
    dist = 1 - dist
    if dist.isnull().values.any():
        warnings.warn("NaN values in correlation of patterns... Aborting")
        return
    clusters = AgglomerativeClustering(affinity="precomputed", linkage="average", n_clusters=cut).fit(dist)
    clustid = []
    for id in set(clusters.labels_):
        if np.count_nonzero(clusters.labels_ == id) >= minNS:
            indices = [a for a, x in enumerate(clusters.labels_) if x == id]
            thislist = allpatterns.iloc[:, indices]
            clustid.append(thislist)
        else:
            warnings.warn("cluster did not meet minNS threshold and will be excluded")
    # print("CORCUT cluster IDs:", clustid)
    return clustid