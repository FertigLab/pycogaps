import warnings

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
    data = helper_functions.toAnndata(path)
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
                paramlst.append([path, params, i, sets[i], None])

            result = pool.map(callInternalCoGAPS, paramlst)
            pool.close()
            pool.join()

            if params.coparams['distributed'] == "genome-wide":
                # unmatched = list(map(lambda x: np.array(x["GapsResult"].Pmean), result))
                unmatched = []
                for r in result:
                    unmatched.append(np.array(r["GapsResult"].Pmean))
            else:
                unmatched = map(lambda x: np.array(x["GapsResult"].Amean), result)
            print("Matching patterns across subsets...\n")
            matched = findConsensusMatrix(unmatched, params)
    else:
        matched = params.gaps.fixedPatterns

    params.gaps.nPatterns = matched["consensus"].shape[1]
    params.gaps.fixedPatterns = pycogaps.Matrix(matched["consensus"])
    if params.coparams["distributed"] == "genome-wide":
        params.gaps.whichMatrixFixed = "P"
    else:
        params.gaps.whichMatrixFixed = "A"

    print("Running final stage...")
    with multiprocessing.get_context("spawn").Pool(processes=len(sets)) as pool:
        paramlst = []
        for i in range(len(sets)):
            paramlst.append([path, params, i, sets[i], None])
        finalresult = pool.map(callInternalCoGAPS, paramlst)
        pool.close()
        pool.join()

    fullresult = stitchTogether(finalresult, params, sets)

    # add diagnostics...
    fullresult.diagnostics.firstPass = result

    return fullresult


def callInternalCoGAPS(paramlst):
    # take out parameters passed as a list to the worker process
    path = paramlst[0]
    params = paramlst[1]
    workerID = paramlst[2]
    subsetIndices = paramlst[3]
    uncertainty = paramlst[4]
    if subsetIndices is None:
        print("No subset indices provided; generating random sets...")
        adata = helper_functions.toAnndata(path)
        subsetIndices = subset_data.createSets(adata, params)
    if params.coparams['distributed'] == "genome-wide":
        genes = np.array(params.coparams['geneNames'])
        params.coparams['geneNames'] = np.take(genes, subsetIndices)
        params.coparams['subsetDim'] = 1
    else:
        samples = np.array(params.coparams['sampleNames'])
        params.coparams['sampleNames'] = np.take(samples, subsetIndices)
        params.coparams['subsetDim'] = 2

    params.coparams['subsetIndices'] = subsetIndices
    params.gaps.workerID = workerID
    params.gaps.asynchronousUpdates = False
    params.gaps.maxThreads = 1
    gapsresult = PyCoGAPS.CoGAPS(path, params, uncertainty)

    return gapsresult


def findConsensusMatrix(unmatched, params):
    allpatterns = pd.DataFrame(np.hstack(unmatched))
    comb = expandgrid(range(params.coparams["nSets"]), range(params.gaps.nPatterns))
    comb = list(comb.values())
    comb = pd.DataFrame(comb)
    comb = pd.DataFrame.transpose(comb)
    comb = comb.to_numpy()
    names = []
    for i in range(comb.shape[0]):
        names.append(str(comb[i, 0]+1) + "." + str(comb[i, 1]+1))
    allpatterns.columns = names
    return patternMatch(allpatterns, params)


def expandgrid(*itrs):
   product = list(itertools.product(*itrs))
   return {'Var{}'.format(i+1):[x[i] for x in product] for i in range(len(itrs))}


def patternMatch(allpatterns, params):
    clusters = corcut(allpatterns, params.coparams["cut"], params.coparams["minNS"])

    def splitcluster(list, index, minNS):
        split = corcut(list.iloc[:, index], 2, minNS)
        list[index] = split['0']
        if len(split) > 1:
            list = np.concatenate(list, split['1'])
        return list

    def toolarge(x):
        return x.shape[1] > params.coparams["maxNS"]

    indx = [c for c in clusters.values() if toolarge(c)]

    while(len(indx) > 0):
        clusters = splitcluster(clusters, indx[0], params.coparams["minNS"])
        indx = [c for c in clusters if toolarge(c)]

    # create matrix of mean patterns - weighted by correlation to mean pattern
    meanpatterns = []
    for cluster in clusters.values():
        cr = corrToMeanPattern(cluster)
        meanpatterns.append(np.average(cluster, weights=[cr[0]**3] * cluster.shape[0], axis=0))
    meanpatterns = pd.DataFrame(meanpatterns)
    print("meanpatterns\n", meanpatterns)
    # returned patterns after scaling max to 1
    result = {
        "clustered": clusters,
        "consensus": meanpatterns.apply(lambda x: x / x.max())
    }
    return result


def corrToMeanPattern(cluster):
    meanpat = np.mean(cluster, axis=1)
    for i in range(cluster.shape[1]):
        x = pd.DataFrame(cluster.iloc[:,i])
        xc = x.corrwith(other=meanpat, drop=True, axis=0, )
        rx = round(xc, ndigits=3)
    return rx


def corcut(allpatterns, cut, minNS):
    dist = allpatterns.corr()
    dist = 1-dist
    if dist.isnull().values.any():
        warnings.warn("NaN values in correlation of patterns... Aborting")
        return
    clusters = AgglomerativeClustering(affinity="precomputed", linkage="average", n_clusters=cut).fit(dist)
    clustid = dict()
    for i in range(cut):
        if np.count_nonzero(clusters.labels_ == clusters.labels_[i]) >= minNS:
            indices = [a for a, x in enumerate(clusters.labels_) if x == i]
            thislist = allpatterns.iloc[:, indices]
            clustid[i] = thislist
        else:
            warnings.warn("cluster did not meet minNS threshold and will be excluded")
    return clustid


def stitchTogether(result, params, sets):
    print("not yet implemented")