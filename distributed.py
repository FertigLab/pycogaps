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
    data = helper_functions.toAnndata(path, hdf_key=params.coparams["hdfKey"])
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
                paramlst.append([data, params, i, sets[i], None])

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
            paramlst.append([data, params, i, sets[i], None])
        finalresult = pool.map(callInternalCoGAPS, paramlst)
        pool.close()
        pool.join()

    stitched = stitchTogether(finalresult, params, sets)
    print("15\n")
    finalresult = finalresult[0]
    print("AMEAN", stitched["Amean"])
    finalresult["GapsResult"].Amean = pycogaps.Matrix(stitched["Amean"])
    finalresult["GapsResult"].Asd = pycogaps.Matrix(stitched["Asd"])
    finalresult["GapsResult"].Pmean = pycogaps.Matrix(stitched["Pmean"])
    finalresult["GapsResult"].Psd = pycogaps.Matrix(stitched["Psd"])
    print("16\n")
    return finalresult


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
        print("gene subset", genesubset)
        params.coparams['geneNames'] = genesubset
        adata = adata[genesubset, :]
        params.coparams['subsetDim'] = 1
    else:
        samples = np.array(params.coparams['sampleNames'])
        print("samples:", samples)
        samplesubset = np.take(samples, subsetIndices)
        print("sample subset:", samplesubset)
        params.coparams['sampleNames'] = samplesubset
        adata = adata[:, samplesubset]
        params.coparams['subsetDim'] = 2
    print("AFTER SUBSETTING:", adata, adata.obs_names, adata.var_names)

    params.coparams['subsetIndices'] = subsetIndices
    params.gaps.workerID = workerID
    params.gaps.asynchronousUpdates = False
    params.gaps.maxThreads = 1
    print("ABOUT TO CALL INTERNAL COGAPS\n")
    gapsresult = PyCoGAPS.CoGAPS(adata, params, uncertainty, transposeData=params.coparams["transposeData"])

    return gapsresult


def findConsensusMatrix(unmatched, params):
    print("IN FIND CONSENSUS MATRIX\n")
    print("UNMATCHED:", unmatched)
    allpatterns = pd.DataFrame(np.hstack(unmatched))
    print("ALLPATTERNS:", allpatterns)
    comb = expandgrid(range(params.coparams["nSets"]), range(params.gaps.nPatterns))
    comb = list(comb.values())
    comb = pd.DataFrame(comb)
    comb = pd.DataFrame.transpose(comb)
    comb = comb.to_numpy()
    names = []
    for i in range(comb.shape[0]):
        names.append(str(comb[i, 0] + 1) + "." + str(comb[i, 1] + 1))
    allpatterns.columns = names
    print(allpatterns)
    return patternMatch(allpatterns, params)


def expandgrid(*itrs):
    print("IN EXPAND GRID\n")
    product = list(itertools.product(*itrs))
    return {'Var{}'.format(i + 1): [x[i] for x in product] for i in range(len(itrs))}


def patternMatch(allpatterns, params):
    print("IN PATTERN MATCH\n")
    clusters = corcut(allpatterns, params.coparams["cut"], params.coparams["minNS"])
    print("CLUSTERS:", clusters)
    def splitcluster(list, index, minNS):
        # list = pd.DataFrame(list)
        print("IN SPLIT CLUSTER*******************", list, index, minNS)
        split = corcut(index, 2, minNS)
        list[0] = split[0]
        if len(split) > 1:
            list.append(split[1])
        return list

    def toolarge(x):
        return x.shape[1] > params.coparams["maxNS"]

    print("CLUSTERS:", clusters)
    indx = [c for c in clusters if toolarge(c)]

    while len(indx) > 0:
        clusters = splitcluster(clusters, indx[0], params.coparams["minNS"])
        indx = [c for c in clusters if toolarge(c)]

    # create matrix of mean patterns - weighted by correlation to mean pattern
    meanpatterns = []
    for cluster in clusters:
        # cluster = pd.DataFrame(cluster)
        cr = corrToMeanPattern(cluster)
        print("CR:", cr)
        print("CLUSTER:", cluster)
        meanpatterns.append(np.average(cluster, weights=[cr[0] ** 3] * cluster.shape[0], axis=0))
    meanpatterns = pd.DataFrame(meanpatterns)
    # returned patterns after scaling max to 1
    result = {
        "clustered": clusters,
        "consensus": meanpatterns.apply(lambda x: x / x.max())
    }
    return result


def corrToMeanPattern(cluster):
    print("IN CORR TO MEAN PATTERN")
    print("CLUSTER:", cluster)
    meanpat = np.mean(cluster, axis=1)
    for i in range(cluster.shape[1]):
        x = pd.DataFrame(cluster.iloc[:, i])
        xc = x.corrwith(other=meanpat, drop=True, axis=0, )
        rx = round(xc, ndigits=3)
    return rx


def corcut(allpatterns, cut, minNS):
    print("IN CORCUT\n")
    dist = allpatterns.corr()
    dist = 1 - dist
    if dist.isnull().values.any():
        warnings.warn("NaN values in correlation of patterns... Aborting")
        return
    clusters = AgglomerativeClustering(affinity="precomputed", linkage="average", n_clusters=cut).fit(dist)
    print("clustering result", clusters.labels_)
    clustid = [None] * cut
    for id in set(clusters.labels_):
        if np.count_nonzero(clusters.labels_ == id) >= minNS:
            indices = [a for a, x in enumerate(clusters.labels_) if x == id]
            print("indices", indices)
            thislist = allpatterns.iloc[:, indices]
            print(thislist)
            clustid[id] = thislist
            # clustid.append(np.array(flatlist))
        else:
            warnings.warn("cluster did not meet minNS threshold and will be excluded")
    return clustid


def stitchTogether(result, params, sets):
    """
    concatenate final results across subsets
    @param result: list of CogapsResult objects
    @param params: CoParams object (params used to generate these results)
    @param sets: sets used to break apart data
    @return final GapsResult object
    """
    print("IN STITCH TOGETHER\n")
    setindices = [item for sublist in sets for item in sublist]
    # step 1:
    print("1\n")
    if params.coparams["distributed"] == "genome-wide":
        print("2\n")
        # combine A matrices
        Amean = pd.DataFrame()
        Asd = pd.DataFrame()
        for r in result:
            df1 = pd.DataFrame(np.array(r["GapsResult"].Amean))
            df2 = pd.DataFrame(np.array(r["GapsResult"].Asd))
            Amean = pd.concat([df1, Amean], axis=1)
            Asd = pd.concat([df2, Asd], axis=1)
        print("4\n")
        # # Amean = [item for sublist in Amean for item in sublist]
        # Amean = np.array(Amean)
        # # Asd = [item for sublist in Asd for item in sublist]
        # Asd = np.array(Asd)
        print("5\n")
        # copy P matrix.. same for all sets
        Pmean = pd.DataFrame(np.array(r["GapsResult"].Pmean))
        print("5.5 PMEAN:\n", Pmean)
        Psd = pd.DataFrame(np.array(r["GapsResult"].Psd))

        # reorder features to match data
        if len(Amean) == len(setindices):
            print("7\n")
            indices = range(Amean.shape[0])
            if indices.sort() == setindices.sort():
                print("8\n")
                reorder = set(indices) & set(setindices)
                Amean = [Amean[i] for i in reorder]
                Asd = [Asd[i] for i in reorder]
    else:
        print("9\n")
        Pmean = pd.DataFrame()
        Psd = pd.DataFrame()
        for r in result:
            df1 = pd.DataFrame(np.array(r["GapsResult"].Pmean))
            df2 = pd.DataFrame(np.array(r["GapsResult"].Psd))
            Pmean = pd.concat([df1, Pmean], axis=1)
            Psd = pd.concat([df2, Psd], axis=1)
        # Pmean = [item for sublist in Pmean for item in sublist]
        # Pmean = np.array(Pmean)
        # Psd = [item for sublist in Psd for item in sublist]
        # Psd = np.array(Psd)
        # copy P matrix.. same for all sets
        Amean = np.array(result[0]["GapsResult"].Amean)
        # Asd = pycogaps.Matrix(Amean.shape[0], Amean.shape[1])
        print("11\n")
        if Pmean.shape[0] == len(setindices):
            print("12\n")
            indices = range(Pmean.shape[0])
            if indices.sort() == setindices.sort():
                print("13\n")
                reorder = set(indices) & set(setindices)
                Pmean = [Pmean[i] for i in reorder]
                Psd = [Psd[i] for i in reorder]
    print("13.5\n")
    reslist = {
        "Amean": Amean,
        "Asd": Asd,
        "Pmean": Pmean,
        "Psd": Psd
    }
    print("14\n")
    return reslist
