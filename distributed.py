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
    if isinstance(path, str):
        data = helper_functions.toAnndata(path, hdf_counts_key=params.coparams["hdfKey"], hdf_dim1_key=params.coparams["hdfRowKey"], hdf_dim2_key=params.coparams["hdfColKey"], transposeData=params.coparams["transposeData"])
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
                paramlst.append([data, params, i, sets[i], None])

            result = pool.imap(callInternalCoGAPS, paramlst)
            pool.close()
            pool.join()
            result=list(result)

            if params.coparams['distributed'] == "genome-wide":
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
        finalresult = pool.imap(callInternalCoGAPS, paramlst)
        pool.close()
        pool.join()
        finalresult=list(finalresult)

    stitched = stitchTogether(finalresult, params, sets)
    finalresult = finalresult[0]
    if params.coparams["distributed"] == "genome-wide":
        print("AMEAN",np.array(result[0]["GapsResult"].Amean))
        finalresult["GapsResult"].Amean = pycogaps.Matrix(np.array(result[0]["GapsResult"].Amean))
        finalresult["GapsResult"].Asd = pycogaps.Matrix(np.array(result[0]["GapsResult"].Asd))
        finalresult["anndata"].obs = np.array(result[0]["GapsResult"].Amean)
        finalresult["anndata"].uns["asd"] = np.array(result[0]["GapsResult"].Asd)
        finalresult["GapsResult"].Pmean = pycogaps.Matrix(np.array(stitched["Pmean"]))
        finalresult["GapsResult"].Psd = pycogaps.Matrix(np.array(stitched["Psd"]))
    else:
        finalresult["GapsResult"].Amean = pycogaps.Matrix(np.array(stitched["Amean"]))
        finalresult["GapsResult"].Asd = pycogaps.Matrix(np.array(stitched["Asd"]))
        finalresult["anndata"].var = np.array(result[0]["GapsResult"].Pmean)
        finalresult["anndata"].uns["psd"] = np.array(result[0]["GapsResult"].Psd)
        finalresult["GapsResult"].Pmean = pycogaps.Matrix(np.array(result[0]["GapsResult"].Pmean))
        finalresult["GapsResult"].Psd = pycogaps.Matrix(np.array(result[0]["GapsResult"].Psd))
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
    print("COMB", comb)
    for i in range(comb.shape[0]):
        names.append(str(comb[i, 0] + 1) + "." + str(comb[i, 1] + 1))
    allpatterns.columns = names
    print("NAMES", names)
    return patternMatch(allpatterns, params)


def expandgrid(*itrs):
    product = list(itertools.product(*itrs))
    return {'Var{}'.format(i + 1): [x[i] for x in product] for i in range(len(itrs))}


def patternMatch(allpatterns, params):
    clusters = corcut(allpatterns, params.coparams["cut"], params.coparams["minNS"])
    maxNS = params.coparams["maxNS"]
    print("MAXNS", maxNS)
    def splitcluster(list, index, minNS):
        #try:
        idx = list.index(index)
        #except ValueError:
         #   print("cluster not found...")
          #  return None
        split = corcut(list[idx], 2, minNS)
        print("Length of split", len(split))
        list[idx] = split[0]
        if len(split) > 1:
            list.append(split[1])
        return list

    def toolarge(x):
        if x is None:
            return False
        return x.shape[1] > maxNS

    indx = [c for c in clusters if toolarge(c)]

    while len(indx) > 0:
        clusters = splitcluster(clusters, indx[0], params.coparams["minNS"])
        #if clusters is None:
        #    break
        #clusters = [c for c in clusters if c is not None]
        indx = [c for c in clusters if toolarge(c)]
        print("SHAPE OF INDX:", len(indx))
        # # print("INDX:", indx)
        # print("len(clusters)", len(clusters))

    # create matrix of mean patterns - weighted by correlation to mean pattern
    meanpatterns = []
    for cluster in clusters:
        if cluster is not None:
            cr = corrToMeanPattern(cluster)
            meanpatterns.append(np.average(cluster, weights=[cr[0] ** 3] * cluster.shape[0], axis=0))
    meanpatterns = pd.DataFrame(meanpatterns)
    # returned patterns after scaling max to 1
    result = {
        "clustered": clusters,
        "consensus": meanpatterns.apply(lambda x: x / x.max())
    }
    return result


def corrToMeanPattern(cluster):
    # print("cluster:", cluster)
    meanpat = np.mean(cluster, axis=1)
    for i in range(cluster.shape[1]):
        x = pd.DataFrame(cluster.iloc[:, i])
        xc = x.corrwith(other=meanpat, drop=True, axis=0, )
        rx = round(xc, ndigits=3)
    return rx


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


def stitchTogether(result, params, sets):
    """
    concatenate final results across subsets
    @param result: list of CogapsResult objects
    @param params: CoParams object (params used to generate these results)
    @param sets: sets used to break apart data
    @return final GapsResult object
    """
    print("Stitching results together...")
    if params.coparams["distributed"] == "genome-wide":
        # combine A matrices
        Amean = pd.DataFrame()
        Asd = pd.DataFrame()
        for r in result:
            df1 = r["anndata"].obs
            df2 = r["anndata"].uns["asd"]
            Amean = pd.concat([df1, Amean])
            Asd = pd.concat([df2, Asd])
        # copy P matrix.. same for all sets
        Pmean = result[0]["anndata"].var
        Psd = result[0]["anndata"].uns["psd"]
    else:
        Pmean = pd.DataFrame()
        Psd = pd.DataFrame()
        for r in result:
            df1 = pd.DataFrame(np.array(r["GapsResult"].Pmean))
            df2 = pd.DataFrame(np.array(r["GapsResult"].Psd))
            Pmean = pd.concat([df1, Pmean], axis=1)
            Psd = pd.concat([df2, Psd], axis=1)
        Amean = np.array(result[0]["GapsResult"].Amean)
        Asd = np.array(result[0]["GapsResult"].Asd)
    reslist = {
        "Amean": Amean,
        "Asd": Asd,
        "Pmean": Pmean,
        "Psd": Psd
    }
    return reslist
