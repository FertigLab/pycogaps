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
                paramlst.append([data, params, i, sets[i], None])

            result = pool.imap(callInternalCoGAPS, paramlst)
            pool.close()
            pool.join()
            result = list(result)

            if params.coparams['distributed'] == "genome-wide":
                # unmatched = []
                # for r in result:
                #     unmatched.append(np.array(r["GapsResult"].Pmean))
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
            paramlst.append([data, params, i, sets[i], None])
        finalresult = pool.imap(callInternalCoGAPS, paramlst)
        pool.close()
        pool.join()
        finalresult = list(finalresult)

    stitched = stitchTogether(finalresult, result, params, sets)
    gapsresult = pycogaps.GapsResult
    if params.coparams["distributed"] == "genome-wide":
        print('~~~~~~~ In distributed genome-wide ~~~~~~~~~~')
        '''
        gapsresult.Amean = stitched["Amean"]
        gapsresult.Asd = stitched["Asd"]
        gapsresult.Pmean = stitched["Pmean"]
        gapsresult.Psd = stitched["Psd"]
        # print("AMEAN",np.array(result[0]["GapsResult"].Amean))
        # finalresult["GapsResult"].Amean = pycogaps.Matrix(np.array(result[0]["GapsResult"].Amean))
        # finalresult["GapsResult"].Asd = pycogaps.Matrix(np.array(result[0]["GapsResult"].Asd))
        # finalresult["anndata"].obs = np.array(result[0]["GapsResult"].Amean)
        # finalresult["anndata"].uns["asd"] = np.array(result[0]["GapsResult"].Asd)
        # finalresult["GapsResult"].Pmean = pycogaps.Matrix(np.array(stitched["Pmean"]))
        # finalresult["GapsResult"].Psd = pycogaps.Matrix(np.array(stitched["Psd"]))
        '''
        # debug anndata here
        '''
        print('obs names: ', finalresult[0]['anndata'].obs_names)
        print('var names: ', finalresult[0]['anndata'].var_names)
        print('obs: ', finalresult[0]['anndata'].obs)
        print('var: ', finalresult[0]['anndata'].var)
        print('Amean: ', gapsresult.Amean)
        print('Amean shape: ', gapsresult.Amean.shape)
        print('Pmean: ', gapsresult.Pmean)
        print('Pmean shape: ', gapsresult.Pmean.shape)
        finalresult[0]['anndata'].obs = gapsresult.Amean
        finalresult[0]['anndata'].var = gapsresult.Pmean
        finalresult[0]['anndata'].uns['asd'] = gapsresult.Asd
        finalresult[0]['anndata'].uns['psd'] = gapsresult.Psd
        '''
        gapsresult.Amean = pycogaps.Matrix(np.array(stitched["Amean"]))
        gapsresult.Asd = pycogaps.Matrix(np.array(stitched["Asd"]))
        gapsresult.Pmean = pycogaps.Matrix(np.array(stitched["Pmean"]))
        gapsresult.Psd = pycogaps.Matrix(np.array(stitched["Psd"]))

        # finalresult["GapsResult"].Amean = pycogaps.Matrix(np.array(stitched["Amean"]))
        # finalresult["GapsResult"].Asd = pycogaps.Matrix(np.array(stitched["Asd"]))
        # finalresult["anndata"].var = np.array(result[0]["GapsResult"].Pmean)
        # finalresult["anndata"].uns["psd"] = np.array(result[0]["GapsResult"].Psd)
        # finalresult["GapsResult"].Pmean = pycogaps.Matrix(np.array(result[0]["GapsResult"].Pmean))
        # finalresult["GapsResult"].Psd = pycogaps.Matrix(np.array(result[0]["GapsResult"].Psd))

        # finalresult[0]['anndata'].obs = 



    else:
        # gapsresult.Amean = stitched["Amean"]
        # gapsresult.Asd = stitched["Asd"]
        # gapsresult.Pmean = stitched["Pmean"]
        # gapsresult.Psd = stitched["Psd"]
        
        # finalresult["GapsResult"].Amean = pycogaps.Matrix(np.array(stitched["Amean"]))
        # finalresult["GapsResult"].Asd = pycogaps.Matrix(np.array(stitched["Asd"]))
        # finalresult["anndata"].var = np.array(result[0]["GapsResult"].Pmean)
        # finalresult["anndata"].uns["psd"] = np.array(result[0]["GapsResult"].Psd)
        # finalresult["GapsResult"].Pmean = pycogaps.Matrix(np.array(result[0]["GapsResult"].Pmean))
        # finalresult["GapsResult"].Psd = pycogaps.Matrix(np.array(result[0]["GapsResult"].Psd))

        gapsresult.Amean = pycogaps.Matrix(np.array(stitched["Amean"]))
        gapsresult.Asd = pycogaps.Matrix(np.array(stitched["Asd"]))
        gapsresult.Pmean = pycogaps.Matrix(np.array(stitched["Pmean"]))
        gapsresult.Psd = pycogaps.Matrix(np.array(stitched["Psd"]))

    return {
        "GapsResult": gapsresult,
        # "anndata": finalresult[0]["anndata"]
        # 'anndata': finalresult[0]['anndata']
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

    def splitcluster(allpatterns, index, minNS):
        print("LIST", allpatterns)
        print("INDEX", index)
        for i in np.arange(len(allpatterns)):
            if len(allpatterns[i].columns.intersection(index.columns)) == allpatterns[i].shape[1]:
                    subindex= i
        split = corcut(allpatterns[subindex], 2, minNS)
        print("Length of split", len(split))
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
        # if clusters is None:
        #    break
        # clusters = [c for c in clusters if c is not None]
        indx = [c for c in clusters if toolarge(c)]
        print("SHAPE OF INDX:", len(indx))
        # # print("INDX:", indx)
        # print("len(clusters)", len(clusters))
    print("AFTER SPlITTING--CLUSTERS\n", clusters)
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
    print("MEAN PATTERNS\n", meanpatterns)
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
    # print('Amean result shape: ', np.array(finalresult[0]["anndata"].obs).shape)
    # print('Pmean result shape: ', np.array(finalresult[0]["GapsResult"].Pmean).shape)
    # print('Amean nonzero: ', np.count_nonzero(np.array(finalresult[0]["anndata"].obs)))
    # print('Pmean nonzero: ', np.count_nonzero(np.array(finalresult[0]["GapsResult"].Pmean)))

    print("Stitching results together...")
    if params.coparams["distributed"] == "genome-wide":
        '''
        # combine A matrices
        df1_ann = []
        df2_ann = []
        df1 = []
        df2 = []
        '''
        # debug
        '''
        for r in result:
            print('shape of r Amean: ', r['anndata'].obs.shape)
            df1_ann.append(r["anndata"].obs)
            df2_ann.append(r['anndata'].uns['asd'])
            # print('shape of r GapsResult Amean: ', np.array(r['GapsResult'].Amean).shape)
            # df1.append(np.array(r['GapsResult'].Amean))
            # df2.append(np.array(r['GapsResult'].Asd))
        # combadata = anndata.concat(df1)
        Amean = pd.concat(df1_ann, axis=0)
        Asd = pd.concat(df2_ann, axis=0)
        # print('Amean: ', np.count_nonzero(df1))
        # print('Asd: ', np.count_nonzero(df2))
        print('Amean ann: ', np.count_nonzero(Amean))
        print('Asd ann: ', np.count_nonzero(Asd))
        # Amean = pd.DataFrame(np.array(result[0]["GapsResult"].Amean))
        # Asd = pd.DataFrame(np.array(result[0]["GapsResult"].Asd))
        setindices = [s for st in sets for s in st]
        '''
        # if (Amean.shape[0] == len(setindices)):
        #     indices = np.arange(len(setindices))
        #     if (indices.sort() == setindices.sort()):
        #         reorder = list(set(indices).intersection(setindices))
        #         print('reorder: ', reorder)
        #         print('len of reorder: ', len(reorder))
        #         print('duplicated Amean: ', Amean[Amean.index.duplicated()])
        #         Amean.reindex(reorder, axis=0)
        #         print('duplicated Amean: ', Amean[Amean.index.duplicated()])
        #         Asd.reindex(reorder)
        '''
        print('finished Amean stitch... starting Pmean stitch')
        # copy P matrix.. same for all sets
        Pmean = pd.DataFrame(np.array(finalresult[0]['GapsResult'].Pmean))
        Psd = pd.DataFrame(np.array(finalresult[0]["GapsResult"].Psd))
        print('Pmean: ', np.count_nonzero(Pmean))
        print('Psd: ', np.count_nonzero(Psd))
        # Pmean_ann = finalresult[0]['anndata'].var
        # Psd_ann = finalresult[0]['anndata'].uns['psd']
        # print('Pmean ann: ', np.count_nonzero(Pmean_ann))
        # print('Psd ann: ', np.count_nonzero(Psd_ann))
        # Amean = pycogaps.Matrix(np.array(Amean))
        # Asd = pycogaps.Matrix(np.array(Asd))


        '''

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