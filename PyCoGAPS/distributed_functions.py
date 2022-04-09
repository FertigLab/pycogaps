from PyCoGAPS.config import *
from PyCoGAPS.helper_functions import *
from PyCoGAPS.subset_data import *

import itertools
from sklearn.cluster import AgglomerativeClustering
from scipy.stats.stats import pearsonr

def findConsensusMatrix(unmatched, params):
    # print("FINDING CONSENSUS MATRIX")
    # allpatterns = pd.DataFrame(np.hstack(unmatched))
    allpatterns = pd.DataFrame(unmatched)
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
    # print("IN PATTERNMATCH")
    clusters = corcut(allpatterns, params.coparams["cut"], params.coparams["minNS"])
    maxNS = params.coparams["maxNS"]
    # print("MAXNS", maxNS)

    def splitcluster(allpatterns, index, minNS):
        # print("IN SPLIT CLUSTER")
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
    
    meanpatterns = np.empty((len(clusters[0]), len(clusters)))
    for i in range(len(clusters)):
        cluster = clusters[i]
        if cluster is not None:
            cr = np.array(corrToMeanPattern(cluster))**3
            meanpat = np.empty((cluster.shape[0]))
            for row in range(cluster.shape[0]):
                meanpat[row] = np.average(cluster.iloc[row], weights=cr)
            meanpatterns[:,i] = meanpat
    meanpatterns = np.divide(meanpatterns, np.max(meanpatterns, axis=0))
    meanpatterns = pd.DataFrame(data=meanpatterns)
    
    # returned patterns after scaling max to 1
    result = {
        "clustered": clusters,
        'consensus': meanpatterns
    }
    return result


def corrToMeanPattern(cluster):
    # print("IN CORR TO MEAN PATTERN")
    # print("cluster:", cluster)
    meanpat = cluster.mean(axis=1, skipna=True)
    corrmat = []
    for column in cluster:
        corrmat.append(pearsonr(cluster[column], meanpat)[0])
    return corrmat



def stitchTogether(finalresult, result, params, sets, adata):
    """
    concatenate final results across subsets
    @param result: list of CogapsResult objects
    @param params: CoParams object (params used to generate these results)
    @param sets: sets used to break apart data
    @return final GapsResult object
    """

    print("Stitching results together...")
    if params.coparams["distributed"] == "genome-wide":
        Amean = finalresult[0].obs
        Asd = finalresult[0].uns['asd']
        for r in finalresult[1:]:
            Amean = Amean.append(r.obs)
            Asd = Asd.append(r.uns["asd"])
        Amean = Amean.reindex(adata.obs_names)
        Asd = Asd.reindex(adata.obs_names)
        Pmean = finalresult[0].var.reindex(adata.var_names)
        Psd = finalresult[0].uns['psd'].reindex(adata.var_names)

    else:
        Pmean = np.array(finalresult[0].var)
        Psd = np.array(finalresult[0].uns['psd'])
        for r in finalresult[1:]:
            df1 = np.array(r.var)
            df2 = np.array(r.uns['psd'])
            Pmean = np.append(Pmean, df1, axis=0)
            Psd = np.append(Psd, df2, axis=0)
        Pmean = np.array(finalresult[0].obs)
        Psd = np.array(finalresult[0].uns['psd'])
        
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

