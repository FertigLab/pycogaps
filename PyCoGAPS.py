import time
import math

import pandas as pd
import pycogaps
import anndata
from helper_functions import *
from subset_data import *
import sys


from scipy.stats.stats import pearsonr

import subset_data
import multiprocessing
import helper_functions
import numpy as np
import itertools
from sklearn.cluster import AgglomerativeClustering

sys.setrecursionlimit(10000)

print("""\

______      _____       _____   ___  ______  _____ 
| ___ \    /  __ \     |  __ \ / _ \ | ___ \/  ___|
| |_/ /   _| /  \/ ___ | |  \// /_\ \| |_/ /\ `--. 
|  __/ | | | |    / _ \| | __ |  _  ||  __/  `--. |
| |  | |_| | \__/\ (_) | |_\ \| | | || |    /\__/ /
\_|   \__, |\____/\___/ \____/\_| |_/\_|    \____/ 
       __/ |                                       
      |___/             
                                 
                    """)


class CoParams:
    """ Encapsulates all parameters for PyCoGAPS.

    """

    def __init__(self, path=None, matrix=None, transposeData=False, hdfKey=None, hdfRowKey=None, hdfColKey=None):
        """ Initializes CoParams object. 
            self.gaps : GapsParameters object
            self.cogaps : dictionary of additional parameters (not in GapsParameters)

        Args:
            path (str, optional): Path to data. Defaults to None.
            matrix (anndata, optional): AnnData object containing supplied data matrix. Defaults to None.
            transposeData (bool, optional): Expects genes x samples. Defaults to False.
            hdfKey (str, optional): For reading .h5 files. Defaults to None.
            hdfRowKey (str, optional): For reading .h5 files. Defaults to None.
            hdfColKey (str, optional): For reading .h5 files. Defaults to None.

        Raises:
            Exception: If path or params not passed as an argument.
        """        

        if matrix is not None:
            self.gaps = GapsParameters(pycogaps.Matrix(matrix.X))
            adata = matrix
        elif path is not None:
            if path.lower().endswith(".h5"):
                adata = toAnndata(path, hdfKey, hdfRowKey, hdfColKey, transposeData=transposeData)
            elif path.lower().endswith(".txt"):
                table = pd.read_table(path)
                adata = anndata.AnnData(table.iloc[:, 2:])
                adata.obs_names = table["symbol"]
            else:
                adata = toAnndata(path, transposeData=transposeData)
            matrix = pycogaps.Matrix(adata.X)
            self.gaps = GapsParameters(matrix)
        else:
            raise Exception('initialize with path= or params=')

        self.coparams = {'cut': self.gaps.nPatterns,
                         'nSets': 4,
                         'minNS': None,
                         'maxNS': None,
                         'explicitSets': None,
                         'samplingAnnotation': None,
                         'samplingWeight': None,
                         'subsetIndices': None,
                         'subsetDim': 0,
                         'geneNames': adata.obs_names,
                         'sampleNames': adata.var_names,
                         'fixedPatterns': None,
                         'distributed': None,
                         'hdfKey': hdfKey,
                         'hdfRowKey': hdfRowKey,
                         'hdfColKey': hdfColKey,
                         'useSparseOptimization': None,
                         'transposeData': transposeData,
                         }
        self.coparams['minNS'] = math.ceil(self.coparams['cut'] / 2)
        self.coparams['maxNS'] = self.coparams['minNS'] + self.coparams['nSets']

    def setDistributedParams(self, nSets=None, cut=None, minNS=None, maxNS=None):
        """ Sets parameters for running distributed CoGAPS.

        Args:
            nSets (int, optional): Number of sets to break data into. Defaults to None.
            cut (int, optional): Number of branches at which to cut dendrogram used in pattern matching. Defaults to None.
            minNS (int, optional): [description]. Minimum of individual set contributions a cluster must contain. Defaults to None.
            maxNS (int, optional): [description]. Maximum of individual set contributions a cluster can contain. Defaults to None.
        """        

        print("setting distributed parameters - call this again if you change nPatterns")
        if self.coparams['distributed'] != "genome-wide":
            print("if you wish to perform genome-wide distributed cogaps, please run setParams(params, "
                  "\"distributed\", ""\"genome-wide\")")
        if nSets is None:
            self.coparams['nSets'] = self.coparams['nSets']
        else:
            self.coparams['nSets'] = nSets
        if cut is None:
            self.coparams['cut'] = self.gaps.nPatterns
        else:
            self.coparams['cut'] = cut
        if minNS is None:
            self.coparams['minNS'] = math.ceil(self.coparams['cut'] / 2)
        else:
            self.coparams['minNS'] = minNS
        if maxNS is None:
            self.coparams['maxNS'] = self.coparams['minNS'] + self.coparams['nSets']
        else:
            self.coparams['maxNS'] = maxNS

    # samplingWeight is a dictionary
    # can use: dict(zip(names, weights))
    def setAnnotationWeights(self, annotation, weight):
        """ Set annotation weights for distributed CoGAPS.

        Args:
            annotation (str list): Specify categories along the rows (cols) to use for weighted sampling.
            weight (int list): Weights associated with samplingAnnotation
        """    

        self.coparams['samplingAnnotation'] = annotation
        self.coparams['samplingWeight'] = weight

    def setFixedPatterns(self, fixedPatterns, whichMatrixFixed):
        """ Fix either 'A' or 'P' matrix to given values.

        Args:
            fixedPatterns (arr): Fix either 'A' or 'P' matrix to these values, 
            in the context of distributed CoGAPS, the first phase is skipped and 
            fixedPatterns is used for all sets - allowing manual pattern matching, 
            as well as fixed runs of standard CoGAPS.
            whichMatrixFixed (str): Either 'A' or 'P', indicating which matrix is fixed
        """ 

        self.coparams['fixedPatterns'] = fixedPatterns
        self.coparams['whichMatrixFixed'] = whichMatrixFixed
        self.gaps.useFixedPatterns = True
        self.gaps.fixedPatterns = pycogaps.Matrix(fixedPatterns)
        self.gaps.whichMatrixFixed = whichMatrixFixed

    def printParams(self):
        """ Print standard and sparsity parameters, and distributed if set.
        """        
        print('\n-- Standard Parameters --')
        print('nPatterns: ', self.gaps.nPatterns)
        print('nIterations: ', self.gaps.nIterations)
        print('seed: ', self.gaps.seed)
        print('sparseOptimization: ', self.gaps.useSparseOptimization)
        print('\n')
        print('-- Sparsity Parameters --')
        print('alpha: {:0.2f}'.format(self.gaps.alphaA))
        print('maxGibbsMass: ', self.gaps.maxGibbsMassA)
        print('\n')
        if self.gaps.runningDistributed:
            print('-- Distributed Parameters --')
            print('cut: ', self.coparams['cut'])
            print('nSets: ', self.coparams['nSets'])
            print('minNS: ', self.coparams['minNS'])
            print('maxNS: ', self.coparams['maxNS'])
            print('\n')

    def printAllParams(self):
        """ Print all GapsParameters.
        """        
        self.gaps.print()



def CoGAPS(path, params=None, nThreads=1, messages=True,
           outputFrequency=1000, uncertainty=None, checkpointOutFile="",
           checkpointInterval=0, checkpointInFile="", transposeData=False,
           BPPARAM=None, workerID=1, asynchronousUpdates=None, nSnapshots=0,
           snapshotPhase='sampling'):
            """ Python wrapper to run either standardCoGAPS or distributedCoGAPS.

            Args:
                See standardCoGAPS Args.

            Returns:
                CogapsResult: A CogapsResult object.
            """


            if params.coparams['distributed'] == 'genome-wide':
                result = distributedCoGAPS(path, params, uncertainty=None)
            else:
                result = standardCoGAPS(path, params=params, nThreads=nThreads, messages=messages,
                            outputFrequency=outputFrequency, uncertainty=uncertainty, checkpointOutFile=checkpointOutFile,
                            checkpointInterval=checkpointInterval, checkpointInFile=checkpointInFile, transposeData=transposeData,
                            BPPARAM=BPPARAM, workerID=workerID, asynchronousUpdates=asynchronousUpdates, nSnapshots=nSnapshots,
                            snapshotPhase=snapshotPhase)

            return result


def standardCoGAPS(path, params, nThreads, messages,
                    outputFrequency, uncertainty, checkpointOutFile,
                    checkpointInterval, checkpointInFile, transposeData,
                    BPPARAM, workerID, asynchronousUpdates, nSnapshots,
                    snapshotPhase):
    """ Python wrapper to run CoGAPS via bindings

    Args:
        path (str): Path to data. 
        params (CoParams, optional): CoParams object of parameters. Defaults to None.
        nThreads (int, optional): Number of threads to use. Defaults to 1.
        messages (bool, optional): Whether to print messages. Defaults to True.
        outputFrequency (int, optional): How often to output messages. Defaults to 1000.
        uncertainty (arr, optional): Optional uncertainty matrix. Defaults to None.
        checkpointOutFile (str, optional): Path to where checkpoint info should be written. Defaults to "".
        checkpointInterval (int, optional): How often to make a checkpoint. Defaults to 0.
        checkpointInFile (str, optional): Path to existing checkpoint file to run CoGAPS from. Defaults to "".
        transposeData (bool, optional): Expects genes x samples. Defaults to False.
        BPPARAM ([type], optional): BiocParallel backend . Defaults to None.
        workerID (int, optional): If calling CoGAPS in parallel the worker ID can be specified,
        only worker 1 prints output and each worker outputs when it finishes, this
        is not neccesary when using the default parallel methods (i.e. distributed
        CoGAPS) but only when the user is manually calling CoGAPS in parallel. Defaults to 1.
        asynchronousUpdates (bool, optional): Enable asynchronous updating which allows for multi-threaded runs. Defaults to None.
        nSnapshots (int, optional): How many snapshots to take in each phase, setting this to 0 disables snapshots. Defaults to 0.
        snapshotPhase (str, optional): One of "sampling", "equilibration", "all". Defaults to 'sampling'.

    Raises:
        Exception: If transposeData=True is not passed as an argument to both CoParams and CoGAPS.

    Returns:
        CogapsResult: A CogapsResult object.
    """           

    # check OpenMP support
    if isCompiledWithOpenMPSupport() is False:
        if asynchronousUpdates is not None and nThreads > 1:
            print("requesting multi-threaded version of CoGAPS but compiler did not support OpenMP")
        asynchronousUpdates = False
        nThreads = 1
    # convert sampling phase to enum
    if snapshotPhase == "sampling":
        snapshotPhase = pycogaps.GAPS_SAMPLING_PHASE
    elif snapshotPhase == "equilibration":
        snapshotPhase = pycogaps.GAPS_EQUILIBRATION_PHASE
    elif snapshotPhase == "all":
        snapshotPhase = pycogaps.GAPS_ALL_PHASES
    else:
        print("The snapshot phase you indicated is not recognized.")
        print("Please choose one of: sampling, equilibration, all")
        return

    gapsresultobj = None

    # convert data to anndata and matrix obj
    if isinstance(path, str):
        if params is not None:
            adata = toAnndata(path, params.coparams['hdfKey'], params.coparams['hdfRowKey'],
                              params.coparams['hdfColKey'], transposeData=transposeData)
        else:
            adata = toAnndata(path, transposeData=transposeData)
    else:
        adata = path

    matrix = pycogaps.Matrix(adata.X)

    if params is None:
        prm = CoParams(matrix=adata, transposeData=transposeData)
    else:
        prm = params

    opts = {
        'maxThreads': nThreads,
        'printMessages': messages,
        'outputFrequency': outputFrequency,
        'checkpointOutFile': checkpointOutFile,
        'checkpointInterval': checkpointInterval,
        'checkpointFile': checkpointInFile,
        'transposeData': transposeData,
        'workerID': workerID,
        'asynchronousUpdates': asynchronousUpdates,
        'snapshotFrequency': nSnapshots,
        'snapshotPhase': snapshotPhase,
    }
    setParams(prm, opts)

    if uncertainty is not None:
        unc = toAnndata(uncertainty)
        unc = pycogaps.Matrix(unc.X)
    else:
        unc = pycogaps.Matrix()

    if prm.coparams["subsetIndices"] is None:
        prm = getDimNames(adata, prm)

    # check data input
    checkData(adata, prm.gaps, uncertainty)
    checkInputs(uncertainty, prm)

    startupMessage(prm, path)
    gapsresultobj = pycogaps.runCogapsFromMatrix(matrix, prm.gaps, unc)
    prm.gaps.transposeData = transposeData

    if prm.gaps.transposeData != prm.coparams["transposeData"]:
        raise Exception("make sure to pass transposeData=True argument in both CoParams() and CoGAPS()")
    result = {
        "GapsResult": gapsresultobj,
        "anndata": GapsResultToAnnData(gapsresultobj, adata, prm)
    }

    show(result["anndata"])
    return result

### start distributed functions

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

    setParams(params, {'checkpointOutFile': ""})

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
    # print("FIXED PATTERNS\n", matched["consensus"])
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
        "anndata": GapsResultToAnnData(gapsresult, data, params)
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
    gapsresult = standardCoGAPS(adata, params, uncertainty, transposeData=params.coparams["transposeData"])

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
                # print(avg)
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

### end distributed functions


def GapsResultToAnnData(gapsresult, adata, prm: CoParams):
    """ Converts a CogapsResult object to anndata object.

    Args:
        gapsresult (CogapsResult): Dictionary result object.
        adata (anndata): Anndata object populated by CoGAPS.
        prm (CoParams): CoParams object.

    Returns:
        anndata: An anndata object.
    """    
    # need to subset matrices based on which dimension we're in...
    if prm.coparams['subsetDim'] == 1:
        Amean = toNumpy(gapsresult.Amean)[prm.coparams["subsetIndices"], :]
        Pmean = toNumpy(gapsresult.Pmean)
        Asd = toNumpy(gapsresult.Asd)[prm.coparams["subsetIndices"], :]
        Psd = toNumpy(gapsresult.Psd)
    else:
        Amean = toNumpy(gapsresult.Amean)
        Pmean = toNumpy(gapsresult.Pmean)[prm.coparams["subsetIndices"], :]
        Asd = toNumpy(gapsresult.Asd)
        Psd = toNumpy(gapsresult.Psd)[prm.coparams["subsetIndices"], :]
    pattern_labels = ["Pattern" + str(i) for i in range(1, prm.gaps.nPatterns + 1)]
    # load adata obs and var from Amean and Pmean results
    if len(Pmean.shape) > 2:
        Pmean = Pmean[0, :, :]
        Psd = Psd[0, :, :]
    adata.obs = pd.DataFrame(data=Amean, index=adata.obs_names, columns=pattern_labels)
    adata.var = pd.DataFrame(data=Pmean, index=adata.var_names, columns=pattern_labels)
    adata.uns["asd"] = pd.DataFrame(data=Asd, index=adata.obs_names, columns=pattern_labels)
    adata.uns["psd"] = pd.DataFrame(data=Psd, index=adata.var_names, columns=pattern_labels)
    return adata


def GapsParameters(path):
    """ Returns C++ GapsParameters object.

    Args:
        path (str): Path to data.

    Returns:
        GapsParameters: A GapsParameters object.
    """    
    return pycogaps.GapsParameters(path)


def getBuildReport():
    """ Returns information about how the package was compiled, i.e. which
    compiler/version was used, which compile time options were enabled, etc...

    Returns:
        str: String containing build report.
    """    
    return pycogaps.getBuildReport()


def isCheckpointsEnabled():
    """ Check if package was built with checkpoints enabled

    Returns:
        bool: true/false if checkpoints are enabled
    """    
    return pycogaps.isCheckpointsEnabled()


def isCompiledWithOpenMPSupport():
    """ Check if compiler supported OpenMP

    Returns:
        bool: true/false if OpenMP was supported
    """    
    return pycogaps.isCompiledWithOpenMPSupport()


def getFileInfo(path):
    """ Get info of inputted file.

    Args:
        path (str): Path to data.

    Returns:
        str: string of file info.
    """    
    return pycogaps.getFileInfo(path)


def current_milli_time():
    """ Return current time in milliseconds.

    Returns:
        int: Current time in milliseconds.
    """    
    return round(time.time() * 1000)


def setParams(paramobj: CoParams, list):
    """ Set CoParams from a list.

    Args:
        paramobj (CoParams): CoParams object.
        list (dict): Dictionary of parameter, value pairings for each parameter you wish to set.
    """    

    for (k, v) in list.items():
        setParam(paramobj, k, v)


# class CogapsParams
# constructor has default values for each parameter
def setParam(paramobj: CoParams, whichParam, value):
    """ Sets CoParams parameters.

    Args:
        paramobj (CoParams): a CoParams object
        whichParam (str): the name of the parameter you wish to change
        value ([type]): the value to set whichParam as

    Returns:
        CoParams: the modified CoParams object.
    """    

    coparam_params = ['hdfKey', 'hdfRowKey', 'hdfColKey', 'explicitSets', 'subsetDim', 'geneNames', 'sampleNames']
    if whichParam == "alpha":
        paramobj.gaps.alphaA = value
        paramobj.gaps.alphaP = value
    elif whichParam == "maxGibbsMass":
        paramobj.gaps.maxGibbsMassA = value
        paramobj.gaps.maxGibbsMassP = value
    elif whichParam in coparam_params:
        paramobj.coparams[whichParam] = value
    elif whichParam in "distributed":
        if value == "genome-wide":
            paramobj.gaps.runningDistributed = True
            paramobj.coparams['distributed'] = value
        elif (value is not None) and (value is not False):
            # print("if you wish to perform genome-wide distributed cogaps, please run setParams(params, "
            #       "\"distributed\", ""\"genome-wide\")")
            paramobj.coparams['distributed'] = 'genome-wide'
            paramobj.gaps.runningDistributed = True
        else:
            paramobj.gaps.runningDistributed = False
            paramobj.coparams['distributed'] = None
    elif whichParam in ("nSets", "cut", "minNS", "maxNS"):
        paramobj.gaps.runningDistributed = True
        print("please set \'", whichParam, "\' with setDistributedParams")
        return
    elif whichParam in ("samplingAnnotation", "samplingWeight"):
        print("please set \'", whichParam, "\' with setAnnotationWeights")
        return
    elif whichParam in ("fixedPatterns", "whichMatrixFixed"):
        print("please set \'", whichParam, "\' with setFixedPatterns")
        return
    elif whichParam in "singleCell":
        print(whichParam, " has been deprecated, this parameter will be ignored")
        return
    else:
        setattr(paramobj.gaps, whichParam, value)
    return paramobj


def getParam(paramobj, whichParam):
    """ Get parameter info and values.

    Args:
        paramobj (CoParams): a CoParams object.
        whichParam (str): which parameter to get the info of.

    Returns:
        [type]: the value of the parameter.
    """    
    return getattr(paramobj, whichParam)
