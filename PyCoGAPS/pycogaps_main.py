import anndata

from PyCoGAPS.config import *
from PyCoGAPS.helper_functions import *
from PyCoGAPS.subset_data import *
from PyCoGAPS.parameters import *
from PyCoGAPS.distributed_functions import *

import multiprocessing


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


def standardCoGAPS(path, params=None, nThreads=1, messages=True,
           outputFrequency=1000, uncertainty=None, checkpointOutFile="",
           checkpointInterval=0, checkpointInFile="", transposeData=False,
           BPPARAM=None, workerID=1, asynchronousUpdates=None, nSnapshots=0,
           snapshotPhase='sampling'):
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
    # # convert sampling phase to enum
    # if snapshotPhase == "sampling":
    #     snapshotPhase = pycogaps.GAPS_SAMPLING_PHASE
    # elif snapshotPhase == "equilibration":
    #     snapshotPhase = pycogaps.GAPS_EQUILIBRATION_PHASE
    # elif snapshotPhase == "all":
    #     snapshotPhase = pycogaps.GAPS_ALL_PHASES
    # else:
    #     print("The snapshot phase you indicated is not recognized.")
    #     print("Please choose one of: sampling, equilibration, all")
    #     return

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

    '''
    make sure uncertainty matrix processed the same way as adata input
    '''
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
    # no longer returning the legacy formatted object
    result = GapsResultToAnnData(gapsresultobj, adata, prm)
    show(result)
    return result


def distributedCoGAPS(path, params, uncertainty=None):
    if isinstance(path, str):
        data = toAnndata(path, hdf_counts_key=params.coparams["hdfKey"],
                                          hdf_dim1_key=params.coparams["hdfRowKey"],
                                          hdf_dim2_key=params.coparams["hdfColKey"],
                                          transposeData=params.coparams["transposeData"])
    else:
        data = path
    sets = createSets(data, params)
    if min(map(len, sets)) < params.gaps.nPatterns:
        warnings.warn("Data subset dimension less than nPatterns. Aborting.")
        return 1

    # setParams(params, {'checkpointOutFile': ""})

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
            # print("POOL IS NOW CLOSED")
            if params.coparams['distributed'] == "genome-wide":
                unmatched = np.array(result[0].var)
                for i in range(1, len(result)):
                    arr = np.array(result[i].var)
                    unmatched = np.append(unmatched, arr, axis=1)
            else:
                unmatched = np.array(result[0].obs)
                for i in range(1, len(result)):
                    arr = np.array(result[i].obs)
                    unmatched = np.append(unmatched, arr, axis=1)
            print("Matching patterns across subsets...\n")
            matched = findConsensusMatrix(unmatched, params)
    else:
        matched = params.gaps.fixedPatterns

    params.gaps.nPatterns = matched["consensus"].shape[1]
    params.gaps.useFixedPatterns = True
    params.gaps.fixedPatterns = pycogaps.Matrix(matched["consensus"])

    # print('=== DEBUG MATRIX FP ===')
    # print('np FP: ', matched["consensus"])
    # print('Matrix FP: ', toNumpy(params.gaps.fixedPatterns))
    # print('=== END DEBUG ===')

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

    stitched = stitchTogether(finalresult, result, params, sets, data)


    adata = data
    adata.obs = stitched["Amean"]
    adata.var = stitched["Pmean"]
    adata.uns["asd"] = stitched["Asd"]
    adata.uns["psd"] = stitched["Psd"]
    adata.uns["atomhistoryA"] = finalresult[0].uns["atomhistoryA"]
    adata.uns["atomhistoryP"] = finalresult[0].uns["atomhistoryP"]
    adata.uns["averageQueueLengthA"] = finalresult[0].uns["averageQueueLengthA"]
    adata.uns["averageQueueLengthP"] = finalresult[0].uns["averageQueueLengthP"]
    adata.uns["chisqHistory"] = finalresult[0].uns["chisqHistory"]
    adata.uns["equilibrationSnapshotsA"] = finalresult[0].uns["equilibrationSnapshotsA"]
    adata.uns["equilibrationSnapshotsP"] = finalresult[0].uns["equilibrationSnapshotsP"]
    adata.uns["meanChiSq"] = finalresult[0].uns["meanChiSq"]
    adata.uns["meanPatternAssignment"] = finalresult[0].uns["meanPatternAssignment"]
    adata.uns["pumpMatrix"] = finalresult[0].uns["pumpMatrix"]
    adata.uns["samplingSnapshotsA"] = finalresult[0].uns["samplingSnapshotsA"]
    adata.uns["samplingSnapshotsP"] = finalresult[0].uns["samplingSnapshotsP"]
    adata.uns["seed"] = finalresult[0].uns["seed"]
    adata.uns["totalRunningTime"] = finalresult[0].uns["totalRunningTime"]
    adata.uns["totalUpdates"] = finalresult[0].uns["totalUpdates"]

    # gapsresult = pycogaps.GapsResult
    # if params.coparams["distributed"] == "genome-wide":
    #     gapsresult.Amean = pycogaps.Matrix((stitched["Amean"]))
    #     gapsresult.Asd = pycogaps.Matrix((stitched["Asd"]))
    #     gapsresult.Pmean = pycogaps.Matrix((stitched["Pmean"]))
    #     gapsresult.Psd = pycogaps.Matrix((stitched["Psd"]))
    #
    # else:
    #     gapsresult.Amean = pycogaps.Matrix((stitched["Amean"]))
    #     gapsresult.Asd = pycogaps.Matrix((stitched["Asd"]))
    #     gapsresult.Pmean = pycogaps.Matrix((stitched["Pmean"]))
    #     gapsresult.Psd = pycogaps.Matrix((stitched["Psd"]))

    return adata



def callInternalCoGAPS(paramlst):
    # take out parameters passed as a list to the worker process
    # print("IN CALL INTERNAL COGAPS")
    path = paramlst[0]
    params = paramlst[1]
    workerID = paramlst[2]
    subsetIndices = paramlst[3]
    uncertainty = paramlst[4]
    if isinstance(path, str):
        adata = toAnndata(path)
    else:
        adata = path
    if subsetIndices is None:
        print("No subset indices provided; generating random sets...")
        subsetIndices = createSets(adata, params)
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
        adata = adata[subsetIndices, :]
        params.coparams['subsetDim'] = 2
       
    params.coparams['subsetIndices'] = subsetIndices
    params.gaps.workerID = workerID
    params.gaps.asynchronousUpdates = False
#     params.gaps.maxThreads = 1
    print("Calling internal CoGAPS...\n")
    gapsresult = standardCoGAPS(adata, params, uncertainty, transposeData=params.coparams["transposeData"])

    return gapsresult

