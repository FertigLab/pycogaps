''' 
this script reads parameters from the command line to run CoGAPS
supports integration with genepattern notebook
'''

if __name__ == '__main__':
    from parameters import *
    from pycogaps_main import CoGAPS
    import pickle
    import argparse
    
    print("This vignette was built using pycogaps version", getVersion())

    ''' 
    command line args which are all parameters to CoGAPS
        - only --path arg is required
        - all other args are optional, have default values
    '''
    parser = argparse.ArgumentParser()

    ## initial params ##
    # path to data
    parser.add_argument('--path', type=str, required=True)
    # result output file name (output saved as a .pkl file)
    parser.add_argument('--resultFile', type=str, default='result.pkl')
    
    ## standard params ##
    # number of patterns CoGAPS will learn
    parser.add_argument('--nPatterns', type=int, default=3)
    # number of iterations for each phase of the algorithm
    parser.add_argument('--nIterations', type=int, default=1000)
    # random number generator seed
    parser.add_argument('--seed', type=int, default=0)
    # speeds up performance with sparse data (roughly >80% of data is zero), note this can only be used with the default uncertainty
    parser.add_argument('--useSparseOptimization', type=bool, default=False)
    
    ## run params ##
    # maximum number of threads to run on
    parser.add_argument('--nThreads', type=bool, default=1)
    # T/F for displaying output
    parser.add_argument('--messages', type=bool, default=True)
    # number of iterations between each output (set to 0 to disable status updates)
    parser.add_argument('--outputFrequency', type=int, default=500)
    # uncertainty matrix - either a matrix or a supported file type
    parser.add_argument('--uncertainty', type=str, default=None)
    # name of the checkpoint file to create
    parser.add_argument('--checkpointOutFile', type=str, default='gaps_checkpoint.out')
    # if this is provided, CoGAPS runs from the checkpoint contained in this file
    parser.add_argument('--checkpointInFile', type=str, default="")
    # T/F for transposing data while reading it in - useful for data that is stored as samples x genes since CoGAPS requires data to be genes x samples    
    parser.add_argument('--transposeData', type=bool, default=False)
    # if calling CoGAPS in parallel the worker ID can be specified
    parser.add_argument('--workerID', type=int, default=1)
    # enable asynchronous updating which allows for multi-threaded runs
    parser.add_argument('--asynchronousUpdates', type=bool, default=False)
    # how many snapshots to take in each phase, setting this to 0 disables snapshots
    parser.add_argument('--nSnapshots', type=int, default=0)
    # which phase to take snapsjots in e.g. "equilibration", "sampling", "all"
    parser.add_argument('--snapshotPhase', type=str, default='sampling', choices=['sampling', 'equilibration', 'all'])
    
    ## sparsity params ##
    # sparsity parameter for feature matrix
    parser.add_argument('--alphaA', type=float, default=0.01)
    # sparsity parameter for sample matrix
    parser.add_argument('--alphaP', type=float, default=0.01)
    # atomic mass restriction for feature matrix
    parser.add_argument('--maxGibbsMassA', type=float, default=100)
    # atomic mass restriction for sample matrix
    parser.add_argument('--maxGibbsMassP', type=float, default=100)
    
    ## distributed params ##
    #  either null or genome-wide
    parser.add_argument('--distributed', type=str, default=None)
    # number of sets to break data into
    parser.add_argument('--nSets', type=int, default=4)
    # number of branches at which to cut dendrogram used in pattern matching
    # default: nPatterns
    parser.add_argument('--cut', type=int, default=None)
    # minimum of individual set contributions a cluster must contain
    # default: math.ceil(cut / 2)
    parser.add_argument('--minNS', type=int, default=None)
    # maximum of individual set contributions a cluster can contain
    # default: minNS + nSets
    parser.add_argument('--maxNS', type=int, default=None)
    # specify subsets by index or name
    parser.add_argument('--explicitSets', type=list, default=None)
    # specify categories along the rows (cols) to use for weighted sampling
    parser.add_argument('--samplingAnnotation', type=list, default=None)
    # weights associated with  samplingAnnotation
    parser.add_argument('--samplingWeight', type=list, default=None)
    
    ## additional params ##
    # set of indices to use from the data
    parser.add_argument('--subsetIndices', type=set, default=None)
    # which dimension (0=rows, 1=cols) to subset
    parser.add_argument('--subsetDim', type=int, default=0, choices=[0,1])
    # vector of names of genes in data
    parser.add_argument('--geneNames', type=list, default=None)
    # vector of names of samples in data
    parser.add_argument('--sampleNames', type=list, default=None)
    # fix either 'A' or 'P' matrix to these values, in the context of distributed CoGAPS, the first phase is skipped and `fixedPatterns: 
    # is used for all sets allowing manual pattern matching, as well as fixed runs of standard CoGAPS
    parser.add_argument('--fixedPatterns', default=None)
    # either 'A' or 'P', indicating which matrix is fixed
    parser.add_argument('--whichMatrixFixed', type=str, default=None, choices=['A', 'P'])
    # whether or not to take PUMP samples
    parser.add_argument('--takePumpSamples', type=bool, default=False)
    # for reading .h5 files
    parser.add_argument('--hdfKey', type=str, default=None)
    # for reading .h5 files
    parser.add_argument('--hdfRowKey', type=str, default=None)
    # for reading .h5 files
    parser.add_argument('--hdfColKey', type=str, default=None)
    
    initial_params = ["path", "resultFile"]
    
    standard_params = ["nPatterns", "nIterations", "seed", "useSparseOptimization"]

    run_params = ["nThreads", "messages", "outputFrequency", "uncertainty", "checkpointOutFile", "checkpointInterval", 
                "checkpointInFile", "transposeData", "workerID", "asynchronousUpdates", 
                "nSnapshots", "snapshotPhase"]

    sparsity_params = ["alphaA", "alphaP", "maxGibbsMassA", "maxGibbsMassP"]
    
    distributed_params = ["distributed", "nSets", "cut", "minNS", "maxNS", 
                "explicitSets", "samplingAnnotation", "samplingWeight"]
   
    additional_params = ["subsetIndices", "subsetDim", "geneNames", "sampleNames",   
                "fixedPatterns", "whichMatrixFixed", "takePumpSamples", 
                "hdfKey", "hdfRowKey", "hdfColKey"]
    
    '''
    parse all args and set as parameters for CoGAPS
    '''
    args = parser.parse_args()

    data_path = args.path
    
    params = CoParams(path=data_path, transposeData=args.transposeData, 
                      hdfKey=args.hdfKey, hdfRowKey=args.hdfRowKey,
                      hdfColKey=args.hdfColKey)

    prm_dict = vars(args)

    for k,v in prm_dict.items():
        if ((k not in initial_params) and (k not in distributed_params) and (k not in ("fixedPatterns", "uncertainty"))):
            setParam(params, k, v)
    
     # set fixed patterns from additional params
    if args.fixedPatterns is not None:
        params.setFixedPatterns(fixedPatterns=args.fixedPatterns, whichMatrixFixed=args.whichMatrixFixed)

    # set distributed parameters
    setParam(params, 'distributed', args.distributed)
    if args.distributed is not None:
        params.setAnnotationWeights(annotation=args.samplingAnnotation, weight=args.samplingWeight)
        params.setDistributedParams(nSets=args.nSets, cut=args.cut, minNS=args.minNS, maxNS=args.maxNS)

    '''
    run CoGAPS, save result
    '''
    result = CoGAPS(data_path, params, uncertainty=args.uncertainty)

    # save CoGAPS result
    print("Pickling...", end='\r')
    pickle.dump(result, open(args.resultFile, "wb"))
    print("Pickling complete!")
