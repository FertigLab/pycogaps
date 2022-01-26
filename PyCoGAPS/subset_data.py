from PyCoGAPS.config import *
from PyCoGAPS.helper_functions import nrowHelper, ncolHelper, getDimNames

# explicitSets either list of indices or names
def sampleWithExplicitSets(allParams, data):
    """ Sample with user provided explicit sets

    Args:
        allParams (CoParams): a CoParams object
        data (anndata): anndata object of data

    Raises:
        Exception: If some named genes in explicitSets not found

    Returns:
        list: list of subsets
    """    
    explicit_sets = allParams.coparams['explicitSets']
    if all(isinstance(item, int) for item in explicit_sets):
        print("using provided indexed subsets")
        return explicit_sets

    if all(isinstance(item, str) for item in explicit_sets): 
        print("using provided named subsets")
        getDimNames(data, allParams)
        if allParams.coparams['distributed'] == "genome-wide":
            allNames = allParams.coparams['geneNames']
        else:
            allNames = allParams.coparams['sampleNames']
        
        for item in explicit_sets:
            if item not in allNames:
                raise Exception("some named genes in explicitSets not found")
        
        return [list(allNames).index(i) for i in explicit_sets]


def sampleWithAnnotationWeights(allParams, setSize):
    """ subset rows (cols) proportional to the user provided weights

    Args:
        allParams (CoParams): CoParams object
        setSize (int): size of each subset of the total

    Returns:
        list: list of subsets
    """    

    # samplingWeight is a dictionary (name: weight)
    weight = allParams.coparams['samplingWeight']
    sorted_weight = []
    for key in sorted(weight):
        sorted_weight.append(weight[key])
    groups = np.unique(allParams.coparams['samplingAnnotation'])
    groups = np.sort(groups)

    sets = []
    
    for i in range(allParams.coparams['nSets']):
        groupCount = np.random.choice(groups, size=setSize, replace=True, p=(sorted_weight/np.array(sorted_weight).sum()))
        subset = []
        for g in groups:
            groupNdx = np.argwhere(g == np.array(allParams.coparams['samplingAnnotation']))
            sub = np.random.choice(groupNdx.flatten(), size=sum(groupCount == g), replace=True)
            if sub.size != 0:
                subset.append(sub)
        subset = np.sort(np.concatenate(subset))
        sets.append(subset)

    return sets


def sampleUniformly(allParams, total, setSize):
    """ subset data by uniformly partitioning rows (cols)

    Args:
        allParams (CoParams): CoParams object
        total (int): total number of rows (cols) that are being paritioned
        setSize (int): size of each subset of the total

    Returns:
        list: list of subsets
    """    

    sets = [None] * (allParams.coparams['nSets'])
    remaining = np.arange(0,total)
    for n in range(allParams.coparams['nSets'] - 1):
        selected = np.random.choice(list(remaining), setSize, replace=False)
        sets[n] = np.sort(selected)
        remaining = set(remaining).difference(set(selected))
    sets[allParams.coparams['nSets'] - 1] = np.sort(list(remaining))
    return sets


def createSets(data, allParams):
    """ either genes or samples or partitioned depending on the type
    of distributed CoGAPS (i.e. genome-wide or single-cell)

    Args:
        data (anndata): anndata object of data
        allParams (CoParams): CoParams object

    Raises:
        Exception: If nSets does not match number of explicit sets given

    Returns:
        list: list of sets
    """    
    subsetRows = allParams.gaps.transposeData != allParams.coparams['distributed'] == "genome-wide"
    if subsetRows:
        total = nrowHelper(data)
    else:
        total = ncolHelper(data)
    setSize = math.floor(total / allParams.coparams['nSets'])

    print('Creating subsets...')

    if allParams.coparams['explicitSets'] is not None:
        if len(allParams.coparams['explicitSets']) != allParams.coparams['nSets']:
            raise Exception('nSets does not match number of explicit sets given')
        sets = sampleWithExplicitSets(allParams, data)

    elif allParams.coparams['samplingAnnotation'] is not None:
        print('sampling with annotation weights')
        sets = sampleWithAnnotationWeights(allParams, setSize)

    else:
        sets = sampleUniformly(allParams, total, setSize)

    # print('shape: ', sets.shape)
    # print('set sizes (min, mean, max): (',
    #     min(len(sets)), ', ',
    #     np.mean(len(sets)), ', ',
    #     max(len(sets)), ')\n')

    return sets