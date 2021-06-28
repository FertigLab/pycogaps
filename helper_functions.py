import pycogaps

# we can use this for testing later 
def getRetinaSubset(n=1):
    if not (1 <= n <= 4):
        raise Exception("invalide number of subsets requested")

def nrowHelper(data):
    if isinstance(data, str):
        return int(pycogaps.getFileInfo(data)["dimensions"][0])
    return data.shape[0] # assuming data is pandas dataframe

def ncolHelper(data):
    if isinstance(data, str):
        return int(pycogaps.getFileInfo(data)["dimensions"][1])
    return data.shape[1] # assuming data is pandas dataframe

def getGeneNames(data, transpose):
    if transpose:
        return getSampleNames(data, False)
    if isinstance(data, str):
        names = pycogaps.getFileInfo(data)["rowNames"]
    else:
        names = list(data.head().index) # assuming data is pandas dataframe
    if names == None or len(names) == 0:
        return ["Gene" + str(i) for i in range(1, nrowHelper(data))]
    return names

def getSampleNames(data, transpose):
    if transpose:
        return getGeneNames(data, transpose)
    if isinstance(data, str):
        names = pycogaps.getFileInfo(data)["colNames"]
    else:
        names = list(data.columns) # assuming data is pandas dataframe
    if names == None or len(names) == 0:
        return ["Sample" + str(i) for i in range(1, ncolHelper(data))]
    return names

def getDimNames(data, allParams):
    geneNames = allParams.gaps.geneNames
    sampleNames = allParams.gaps.sampleNames

    if allParams.gaps.geneNames == None:
        geneNames = getGeneNames(data, allParams.transposeData)
    if allParams.gaps.sampleNames == None:
        sampleNames = getSampleNames(data, allParams.transposeData)
    
    if allParams.transposeData:
        nGenes = ncolHelper(data)
        nSamples = nrowHelper(data)
    else:
        nGenes = nrowHelper(data)
        nSamples = ncolHelper(data)

    if allParams.gaps.subsetDim == 1:
        nGenes = len(allParams.gaps.subsetIndices)
        geneNames = geneNames[allParams.subsetIndices]
    elif allParams.gaps.subsetDim == 2:
        nSamples = len(allParams.subsetIndices)
        sampleNames = sampleNames[allParams.subsetIndices]

    if len(geneNames) != nGenes:
        raise Exception(len(geneNames), " != ", nGenes, " incorrect number of gene names given")
    if len(sampleNames) != nSamples:
        raise Exception(len(sampleNames), " != ", nSamples, " incorrect number of sample names given")

    # store processed gene/sample names directly in allParams list
    # this is an important distinction - allParams@gaps contains the
    # gene/sample names originally passed by the user, allParams contains
    # the procseed gene/sample names to be used when labeling the result
    allParams.geneNames <- geneNames
    allParams.sampleNames <- sampleNames
    return(allParams)

##### testing #####

path = "./data/GIST.csv"
prm = pycogaps.GapsParameters(path)
getDimNames(path, prm)