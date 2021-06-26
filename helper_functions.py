import pycogaps

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
