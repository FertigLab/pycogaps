import time
import math


def current_milli_time():
    return round(time.time() * 1000)


# class CogapsParams
# constructor has default values for each parameter
class CogapsParams(object):
    def __init__(self) -> object:
        self.nPatterns = 7
        self.nIterations = 50000
        self.alphaA = 0.01
        self.alphaP = 0.01
        self.maxGibbsMassA = 100
        self.massGibbMassP = 100
        self.seed = current_milli_time()
        self.sparseOptimization = False
        self.cut = self.nPatterns
        self.nSets = 4
        self.minNS = math.ceil(self.nSets / 2)
        self.maxNS = self.minNS + self.nSets
        self.explicitSets = None
        self.samplingAnnotation = None
        self.samplingWeight = None
        self.subsetIndcies = None
        self.subsetDim = 0
        self.geneNames = None
        self.sampleNames = None
        self.fixedPatterns = None
        self.whichMatrixFixed = 'N'
        self.takePumpSamples = False

    def __setParams__(self, whichParam, value):
        if(whichParam == "nPatterns"):
            self.nPatterns=value
        elif(whichParam == "nIterations"):
            self.nIterations = value
        return self







