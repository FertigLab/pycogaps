import time
import math


def current_milli_time():
    return round(time.time() * 1000)


# class CogapsParams
# constructor has default values for each parameter
class CogapsParams(object):
    def __init__(self) -> object:
        self.distributed = None
        self.nPatterns = 7
        self.nIterations = 50000
        self.alphaA = 0.01
        self.alphaP = 0.01
        self.maxGibbsMassA = 100
        self.maxGibbMassP = 100
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

    def setParam(self, whichParam, value):
        if whichParam == "nPatterns":
            self.nPatterns = value
        elif whichParam == "nIterations":
            self.nIterations = value
        elif whichParam == "alpha":
            self.alphaA = value
            self.alphaP = value
        elif whichParam == "maxGibbsMass":
            self.maxGibbsMassA = value
            self.maxGibbsMassP = value
        elif whichParam in ("nSets", "cut", "minNS", "maxNS"):
            print("please set \'", whichParam, "\' with setDistributedParams")
            return
        elif whichParam in ("samplingAnnotation", "samplingWeight"):
            print("please set \'", whichParam, "\' with setAnnotationWeights")
            return
        elif whichParam in ("fixedPatterns", "whichMatrixFixed"):
            print("please set \'", whichParam, "\' with setFixedPatterns")
            return
        elif whichParam == "nPatterns":
            self.nPatterns = value
            self.cut = min(self.cut, self.nPatterns)
        elif whichParam == "distributed":
            if value == "none":
                self.distributed = None
            else:
                self.distributed = value
        elif whichParam in "singleCell":
            print(whichParam, " has been deprecated, this parameter will be ignored")
            return
        else:
            setattr(self, whichParam, value)
        return self

    def show(self):
        print("-- Standard Parameters --\n")
        print("nPatterns           ", self.nPatterns, "\n")
        print("nIterations         ", self.nIterations, "\n")
        print("seed                ", self.seed, "\n")
        print("sparseOptimization  ", self.sparseOptimization, "\n")
        if self.distributed is not None:
            print("distributed         ", self.distributed, "\n")
        print("\n")
        print("-- Sparsity Parameters --\n")
        if self.alphaA==self.alphaP:
            print("alpha         ", self.alphaA, "\n")
        else:
            print("alphaA        ", self.alphaA, "\n")
            print("alphaP        ", self.alphaP, "\n")
        if self.maxGibbsMassA==self.maxGibbsMassP:
            print("maxGibbsMass  ", self.maxGibbsMassA, "\n")
        else:
            print("maxGibbsMassA ", self.maxGibbsMassA, "\n")
            print("maxGibbsMassP ", self.maxGibbsMassP, "\n")
        if self.distributed is not None:
            print("\n")
            print("-- Distributed CoGAPS Parameters --", "\n")
            print("nSets         ", self.nSets, "\n")
            print("cut           ", self.cut, "\n")
            print("minNS         ", self.minNS, "\n")
            print("maxNS         ", self.maxNS, "\n")
        if self.geneNames is not None:
            print("\n")
            print(len(self.geneNames), "gene names provided\n")
            print("first gene name:", self.geneNames[1], "\n")
        if self.sampleNames is not None:
            print("\n")
            print(len(self.sampleNames), "sample names provided\n")
            print("first sample name:", self.sampleNames[1], "\n")