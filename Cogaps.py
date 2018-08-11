import CogapsPy
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
import os
from zipfile import ZipFile

def geneThreshold(rankMat, patt):
    cutRank = rankMat.shape[0]
    for i in range(rankMat.shape[0]):
        for j in range(rankMat.shape[1]):
            if (j != patt) and (rankMat[i,j] <= rankMat[i,patt]):
                cutRank = min(cutRank, max(0.0, rankMat[i,patt]-1))
    return cutRank

def patternMarkers(Amat, Pmat, pumpThreshold):
    numGenes, numPatterns = Amat.shape
    statMat = np.zeros(shape=(numGenes, numPatterns))

    #normalize matrices
    factors = np.sum(Pmat, axis=1)
    factors[factors == 0] = 1.0
    Pmat = Pmat / factors[:,None]
    scales = np.amax(Pmat, axis=1)
    Amat = Amat * factors * scales

    #compute sstat (can be further numpy optimized)
    sstat = np.zeros(shape=(numGenes, numPatterns))
    lp = np.zeros(shape=(numPatterns))
    diff = np.zeros(shape=(numPatterns))
    for j in range(numPatterns):
        lp[j] = 1.0
        for i in range(numGenes):
            geneMax = np.amax(Amat[i])
            if geneMax > 0.0:
                diff = Amat[i]/geneMax - lp
            else:
                diff = lp * -1.0
                sstat[i,j] = np.asscalar(np.sqrt(np.dot(diff,diff)))
                lp[j] = 0.0

    #update PUMP matrix
    if pumpThreshold == 0: # == PUMP_UNIQUE
        for i in range(numGenes):
            minIndex = np.asscalar(np.argmin(sstat[i]))
            statMat[i,minIndex] += 1
    elif pumpThreshold == 1: # == PUMP_CUT
        rankMat = np.zeros(shape=(numGenes, numPatterns))
        for j in range(numPatterns):
            rankMat[:,j] = np.argsort(sstat[:,j])
        for j in range(numPatterns):
            cutRank = geneThreshold(rankMat, j)
            for i in range(numGenes):
                if (rankMat[i,j] <= cutRank):
                    statMat[i,j] += 1
    return (statMat, sstat)

class Cogaps:
    def __init__(self, load=False, dataPath='', numPatterns=0, maxIterations=0, runName=None):
        '''
        if isinstance(data, str):
             data needs to be parsed from file
        elif isinstance(data, np.ndarray):
             data is already in numpy array
        else:
             throw input error
        '''
        if load is False:
            if runName is None:
                self.runName = os.path.splitext(os.path.basename(dataPath))[0] + ".cogaps_" + str(numPatterns)
            else:
                self.runName = runName
            self.Amean, self.Asd, self.Pmean, self.Psd = CogapsPy.CoGAPS(dataPath, numPatterns, maxIterations)
        else:
            self.Amean = None
            self.Asd = None
            self.Pmean = None
            self.Psd = None
            self.runName = os.path.splitext(os.path.basename(dataPath))[0]
            zip = ZipFile(dataPath, 'r')
            zip.extractall('.')
            zip.close()
            self.Amean = np.loadtxt(self.runName + "_Amean.tsv", delimiter='\t')
            self.Asd = np.loadtxt(self.runName + "_Asd.tsv", delimiter='\t')
            self.Pmean = np.loadtxt(self.runName + "_Pmean.tsv", delimiter='\t')
            self.Psd = np.loadtxt(self.runName + "_Psd.tsv", delimiter='\t')
            os.remove(self.runName + "_Amean.tsv")
            os.remove(self.runName + "_Asd.tsv")
            os.remove(self.runName + "_Pmean.tsv")
            os.remove(self.runName + "_Psd.tsv")
            print('Cogaps object loaded with results from ' + dataPath + '.')
            

    def saveResults(self):
        files = ['./'+ self.runName + "_Amean.tsv", './'+ self.runName + "_Asd.tsv", './'+ self.runName + "_Pmean.tsv", './'+ self.runName + "_Psd.tsv"]
        np.savetxt(self.runName + "_Amean.tsv", self.Amean, delimiter='\t')
        np.savetxt(self.runName + "_Asd.tsv", self.Asd, delimiter='\t')
        np.savetxt(self.runName + "_Pmean.tsv", self.Pmean, delimiter='\t')
        np.savetxt(self.runName + "_Psd.tsv", self.Psd, delimiter='\t')

        with ZipFile(self.runName + '.zip', 'w') as zip:
            for file in files:
                zip.write(file)
                os.remove(file)

        print('Cogaps results saved and zipped in ' + self.runName + '.zip.')
        
        
        
    def graphPatternMarkerStats(self, pumpThreshold):
        statMat, sstat = patternMarkers(self.Amean, self.Pmean, pumpThreshold)

        origMat = np.matmul(self.Amean, self.Pmean)
        numGenes, numSamples = origMat.shape
        numPatterns = self.Amean.shape[1]
        #geneNames = np.array(AmeanRowNames, dtype=object)

        rearrOrigMat = None
        #allRearrRowNames = None
        for j in range(numPatterns):
            whichGenes = statMat[:,j]
            selectedGenes = origMat[whichGenes==1,:]
            selectedSstat = sstat[whichGenes==1,j]
            #selectedRowNames = geneNames[whichGenes==1,0]

            rearrIndices = np.argsort(selectedSstat)
            rearrGenes = np.empty(shape=(selectedGenes.shape[0], numSamples))
            #rearrRowNames = np.empty(shape=(selectedRowNames.shape[0],1), dtype=object)
            for i in range(selectedGenes.shape[0]):
                rearrGenes[i] = selectedGenes[rearrIndices[i]]
                #rearrRowNames[i] = selectedRowNames[rearrIndices[i]]
            if rearrOrigMat is None:
                rearrOrigMat = rearrGenes
                #allRearrRowNames = rearrRowNames
            else:
                rearrOrigMat = np.concatenate((rearrOrigMat, rearrGenes), axis=0)
                #allRearrRowNames = np.concatenate((allRearrRowNames, rearrRowNames), axis=0)

        normRearrOrigMat = rearrOrigMat / (np.abs(rearrOrigMat).sum(axis=1)[:,None])

        fig, ax = plt.subplots(figsize=(100,3))
        im = ax.imshow(normRearrOrigMat.T, cmap='afmhot', interpolation=None, aspect='auto')
        ax.set_yticks(np.arange(self.Pmean.shape[1]))
        ax.set_xticks(np.arange(self.Amean.shape[0])[0::10])
        #ax.set_yticklabels(PmeanColNames)
        #ax.set_xticklabels(allRearrRowNames[0::10])
        ax.tick_params(top=True, bottom=False,
                       labeltop=True, labelbottom=False)
        plt.setp(ax.get_xticklabels(), rotation=-90, ha="right",
             rotation_mode="anchor")
        plt.show()
        
    def graphPmean(self):
        plt.figure(figsize=(10,10))
        lines = []
        for i in range(self.Pmean.shape[0]):
            line, = plt.plot(self.Pmean[i])
            lines.append(line)
        names = []
        for i in range(len(lines)):
            names.append("Pattern " + str(i+1))
        plt.legend(lines, names, loc='best', shadow=True, fancybox=True)
        plt.xlabel('Index')
        plt.ylabel('Relative Amplitude')
        plt.show()
