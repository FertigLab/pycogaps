import CogapsPy
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter

def gene_threshold(rank_mat, pat):
    cut_rank = rank_mat.shape[0]
    for i in range(rank_mat.shape[0]):
        for j in range(rank_mat.shape[1]):
            if (j != pat) and (rank_mat[i,j] <= rank_mat[i,pat]):
                cut_rank = min(cut_rank, max(0.0, rank_mat[i,pat]-1))
    return cut_rank

def pattern_markers(P_mat, A_mat, pump_threshold):
    num_genes, num_patterns = A_mat.shape
    stat_mat = np.zeros(shape=(num_genes, num_patterns))

    #normalize matrices
    factors = np.sum(P_mat, axis=1)
    factors[factors == 0] = 1.0
    P_mat = P_mat / factors[:,None]
    scales = np.amax(P_mat, axis=1)
    A_mat = A_mat * factors * scales

    #compute sstat (can be further numpy optimized)
    sstat = np.zeros(shape=(num_genes, num_patterns))
    lp = np.zeros(shape=(num_patterns))
    diff = np.zeros(shape=(num_patterns))
    for j in range(num_patterns):
        lp[j] = 1.0
        for i in range(num_genes):
            gene_max = np.amax(A_mat[i])
            if gene_max > 0.0:
                diff = A_mat[i]/gene_max - lp
            else:
                diff = lp * -1.0
                sstat[i,j] = np.asscalar(np.sqrt(np.dot(diff,diff)))
                lp[j] = 0.0

    #update PUMP matrix
    if pump_threshold == 0: # == PUMP_UNIQUE
        for i in range(num_genes):
            min_index = np.asscalar(np.argmin(sstat[i]))
            stat_mat[i,min_index] += 1
    elif pump_threshold == 1: # == PUMP_CUT
        rank_mat = np.zeros(shape=(num_genes, num_patterns))
        for j in range(num_patterns):
            rank_mat[:,j] = np.argsort(sstat[:,j])
        for j in range(num_patterns):
            cut_rank = gene_threshold(rank_mat, j)
            for i in range(num_genes):
                if (rank_mat[i,j] <= cut_rank):
                    stat_mat[i,j] += 1
    return (stat_mat, sstat)



class Cogaps:
    def __init__(self, data):
        '''
        if isinstance(data, str):
             data needs to be parsed from file
        elif isinstance(data, np.ndarray):
             data is already in numpy array
        else:
             throw input error
        '''
        self.Amean, self.Asd, self.Pmean, self.Psd = CogapsPy.CoGAPS(data)
        
        
    def graphPatternMarkerStats(self, pumpThreshold):
        stat_mat, sstat = pattern_markers(self.Pmean, self.Amean, pumpThreshold)

        orig_mat = np.matmul(self.Amean, self.Pmean)
        num_genes, num_samples = orig_mat.shape
        num_patterns = self.Amean.shape[1]
        #gene_names = np.array(Amean_rowNames, dtype=object)

        rearr_orig_mat = None
        #all_rearr_rowNames = None
        for j in range(num_patterns):
            which_genes = stat_mat[:,j]
            selected_genes = orig_mat[which_genes==1,:]
            selected_sstat = sstat[which_genes==1,j]
            #selected_rowNames = gene_names[which_genes==1,0]

            rearr_indices = np.argsort(selected_sstat)
            rearr_genes = np.empty(shape=(selected_genes.shape[0], num_samples))
            #rearr_rowNames = np.empty(shape=(selected_rowNames.shape[0],1), dtype=object)
            for i in range(selected_genes.shape[0]):
                rearr_genes[i] = selected_genes[rearr_indices[i]]
                #rearr_rowNames[i] = selected_rowNames[rearr_indices[i]]
            if rearr_orig_mat is None:
                rearr_orig_mat = rearr_genes
                #all_rearr_rowNames = rearr_rowNames
            else:
                rearr_orig_mat = np.concatenate((rearr_orig_mat, rearr_genes), axis=0)
                #all_rearr_rowNames = np.concatenate((all_rearr_rowNames, rearr_rowNames), axis=0)

        norm_rearr_orig_mat = rearr_orig_mat / (np.abs(rearr_orig_mat).sum(axis=1)[:,None])

        fig, ax = plt.subplots(figsize=(100,3))
        im = ax.imshow(norm_rearr_orig_mat.T, cmap='afmhot', interpolation=None, aspect='auto')
        ax.set_yticks(np.arange(self.Pmean.shape[1]))
        ax.set_xticks(np.arange(self.Amean.shape[0])[0::10])
        #ax.set_yticklabels(Pmean_colNames)
        #ax.set_xticklabels(all_rearr_rowNames[0::10])
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
