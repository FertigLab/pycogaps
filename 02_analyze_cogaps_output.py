import pickle

with open('./data/PDACCoGAPSResult031422.pkl', 'rb') as fp:
    cogaps_result = pickle.load(fp)

# cogaps_result
#
# AnnData object with n_obs × n_vars = 25442 × 15219 obs: 0, 1, 2, 3, 4, 5, 6, 7 var: 0, 1, 2, 3, 4, 5,
# 6, 7 uns: 'asd', 'psd', 'atomhistoryA', 'atomhistoryP', 'averageQueueLengthA', 'averageQueueLengthP',
# 'chisqHistory', 'equilibrationSnapshotsA', 'equilibrationSnapshotsP', 'meanChiSq', 'meanPatternAssignment',
# 'pumpMatrix', 'samplingSnapshotsA', 'samplingSnapshotsP', 'seed', 'totalRunningTime', 'totalUpdates'

import scanpy as sc
adata = sc.read_h5ad("/Users/jeanette/fertiglab/PDAC_Atlas_Pipeline/PDAC.h5ad")
# copy over gene names (should probably do this earlier.. oh well)
cogaps_result.var_names = cogaps_result.var_names = adata.var["gene_short_name"]
from PyCoGAPS import *
from PyCoGAPS.analysis_functions import *
from PyCoGAPS.helper_functions import *

adata = cogaps_result
sc.pp.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable]

# scale data and compute PCA
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

# find neighbor embeddings and run UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# add categorical annotations to observation matrix for umap plot
# adata.obs['majorCluster']=list(majorCluster)

# plot pattern amplitude on UMAP
sc.pl.umap(adata, color="")
