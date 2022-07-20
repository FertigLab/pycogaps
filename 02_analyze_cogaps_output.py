import pickle

with open('./data/testing032522.pkl', 'rb') as fp:
    cogaps_result = pickle.load(fp)

# cogaps_result
#
# AnnData object with n_obs × n_vars = 25442 × 15219 obs: 0, 1, 2, 3, 4, 5, 6, 7 var: 0, 1, 2, 3, 4, 5,
# 6, 7 uns: 'asd', 'psd', 'atomhistoryA', 'atomhistoryP', 'averageQueueLengthA', 'averageQueueLengthP',
# 'chisqHistory', 'equilibrationSnapshotsA', 'equilibrationSnapshotsP', 'meanChiSq', 'meanPatternAssignment',
# 'pumpMatrix', 'samplingSnapshotsA', 'samplingSnapshotsP', 'seed', 'totalRunningTime', 'totalUpdates'

import scanpy as sc
import scipy
# adata = sc.read_h5ad("/Users/jeanette/fertiglab/PDAC_Atlas_Pipeline/PDAC.h5ad")

# copy over gene names (should probably do this earlier.. oh well)
# cogaps_result.var_names = cogaps_result.var_names = adata.var["gene_short_name"]
from PyCoGAPS import *
from PyCoGAPS.analysis_functions import *
from PyCoGAPS.helper_functions import *

adata = cogaps_result
adata.X = scipy.sparse.csr_matrix(adata.X)
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
patterns = list(adata.obs.columns)[0:7]
sc.pl.umap(adata, color=patterns)



# PATTERNS IMPORTED FROM R
samplefactors = pd.read_csv("/Users/jeanette/fertiglab/PDAC_Atlas_Pipeline/PDACsampleFactors.csv", index_col=0)
featureloadings = pd.read_csv("/Users/jeanette/fertiglab/PDAC_Atlas_Pipeline/PDACfeatureLoadings.csv", index_col=0)
epiMtx = sc.read_mtx("/Users/jeanette/Downloads/epiMat.mtx")
epiMtx.X = scipy.sparse.csr_matrix(epiMtx.X)
epiMtx.obs = featureloadings
epiMtx.var = samplefactors
plaindata = epiMtx
epiMtx = epiMtx.T

coldata = pd.read_csv("/Users/jeanette/fertiglab/PDAC_Atlas_Pipeline/PDACcoldata.csv")
epiMtx.obs["cell type"] = list(coldata["TN_assigned_cell_type_immune_specific"])

sc.pp.normalize_total(epiMtx)
# sc.pp.log1p(epiMtx)

sc.pp.scale(epiMtx)
sc.tl.pca(epiMtx, svd_solver='arpack')
sc.pp.neighbors(epiMtx, n_pcs=25, n_neighbors=50)
sc.tl.leiden(epiMtx, resolution=0.1)
sc.tl.umap(epiMtx, min_dist=.1)
sc.pl.umap(epiMtx, color=list(epiMtx.obs.columns), legend_loc="on data")