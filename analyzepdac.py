import anndata
import pandas as pd
import scanpy as sc
# import pandas as pd

cogapsresult = anndata.read_h5ad("data/cogapsresult.h5ad")
# pdac = anndata.read_h5ad("data/PDAC.h5ad")
# pdac_peng_epi = anndata.read_h5ad("data/PDAC_Peng_Epi.h5ad")

# coldata = pd.read_csv("data/PDACcoldata.csv")
# Index(['Unnamed: 0', 'barcode_raw', 'celltype', 'sample_ID',
#        'sample_ID_celltype', 'TN', 'TN_manuscript', 'manuscript', 'nCount_RNA',
#        'nFeature_RNA', 'percent.mt', 'Size_Factor', 'TN_cluster_resolution_5',
#        'TN_assigned_cell_type', 'TN_assigned_cell_type_immune',
#        'TN_assigned_cell_type_immune_specific',
#        'TN_assigned_cell_type_immune_broad', 'cc', 'ccstage',
#        'Classifier_T_duct', 'Classifier_T_Fibroblast_only',
#        'Classifier_T_Fibroblast_Stellate'],
#       dtype='object')
#
# adata = pycogapsresult
# # get readable gene names from original object
# adata_original = sc.read_h5ad("data/PDAC_Peng_Epi.h5ad").T

# adata = pycogapsresult.T

# adata.obs["cell type"] = list(coldata["TN_assigned_cell_type_immune_specific"])
adata = cogapsresult.T
from PyCoGAPS.analysis_functions import *
plotPatternUMAP(adata)

# pm = patternMarkers(adata, threshold="cut")
# add cell type annotations
adata.var["cell_type"] = adata_original.var["TN_assigned_cell_type_immune_broad"]

# from PyCoGAPS.parameters import *
# from PyCoGAPS.pycogaps_main import CoGAPS
pm = patternMarkers(adata, threshold="all")
# trying to get hallmark results
markers = pm["PatternMarkers"]
# colnames = list(markers)
# pattern_names = {sub for sub in colnames if sub.startswith('Pattern')}
p1_markers = list(markers["Pattern1"])