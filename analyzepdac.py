import anndata
import scanpy as sc
import pandas as pd

pycogapsresult = anndata.read_h5ad("data/pdacresult.h5ad")

coldata = pd.read_csv("data/PDACcoldata.csv")
# Index(['Unnamed: 0', 'barcode_raw', 'celltype', 'sample_ID',
#        'sample_ID_celltype', 'TN', 'TN_manuscript', 'manuscript', 'nCount_RNA',
#        'nFeature_RNA', 'percent.mt', 'Size_Factor', 'TN_cluster_resolution_5',
#        'TN_assigned_cell_type', 'TN_assigned_cell_type_immune',
#        'TN_assigned_cell_type_immune_specific',
#        'TN_assigned_cell_type_immune_broad', 'cc', 'ccstage',
#        'Classifier_T_duct', 'Classifier_T_Fibroblast_only',
#        'Classifier_T_Fibroblast_Stellate'],
#       dtype='object')

# adata = pycogapsresult.T
adata = pycogapsresult
adata.obs["cell type"] = list(coldata["TN_assigned_cell_type_immune_specific"])
adata = pycogapsresult.T
from PyCoGAPS.analysis_functions import *
plotPatternUMAP(adata)