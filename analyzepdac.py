import anndata
import scanpy as sc
# import pandas as pd

pycogapsresult = anndata.read_h5ad("data/pdacresult.h5ad")

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
adata = pycogapsresult
# get readable gene names from original object
adata_original = sc.read_h5ad("/Users/jeanette/fertiglab/PDAC_Atlas_Pipeline/PDAC_Peng_Epi.h5ad").T
adata.obs_names = adata_original.obs["gene_short_name"]

# adata = pycogapsresult.T

# adata.obs["cell type"] = list(coldata["TN_assigned_cell_type_immune_specific"])
# adata = pycogapsresult.T
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
for key in markers:
    print(markers[key])
    thispattern_markers = markers[key]




import gseapy as gp
# run enrichr
# if you are only intrested in dataframe that enrichr returned, please set outdir=None
enr = gp.enrichr(gene_list=p1_markers, # or "./tests/data/gene_list.txt",
                 gene_sets=['MSigDB_Hallmark_2020'],
                 organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )

gsea = pd.DataFrame(enr.results)

gsea = gsea[gsea["Adjusted P-value"] < 0.05]

import seaborn as sns



sns.barplot(data= gsea, x="Term", y="Adjusted P-value")

p2_markers = markers["Pattern2"]
p3_markers = markers["Pattern3"]
p4_markers = markers["Pattern4"]
p5_markers = markers["Pattern5"]
p6_markers = markers["Pattern6"]
p7_markers = list(markers["Pattern7"])
p7enr = gp.enrichr(gene_list=p7_markers, # or "./tests/data/gene_list.txt",
                 gene_sets=['MSigDB_Hallmark_2020'],
                 organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )
p7gsea = pd.DataFrame(p7enr.results)
p7gsea = p7gsea[p7gsea["P-value"] < 0.2]
sns.barplot(data= p7gsea, x="Term", y="P-value")

p8_markers = markers["Pattern8"]

