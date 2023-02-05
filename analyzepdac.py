import anndata
import pandas as pd
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
adata_original = sc.read_h5ad("data/PDAC_Peng_Epi.h5ad").T

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


# gsea = pd.DataFrame(enr.results)
#
# gsea = gsea[gsea["Adjusted P-value"] < 0.05]

import seaborn as sns



# sns.barplot(data= gsea, x="Term", y="Adjusted P-value")

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
p7gsea = p7gsea[p7gsea["P-value"] < 0.05]
sns.barplot(data= p7gsea, x="P-value", y="Term")

p8_markers = markers["Pattern8"]


prerank = adata.obs["Pattern7"]



marker_ranks = pm["PatternMarkerRanks"]
p7_markers_prerank = pd.DataFrame(markers["Pattern7"])

p7_ranks = marker_ranks["Pattern7"]

p7prerank = gp.prerank(rnk=p7_scores, # or rnk = rnk,
                     gene_sets='MSigDB_Hallmark_2020',
                     threads=4,
                     min_size=5,
                     max_size=1000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=None, # don't write to disk
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )

p7prerank_df = pd.DataFrame(p7prerank.results)
plotdf = p7prerank_df.T
plotdf = plotdf[plotdf["pval"] < 0.05]

sns.barplot(data= plotdf, x="neg.log.q", y="Term")