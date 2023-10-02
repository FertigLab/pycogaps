# %%
import anndata
import pandas as pd
import scanpy as sc

from PyCoGAPS.analysis_functions import *
from PyCoGAPS.parameters import *
from PyCoGAPS.pycogaps_main import CoGAPS

# %%

# %%
cogapsresult = anndata.read_h5ad("data/visiumresult50k.h5ad")
adata = cogapsresult

# %%
import numpy as np

path = "data/VI_116_4"
adata = sc.read_visium(path)

adata.obs["Pattern1"] = cogapsresult.var["Pattern1"]
adata.obs["Pattern2"] = cogapsresult.var["Pattern2"]
adata.obs["Pattern3"] = cogapsresult.var["Pattern3"]
adata.obs["Pattern4"] = cogapsresult.var["Pattern4"]
adata.obs["Pattern5"] = cogapsresult.var["Pattern5"]

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat")

sc.pp.pca(adata)
sc.pl.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters", resolution=0.5)

import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(
    adata,
    color=["Pattern1", "Pattern2", "Pattern3", "Pattern4", "Pattern5"],
)
