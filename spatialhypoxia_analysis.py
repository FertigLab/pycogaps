import pickle
import scanpy as sc

with open('./data/spatialhypoxia4patterns.pkl', 'rb') as fp:
    hypoxiaresult4pattern = pickle.load(fp)



adata = hypoxiaresult4pattern.T
sc.pp.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable]

# scale data and compute PCA
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

# find neighbor embeddings and run UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

patterns = list(adata.obs.columns)
# plot pattern amplitude on UMAP
sc.pl.umap(adata, color=patterns)
