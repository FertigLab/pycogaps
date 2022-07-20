if __name__ == "__main__":
    import pickle
    from PyCoGAPS.parameters import *
    from PyCoGAPS.pycogaps_main import CoGAPS
    import scanpy as sc

    # load CoGAPS result object
    path = "/Users/jeanette/fertiglab/PDAC_Atlas_Pipeline/PDAC.h5ad"
    adata = sc.read_h5ad(path)

    adata.X = adata.X.todense()
    # new_matrix = np.genfromtxt("/Users/jeanette/fertiglab/pycogaps/data/densematR.csv", delimiter=",")
    # counts = scipy.io.mmread("/Users/jeanette/fertiglab/PDAC_Atlas_Pipeline/CoGAPS/counts.mtx")
    # countsarray = counts.toarray()
    # adata.X = np.transpose(countsarray)
    # adata = adata.T
    # scale data matrix (beginning single cell workflow)
    # sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    params = CoParams(path)

    setParams(params, {
        'nIterations': 1000,
        'seed': 42,
        'nPatterns': 8,
        'useSparseOptimization': True,
        'distributed': "genome-wide",
    })

    params.setDistributedParams(nSets=15, minNS=8, maxNS=23, cut=8)
    params.printParams()
    start = time.time()
    result = CoGAPS(adata, params)
    end = time.time()
    print("TIME:", end - start)

    print("Pickling...")
    pickle.dump(result, open("./data/testing032522.pkl", "wb"))
    print("Pickling complete!")







    # # for making an embedding, only need the most variable genes
    # sc.pp.highly_variable_genes(adata)
    # adata = adata[:, adata.var.highly_variable]
    #
    # # scale data and compute PCA
    # sc.pp.scale(adata)
    # sc.tl.pca(adata, svd_solver='arpack')
    # sc.pl.pca_variance_ratio(adata, log=True)
    #
    # # find neighbor embeddings and run UMAP
    # sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
    # sc.tl.umap(adata)
    #
    # # add categorical annotations to observation matrix for umap plot
    # adata.obs['majorCluster']=list(majorCluster)
    #
    # # plot pattern amplitude on UMAP
    # sc.pl.umap(adata, color=["majorCluster"])
