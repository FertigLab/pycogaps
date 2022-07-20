if __name__ == "__main__":
    import pickle
    from PyCoGAPS.parameters import *
    from PyCoGAPS.pycogaps_main import CoGAPS
    import scanpy as sc
    import anndata

    path = "/Users/jeanette/fertiglab/PDAC_Atlas_Pipeline/PDAC.h5ad"
    adata = sc.read_h5ad(path)
    adata.X = adata.X.todense()

    sc.pp.log1p(adata)

    params = CoParams(path)

    setParams(params, {
        'nIterations': 50000,
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