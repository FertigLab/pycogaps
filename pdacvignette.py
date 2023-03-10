if __name__ == "__main__":
    from PyCoGAPS.parameters import *
    from PyCoGAPS.pycogaps_main import CoGAPS
    import scanpy as sc

    path = "data/inputdata.h5ad"
    adata = sc.read_h5ad(path)
    adata.X = adata.X.todense()
    sc.pp.log1p(adata)
    adata = adata.T
    adata

    params = CoParams(adata=adata)

    setParams(params, {
        'nIterations':1000,
        'seed': 42,
        'nPatterns': 8,
        'useSparseOptimization': True,
        'distributed': "genome-wide"
        # 'transposeData': True
    })

    params.setDistributedParams(nSets=7)
    params.printParams()
    start = time.time()
    result = CoGAPS(adata, params)
    end = time.time()
    print("TIME:", end - start)

    result.write("data/pdacresultTEST.h5ad")