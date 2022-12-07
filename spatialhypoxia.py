if __name__ == "__main__":
    import pickle
    from PyCoGAPS.parameters import *
    from PyCoGAPS.pycogaps_main import CoGAPS
    import scanpy as sc

    # load CoGAPS result object
    path = "/Users/jeanette/fertiglab/spatialhypoxia/outs/filtered_feature_bc_matrix/matrix.mtx"
    rawdata = sc.read_mtx(path)
    rawdata.X=rawdata.X.todense()
    features = pd.read_csv("/Users/jeanette/fertiglab/spatialhypoxia/outs/filtered_feature_bc_matrix/features.tsv", sep="\t", header=None)
    cell_labels = pd.read_csv("/Users/jeanette/fertiglab/spatialhypoxia/outs/filtered_feature_bc_matrix/barcodes.tsv",
                           sep="\t", header=None)
    rawdata.var_names=cell_labels[0]
    rawdata.obs_names=features[1]

    sc.pp.log1p(rawdata)

    params = CoParams(path=path)

    setParams(params, {
        'nIterations': 1000,
        'seed': 42,
        'nPatterns': 4,
        'useSparseOptimization': True,
        'distributed': "genome-wide"
    })

    params.setDistributedParams(nSets=8)
    params.printParams()
    start = time.time()
    result = CoGAPS(rawdata, params)
    end = time.time()
    print("TIME:", end - start)

    print("Pickling...")
    pickle.dump(result, open("./data/spatialhypoxia4patterns.pkl", "wb"))
    print("Pickling complete!")

