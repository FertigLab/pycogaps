if __name__ == '__main__':
    from PyCoGAPS import *
    import pickle
    from PyCoGAPS.parameters import *
    from PyCoGAPS.pycogaps_main import CoGAPS
    print("This vignette was built using pycogaps version", getVersion())
    path = "data/GSE98638_HCC.TCell.S5063.count.txt"
    pd_table = pd.read_table(path)

    # remove any rows and columns that do not contain data, and then create an anndata object
    table = pd.DataFrame(data=pd_table.values, index=pd_table.index, columns=pd_table.columns)
    adata = anndata.AnnData(table.iloc[:, 2:])  # exclude annotation columns
    adata.obs_names = table["symbol"]

    labelfile = "data/pheno.txt"
    annotation_table = pd.read_table(labelfile)
    majorCluster = annotation_table["majorCluster"]
    patientID = annotation_table["Patient"]
    sampleType = annotation_table["sampleType"]

    # path = "data/GSE98638_HCC.TCell.S5063.count.txt"
    params = CoParams(path)

    setParams(params, {
        'nIterations': 100,
        'seed': 42,
        'nPatterns': 10,
        'useSparseOptimization': True,


        'distributed': "single-cell",
    })

    # params.setDistributedParams(nSets=10, minNS=6, maxNS=23)
    params.setDistributedParams()
    start = time.time()
    result = CoGAPS(adata, params)
    end = time.time()
    print("TIME:", end - start)

    print("Pickling...")
    pickle.dump(result, open("./data/030422result.pkl", "wb"))
    print("Pickling complete!")

    with open('./data/030422result.pkl', 'rb') as fp:
        result = pickle.load(fp)

    result.uns["sampleType"] = sampleType
    result.uns["patientID"] = patientID
    result.uns["majorCluster"] = majorCluster

    import scanpy as sc

    # the scanpy package requires genes to be in columns, and cells in rows
    adata = result.transpose()

    sc.pp.log1p(adata)
    sc.tl.pca(adata, svd_solver='arpack')

    sc.pl.highest_expr_genes(adata, n_top=20, )






    patterns = list(result.obs.columns)
    # verbosity: errors (0), warnings (1), info (2), hints (3)
    # set logging options
    sc.settings.verbosity = 3
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')

    sc.pp.highly_variable_genes(result, min_mean=0.0125, max_mean=3, min_disp=0.5)

    result = result[:, result.var.highly_variable]
    sc.pp.scale(result, max_value=10)
    sc.tl.pca(result, svd_solver='arpack')
    sc.pp.neighbors(result)
    sc.tl.umap(result)
    # plot pattern amplitude on UMAP
    sc.pl.umap(adata, color=patterns)
