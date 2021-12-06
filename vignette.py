print("checkpoint 1")
if __name__ == '__main__':
    print("checkpoint 2")
    from PyCoGAPS import *
    import pickle
    import scanpy as sc
    print("checkpoint 3")
    import boto3
    print("This vignette was built using pycogaps version", getVersion())
    
    s3 = boto3.client('s3')
    with open('data/GSE98638_HCC.TCell.S5063.count.txt', 'wb') as f:
        s3.download_fileobj('pycogaps', 'GSE98638_HCC.TCell.S5063.count.txt', f)

    path = "data/GSE98638_HCC.TCell.S5063.count.txt"

    table = pd.read_table(path)
    adata = anndata.AnnData(table.iloc[:, 2:])
    adata.obs_names = table["symbol"]
    sc.pp.log1p(adata)

    labelfile = "data/pheno.txt"
    table = pd.read_table(labelfile)
    majorCluster = table["majorCluster"]
    adata.var_names = majorCluster

    params = CoParams(path)

    setParams(params, {
        'nIterations': 50000,
        'seed': 42,
        'nPatterns': 10,
        'useSparseOptimization': True,
        'distributed': "genome-wide",
    })

    start = time.time()
    params.setDistributedParams()
    result = CoGAPS(path, params)
    end = time.time()
    print("TIME:", end - start)
            
    print("Pickling...")
    pickle.dump(result, open("./data/50kresultfinal.pkl", "wb"))
    print("Pickling complete!")
    with open('./data/50kresult.pkl', 'rb') as data:
        s3.upload_fileobj(data, 'pycogaps', '50kresult.pkl')
