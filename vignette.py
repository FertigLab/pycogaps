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
    
    with open('data/pheno.txt', 'wb') as f:
        s3.download_fileobj('pycogaps', 'pheno.txt', f)

    path = "data/GSE98638_HCC.TCell.S5063.count.txt"

    table = pd.read_table(path)
    adata = anndata.AnnData(table.iloc[:, 2:])
    adata.obs_names = table["symbol"]
    sc.pp.log1p(adata)

    labelfile = "data/pheno.txt"
    table = pd.read_table(labelfile)
    majorCluster = table["majorCluster"]
    adata.var_names = majorCluster

    params = CoParams(adata)

    setParams(params, {
        'nIterations': 30000,
        'seed': 42,
        'nPatterns': 10,
        'useSparseOptimization': True,
        'distributed': "genome-wide",
    })

    start = time.time()
    params.setDistributedParams(nSets=10, minNS=8, maxNS=23)
    result = CoGAPS(adata, params)
    end = time.time()
    print("TIME:", end - start)
            
    print("Pickling...")
    pickle.dump(result, open("./data/30khccresult.pkl", "wb"))
    print("Pickling complete!")
    with open('./data/30khccresult.pkl', 'rb') as data:
        s3.upload_fileobj(data, 'pycogaps', '30khccresult.pkl')
