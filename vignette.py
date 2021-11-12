from PyCoGAPS import *
import pickle
import scanpy as sc
print("This vignette was built using pycogaps version", getVersion())

# Running CoGAPS with default parameters
# replace with the path to your data, or use this provided example
# of scRNAseq/TCRseq data from HCC patients undergoing immunotherapy

# The only required argument to CoGAPS is the path to the data.
# This can be a .csv, .tsv, .mtx, .txt, .h5, or .h5ad file containing the data.
path = "data/GSE98638_HCC.TCell.S5063.count.txt"
table = pd.read_table(path)
adata = anndata.AnnData(table.iloc[:, 2:])
adata.obs_names = table["symbol"]
sc.pp.log1p(adata)
labelfile = "data/pheno.txt"
table = pd.read_table(labelfile)
majorCluster = table["majorCluster"]
# adata.var_names = majorCluster
params = CoParams(path)

setParams(params, {
            'nIterations': 10,
            'seed': 42,
            'nPatterns': 10,
            'useSparseOptimization': True,
            'distributed': 'genome-wide'
        })

if __name__ == '__main__':
    start = time.time()
    # params.setDistributedParams()
    result = CoGAPS(adata, params)
    end = time.time()
    print("TIME:", end - start)
    print("Pickling...")
    pickle.dump(result, open("./data/200result.pkl", "wb"))
    print("Pickling complete!")
    plot(result, groups=majorCluster)

# unpickled = pickle.load(open("./data/testresult.pkl", "rb"))