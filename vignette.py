# import pycogaps
import distributed
from PyCoGAPS import *
print("This vignette was built using pycogaps version", getVersion())

# Running CoGAPS with default parameters
# replace with the path to your data, or use this provided example
# of scRNAseq/TCRseq data from HCC patients undergoing immunotherapy

# The only required argument to CoGAPS is the path to the data.
# This can be a .csv, .tsv, .mtx, .txt, .h5, or .h5ad file containing the data.
path = "data/GSE98638_HCC.TCell.S5063.count.txt"
# run CoGAPS on your dataset
# result = CoGAPS(path)

# Running CoGAPS with custom parameters
# Most of the time we’ll want to set some parameters before running CoGAPS.
# Parameters are managed with a CoParams object. This object will store all
# parameters needed to run CoGAPS and provides a simple interface for viewing
# and setting the parameter values.
# create a CoParams object
params = CoParams(path)

# set desired parameters
setParam(params, "nPatterns", 5)
# and/or:
setParams(params, {
            'nIterations': 100,
            'seed': 42,
            'nPatterns': 5,
        })

result = CoGAPS(path, params)

# Breaking down the result object from CoGAPS
# CoGAPS returns both a C++-bound GapsResult object and an anndata result object.
# Both are equivalent, but using anndata is probably faster in most cases.
# this retrieves the anndata result object
adataresult = result["anndata"]

# Visualizing results

# Default plot: visualize how patterns vary across samples
plot(result)

# Residuals plot: calculates residuals and produces a heatmap.
plotResiduals(result)

# Plot Pattern Markers: plotPatternMarkers plots a heatmap of the
# original data clustered by the pattern markers statistic, which
# computes the most associated pattern for each gene.
plotPatternMarkers(result, legend_pos=None)

# Binary plot: binaryA creates a binarized heatmap of the A matrix
# in which the value is 1 if the value in Amean is greater than
# threshold * Asd and 0 otherwise.
binaryA(result, threshold=3)
# plotting clustered binary plot
binaryA(result, threshold=3, cluster=True)

# subset data / cluster by groups
from PyCoGAPS import *
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
labelfile = "data/pheno.txt"
table = pd.read_table(labelfile)
majorCluster = table["majorCluster"]
adata.var_names = majorCluster
params = CoParams(path)

setParams(params, {
            'nIterations': 10000,
            'seed': 42,
            'nPatterns': 10,
            'useSparseOptimization': True
        })

start = time.time()
if __name__ == '__main__':
    params.setDistributedParams(nSets=10)
    result = CoGAPS(path, params, None)
end = time.time()
print("TIME:", end - start)


result = CoGAPS(path, params)
plot(result, groups=majorCluster)

result
