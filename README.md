```
TODO: add TOC with links
```

# **1. Introduction**
Coordinated Gene Activity in Pattern Sets (CoGAPS) implements a Bayesian MCMC matrix factorization algorithm, GAPS, and links it to gene set statistic methods to infer biological process activity. It can be used to perform sparse matrix factorization on any data, and when this data represents biomolecules, to do gene set analysis.

This package presents a unified python interface, with a parallel, efficient underlying implementation in C++.


# **2. Installation**
```
git clone https://github.com/FertigLab/pycogaps --recursive
cd pycogaps
pip3 install . 
OR
python3 setup.py install
```
# **3. Usage**

To import the python CoGAPS package:
```python
from PyCoGAPS import *
```

## 3.1 Running CoGAPS with Default Parameters
The only required argument to CoGAPS is the path to the data. This can be a *.csv, .tsv, .mtx, .h5, or .h5ad* file containing the data.

```python 
# replace with the path to your data, or use this provided example
path = "./data/GIST.csv" 

# run CoGAPS on your dataset
result = CoGAPS(path)
```
```
This is pycogaps version  0.0.1
Running Standard CoGAPS on GIST.csv (1363 genes and 9 samples) with parameters: 

-- Standard Parameters --
nPatterns:  3
nIterations:  1000
seed:  0
sparseOptimization:  False

-- Sparsity Parameters --
alpha: 0.01
maxGibbsMass:  100.0

GapsResult result object with 1363 features and 9 samples
3 patterns were learned
```

<details>
  <summary> About CoGAPS print status messages </summary>

While CoGAPS is running it periodically prints status messages. For example, 20000 of 25000, Atoms: 2932(80), ChiSq: 9728, time: 00:00:29 / 00:01:19. This message tells us that CoGAPS is at iteration 20000 out of 25000 for this phase, and that 29 seconds out of an estimated 1 minute 19 seconds have passed. It also tells us the size of the atomic domain which is a core component of the algorithm but can be ignored for now. Finally, the ChiSq value tells us how closely the A and P matrices reconstruct the original data. In general, we want this value to go down - but it is not a perfect measurment of how well CoGAPS is finding the biological processes contained in the data. CoGAPS also prints a message indicating which phase is currently happening. There are two phases to the algorithm - Equilibration and Sampling.

</details>

## 3.2 Running CoGAPS with Custom Parameters

Most of the time we’ll want to set some parameters before running CoGAPS. Parameters are managed with a CoParams object. This object will store all parameters needed to run CoGAPS and provides a simple interface for viewing and setting the parameter values.

```python
# create a CoParams object
params = CoParams(path)

# set desired parameters
setParam(params, "nPatterns", 7) 
# and/or:
setParams(params, {
            'nIterations': 1500,
            'seed': 42,
            'nPatterns': 7,
        })

result = CoGAPS(path, params)
```

```
This is pycogaps version  0.0.1
Running Standard CoGAPS on GIST.csv (1363 genes and 9 samples) with parameters: 

-- Standard Parameters --
nPatterns:  7
nIterations:  1500
seed:  0
sparseOptimization:  False

-- Sparsity Parameters --
alpha: 0.01
maxGibbsMass:  100.0

GapsResult result object with 1363 features and 9 samples
7 patterns were learned
```

<details>
  <summary> About all parameters </summary>

```
TODO: list all main params here, add all params as fields in CoParams object

TODO: remove following heading and wording, shorten to say additional run params can be passed directly as arguments - list the args.
```

### 3.2.2 Run Configuration Options
The CoParams class manages the model parameters - i.e. the parameters that affect the result. There are also a few parameters that are passed directly to CoGAPS that control things like displaying the status of the run.

```python
# run config arguments can be passed as an argument to CoGAPS
result = CoGAPS(path, params, nIterations=1000, outputFrequency=250)
```

</details>



## 3.3 Breaking Down the Return Object from CoGAPS
CoGAPS returns a dictionary of the result as two representations: an `anndata` object and `GapsResult` object. For simplicity and relevancy, we will only consider the `anndata` object. CoGAPS stores the lower dimensional representation of the samples (P matrix) in the `.var` slot and the weight of the features (A matrix) in the `.obs` slot. The standard deviation across sample points for each matrix are stored in the `.uns` slots.

```python
# this retrieves the anndata result object
result["anndata"]
```
```
TODO: add image representation of anndata object 
```

## 3.4 Visualizing Output
The result object can be passed on to the analysis and plotting functions provided in the package. By default, the `plot` function displays how the patterns vary across the samples.

```python
# plot result object returned from CoGAPS
plot(result)
```

```
TODO: add image example - res_plot
```

```
TODO: output plotting, stat functions, etc. here
```

## 3.5 Running CoGAPS in Parallel
Non-Negative Matrix Factorization algorithms typically require long computation times and CoGAPS is no exception. In order to scale CoGAPS up to the size of data sets seen in practice we need to take advantage of modern hardware and parallelize the algorithm.

### 3.5.1 Multi-Threaded Parallelization
The simplest way to run CoGAPS in parallel is to provide the `nThreads` argument to CoGAPS. This allows the underlying algorithm to run on multiple threads and has no effect on the mathematics of the algorithm i.e. this is still standard CoGAPS. The precise number of threads to use depends on many things like hardware and data size. The best approach is to play around with different values and see how it effects the estimated time.

```python
result = CoGAPS(path, nThreads=1, nIterations=10000, seed=5)
```

Note this method relies on CoGAPS being compiled with OpenMP support, use `buildReport` to check.
```python
print(getBuildReport())
```

### 3.5.2 Distributed CoGAPS
For large datasets (greater than a few thousand genes or samples) the multi-threaded parallelization isn’t enough. It is more efficient to break up the data into subsets and perform CoGAPS on each subset in parallel, stitching the results back together at the end (Stein-O’Brien et al. (2017)).

In order to use these extensions, some additional parameters are required. We first need to set CoParam's `distributed` parameter to be `genome-wide` using `setParam`. Next, `nSets` specifies the number of subsets to break the data set into. `cut`, `minNS`, and `maxNS` control the process of matching patterns across subsets and in general should not be changed from defaults. More information about these parameters can be found in the original papers. These parameters need to be set with a different function, `setDistributedParameters`, than `setParam` since they depend on each other. Here we only set `nSets` (always required), but we have the option to pass the other parameters as well.

```python
setParam(params, "distributed", "genome-wide")
params.setDistributedParameters(nSets=3)
```

Setting `nSets` requires balancing available hardware and run time against the size of your data. In general, `nSets` should be less than or equal to the number of nodes/cores that are available. If that is true, then the more subsets you create, the faster CoGAPS will run - however, some robustness can be lost when the subsets get too small. The general rule of thumb is to set `nSets` so that each subset has between 1000 and 5000 genes or cells. We will see an example of this on real data in the next two sections.

Once the distributed parameters have been set we can call CoGAPS as usual.

```python
result = CoGAPS(path, params)
```

# **4. Additional Features of CoGAPS**

## 4.1 Checkpoint System -- Saving/Loading CoGAPS Runs
CoGAPS allows the user to save their progress throughout the run, and restart from the latest saved “checkpoint”. This is intended so that if the server crashes in the middle of a long run it doesn’t need to be restarted from the beginning. Set the `checkpointInterval` parameter to save checkpoints and pass a file name as `checkpointInFile` to load from a checkpoint.

## 4.2 Transposing Data
If your data is stored as samples x genes, CoGAPS allows you to pass `transposeData=True` and will automatically read the transpose of your data to get the required genes x samples configuration.

## 4.3 Passing Uncertainty Matrix
In addition to providing the data, the user can also specify an uncertainty measurement - the standard deviation of each entry in the data matrix. By default, CoGAPS assumes that the standard deviation matrix is 10% of the data matrix. This is a reasonable heuristic to use, but for specific types of data you may be able to provide better information.

## 4.4 Distributed CoGAPS

### 4.4.1 Methods of Subsetting Data
The default method for subsetting the data is to uniformly break up the rows (cols) of the data. There is an alternative option where the user provides an annotation vector for the rownames (colnames) of the data and gives a weight to each category in the annotation vector. Equal sized subsets are then drawn by sampling all rows (cols) according to the weight of each category.

```
```
Finally, the user can set `explicitSets` which is a list of character or numeric vectors indicating which names or indices of the data should be put into each set. Make sure to set nSets to the correct value before passing `explicitSets`.

### 4.4.2 Additional Return Information
When running GWCoGAPS or scCoGAPS, some additional metadata is returned that relates to the pattern matching process. This process is how CoGAPS stitches the results from each subset back together.

```
```

### 4.4.3 Manual Pipeline
CoGAPS allows for a custom process for matching the patterns together. If you have a result object from a previous run of Distributed CoGAPS, the unmatched patterns for each subset are found by calling `getUnmatchedPatterns`. Apply any method you like as long as the result is a matrix with the number of rows equal to the number of samples (genes) and the number of columns is equal to the number of patterns. Then pass the matrix to the `fixedPatterns` argument along with the original parameters for the GWCoGAPS/scCoGAPS run.

```
```

# **5. Citing CoGAPS**
If you use the CoGAPS package for your analysis, please cite Fertig et al. (2010)

If you use the gene set statistic, please cite Ochs et al. (2009)

# **References**
Fertig, Elana J., Jie Ding, Alexander V. Favorov, Giovanni Parmigiani, and Michael F. Ochs. 2010. “CoGAPS: An R/C++ Package to Identify Patterns and Biological Process Activity in Transcriptomic Data.” Bioinformatics 26 (21): 2792–3. https://doi.org/10.1093/bioinformatics/btq503.

Ochs, Michael F., Lori Rink, Chi Tarn, Sarah Mburu, Takahiro Taguchi, Burton Eisenberg, and Andrew K. Godwin. 2009. “Detection of Treatment-Induced Changes in Signaling Pathways in Gastrointestinal Stromal Tumors Using Transcriptomic Data.” Cancer Research 69 (23): 9125–32. https://doi.org/10.1158/0008-5472.CAN-09-1709.

Seung, Sebastian, and Daniel D. Lee. 1999. “Learning the Parts of Objects by Non-Negative Matrix Factorization.” Nature 401 (6755): 788–91. https://doi.org/10.1038/44565.

Stein-O’Brien, Genevieve L., Jacob L. Carey, Wai S. Lee, Michael Considine, Alexander V. Favorov, Emily Flam, Theresa Guo, et al. 2017. “PatternMarkers & Gwcogaps for Novel Data-Driven Biomarkers via Whole Transcriptome Nmf.” Bioinformatics 33 (12): 1892–4. https://doi.org/10.1093/bioinformatics/btx058.