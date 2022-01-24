# **PyCoGAPS**

Coordinated Gene Activity in Pattern Sets (CoGAPS) implements a Bayesian MCMC matrix factorization algorithm, GAPS, and links it to gene set statistic methods to infer biological process activity. It can be used to perform sparse matrix factorization on any data, and when this data represents biomolecules, to do gene set analysis.

This package, PyCoGAPS, presents a unified Python interface, with a parallel, efficient underlying implementation in C++.

# **Table of Contents**

1. [ Usage ](#1-usage)  
  1.1 [ Running CoGAPS with Default Parameters ](#11-running-cogaps-with-default-parameters)  
  1.2 [ Running CoGAPS with Custom Parameters ](#12-running-cogaps-with-custom-parameters)  
  1.3 [ Running CoGAPS in Parallel ](#13-running-cogaps-in-parallel)  
2. [ Analyzing the Result ](#2-analyzing-the-result)  
  2.3 [ Breaking Down the Return Object from CoGAPS ](#23-breaking-down-the-return-object-from-cogaps)  
  2.4 [ Visualizing Output ](#24-visualizing-output)        
3. [ Additional Features of CoGAPS ](#3-additional-features-of-cogaps)  
  3.1 [ Checkpoint System: Saving/Loading CoGAPS Runs ](#31-checkpoint-system-savingloading-cogaps-runs)  
  3.2 [ Transposing Data ](#32-transposing-data)  
  3.3 [ Passing Uncertainty Matrix ](#33-passing-uncertainty-matrix)  
  3.4 [ Distributed CoGAPS ](#34-distributed-cogaps)  
4. [ Citing CoGAPS ](#4-citing-cogaps)


# **1. Usage**
Please follow the steps below to run the PyCoGAPS vignette:
1. Install docker at https://docs.docker.com/desktop/mac/install/ 
2. Open docker
3. Copy the commands and paste in terminal (Tested via Mac OX)

```
docker pull ashleyt2000/pycogaps:docker_pycogaps
mkdir PyCoGAPS
cd PyCoGAPS
curl -O https://raw.githubusercontent.com/FertigLab/pycogaps/master/params.yaml
mkdir data
cd data
curl -O https://raw.githubusercontent.com/FertigLab/pycogaps/master/data/GIST.csv
cd ..
docker run -v $PWD:$PWD ashleyt2000/pycogaps:docker_pycogaps $PWD/params.yaml

```
This produces a CoGAPS run on a simple dataset with default parameters.
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

</br>

While CoGAPS is running it periodically prints status messages. For example, `20000 of 25000, Atoms: 2932(80), ChiSq: 9728, time: 00:00:29 / 00:01:19`. This message tells us that CoGAPS is at iteration 20000 out of 25000 for this phase, and that 29 seconds out of an estimated 1 minute 19 seconds have passed. It also tells us the size of the atomic domain which is a core component of the algorithm but can be ignored for now. Finally, the `ChiSq` value tells us how closely the A and P matrices reconstruct the original data. In general, we want this value to go down - but it is not a perfect measurment of how well CoGAPS is finding the biological processes contained in the data. CoGAPS also prints a message indicating which phase is currently happening. There are two phases to the algorithm - Equilibration and Sampling.

</details>

It also sets up your working directory to be in the `PyCoGAPS` folder with the following structure and files:
```
PyCoGAPS
├── data
│   └── GIST.csv
├── params.yaml
├── output
│   └── result.pkl
```

## 1.2 Running CoGAPS with Custom Parameters

In order to analyze your desired data, we'll need to input it and modify the default parameters before running CoGAPS. All parameter values can be modified directly in the `params.yaml` file already downloaded earlier. 

Please follow the steps below to run PyCoGAPS with custom parameters:
1. Open `params.yaml` with any text or code editor
2. Modify the desired parameters
3. Save file
4. Copy the command and paste in terminal:
```
docker run -v $PWD:$PWD ashleyt2000/pycogaps:docker_pycogaps $PWD/params.yaml
```

A snippet of `params.yaml` is shown below, where we have changed some default parameter values to our own specified example values.

```
## This file holds all parameters to be passed into PyCoGAPS.
## To modify default parameters, simply replace parameter values below with user-specified values, and save file. 

# RELATIVE path to data -- make sure to move your data into the created data/ folder
path: data/liver_dataset.txt

# result output file name (output saved in the created output/ folder as a .pkl file)
result_file: result.pkl

standard_params:
  # number of patterns CoGAPS will learn
  nPatterns: 10
  # number of iterations for each phase of the algorithm
  nIterations: 5000
  # random number generator seed
  seed: 0
  # speeds up performance with sparse data (roughly >80% of data is zero), note this can only be used with the default uncertainty
  useSparseOptimization: True
 
... 

```

<details>
  <summary> Standard & Sparsity Parameters </summary>

- Standard Parameters
  - `nPatterns` number of patterns CoGAPS will learn
  - `nIterations` number of iterations for each phase of the algorithm
  - `seed` random number generator seed
  - `sparseOptimization` speeds up performance with sparse data (roughly >80% of data is zero), note this can only be used with the default uncertainty

- Sparsity Parameters
  - `alphaA` sparsity parameter for feature matrix
  - `alphaP` sparsity parameter for sample matrix
  - `maxGibbsMassA` atomic mass restriction for feature matrix
  - `maxGibbsMassP` atomic mass restriction for sample matrix

</details>

<details>
  <summary> Run Configuration Parameters </summary>

- Run Configuration Parameters (these can be passed in to CoGAPS directly)
  - `nThreads` maximum number of threads to run on
  - `messages` T/F for displaying output
  - `outputFrequency` number of iterations between each output (set to 0 to disable status updates)
  - `uncertainty` uncertainty matrix - either a matrix or a supported file type
  - `checkpointOutFile` name of the checkpoint file to create
  - `checkpointInterval` number of iterations between each checkpoint (set to 0 to disable checkpoints)
  - `checkpointInFile` if this is provided, CoGAPS runs from the checkpoint contained in this file
  - `transposeData` T/F for transposing data while reading it in - useful for data that is stored as samples x genes since CoGAPS requires data to be genes x samples
  - `workerID` if calling CoGAPS in parallel the worker ID can be specified
  - `asynchronousUpdates` enable asynchronous updating which allows for multi-threaded runs
  - `nSnapshots` how many snapshots to take in each phase, setting this to 0 disables snapshots
  - `snapshotPhase` which phase to take snapsjots in e.g. "equilibration", "sampling", "all"

</details>

<details>
  <summary> Distributed Parameters </summary>

- Distributed Parameters
  - `distributed`  must set to `"genome-wide"`
  - `nSets` [distributed parameter] number of sets to break data into
  - `cut` [distributed parameter] number of branches at which to cut dendrogram used in pattern matching
  - `minNS` [distributed parameter] minimum of individual set contributions a cluster must contain
  - `maxNS` [distributed parameter] maximum of individual set contributions a cluster can contain
  - `explicitSets` [distributed parameter] specify subsets by index or name
  - `samplingAnnotation` [distributed parameter] specify categories along the rows (cols) to use for weighted sampling
  - `samplingWeight` [distributed parameter] weights associated with  samplingAnnotation

</details>

<details>
  <summary> Additional Parameters </summary>

- Additional Parameters
  - `subsetIndices` set of indices to use from the data
  - `subsetDim` which dimension (0=rows, 1=cols) to subset
  - `geneNames` vector of names of genes in data
  - `sampleNames` vector of names of samples in data
  - `fixedPatterns` fix either 'A' or 'P' matrix to these values, in the context of distributed CoGAPS, the first phase is skipped and `fixedPatterns` is used for all sets allowing manual pattern matching, as well as fixed runs of standard CoGAPS
  - `whichMatrixFixed` either 'A' or 'P', indicating which matrix is fixed
  - `takePumpSamples` whether or not to take PUMP samples

</details>


## 2.3 Running CoGAPS in Parallel
Non-Negative Matrix Factorization algorithms typically require long computation times and CoGAPS is no exception. In order to scale CoGAPS up to the size of data sets seen in practice we need to take advantage of modern hardware and parallelize the algorithm.

### 2.3.1 Multi-Threaded Parallelization
The simplest way to run CoGAPS in parallel is to modify the `nThreads` parameter in `params.yaml`. This allows the underlying algorithm to run on multiple threads and has no effect on the mathematics of the algorithm i.e. this is still standard CoGAPS. The precise number of threads to use depends on many things like hardware and data size. The best approach is to play around with different values and see how it effects the estimated time.

A snippet of `params.yaml` is shown below where `nThreads` parameter is modified.

```
## This file holds all parameters to be passed into PyCoGAPS.
...

run_params:
  # maximum number of threads to run on
  nThreads: 4
```

Note this method relies on CoGAPS being compiled with OpenMP support, use `buildReport` to check.
```python
print(getBuildReport())
```

### 2.3.2 Distributed CoGAPS
For large datasets (greater than a few thousand genes or samples) the multi-threaded parallelization isn’t enough. It is more efficient to break up the data into subsets and perform CoGAPS on each subset in parallel, stitching the results back together at the end (Stein-O’Brien et al. (2017)).

In order to use these extensions, some additional parameters are required, specifically modifying the `distributed_params` in `params.yaml`. We first need to set`distributed`  to be `genome-wide.` Next, `nSets` specifies the number of subsets to break the data set into. `cut`, `minNS`, and `maxNS` control the process of matching patterns across subsets and in general should not be changed from defaults. More information about these parameters can be found in the original papers. 

A snippet of `params.yaml` is shown below where `distributed_params` parameters are modified.

```
## This file holds all parameters to be passed into PyCoGAPS.
...

distributed_params:
  #  either null or genome-wide
  distributed: genome-wide
  # number of sets to break data into
  nSets: 4
  # number of branches at which to cut dendrogram used in pattern matching
  cut: null
  # minimum of individual set contributions a cluster must contain
  minNS: null
  # maximum of individual set contributions a cluster can contain
  maxNS: null
```

Setting `nSets` requires balancing available hardware and run time against the size of your data. In general, `nSets` should be less than or equal to the number of nodes/cores that are available. If that is true, then the more subsets you create, the faster CoGAPS will run - however, some robustness can be lost when the subsets get too small. The general rule of thumb is to set `nSets` so that each subset has between 1000 and 5000 genes or cells. 


# **2. Analyzing the Result**

## 2.1 Breaking Down the Return Object from CoGAPS
CoGAPS saves the result in a pickle file, which is a serialized Python object. It stores a dictionary of the result as two representations: an `anndata` object and `GapsResult` object. For simplicity and relevancy, we will only consider the `anndata` object. CoGAPS stores the lower dimensional representation of the samples (P matrix) in the `.var` slot and the weight of the features (A matrix) in the `.obs` slot. The standard deviation across sample points for each matrix are stored in the `.uns` slots.

```python
import pickle

# path to your result file
pkl_path = "./output/result.pkl"

# this unpickles the result object for use
result = pickle.load(open(pkl_path, "rb"))

# this retrieves the anndata result object
result["anndata"]
```

![alt text][anndata result] 

[anndata result]: https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/anndata-result.png "anndata result object"


## 2.2 Visualizing Output
The result object can be passed on to the analysis and plotting functions provided in the package. 

2.2.1 [ Default Plot ](#241-default-plot)  
2.2.2 [ Residuals Plot ](#242-residuals-plot)  
2.2.3 [ Pattern Markers Plot ](#243-pattern-markers-plot)  
2.2.4 [ Binary Plot ](#244-binarya-plot)  
2.2.5 [ Calculate CoGAPS Statistics ](#245-calculate-cogaps-statistics)  
2.2.6 [ Calculate Gene GSS Statistic ](#246-calculate-gene-gss-statistic)  
2.2.7 [ Calculate Gene GS Probability ](#247-calculate-gene-gs-probability)  



### 2.2.1 Default Plot
By default, the `plot` function displays how the patterns vary across the samples.

```python
# plot result object returned from CoGAPS
plot(result)
```

![alt text][show] 

[show]: https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/res_show.png "show result function"

### 2.2.2 Residuals Plot
`plotResiduals` calculates residuals and produces a heatmap.

```python
plotResiduals(result)
```

![alt text][plot residuals] 

[plot residuals]: https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/plot_residuals.png "plot residuals"


### 2.2.3 Pattern Markers Plot
`plotPatternMarkers` plots a heatmap of the original data clustered by the pattern markers statistic, which computes the most associated pattern for each gene.

```python
plotPatternMarkers(result, legend_pos=None)
```

![alt text][plot pm] 

[plot pm]: https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/plot_pm.png "plot pattern markers"

### 2.2.4 Binary Plot
`binaryA` creates a binarized heatmap of the A matrix in which the value is 1 if the value in Amean is greater
than `threshold * Asd` and 0 otherwise.

```python
binaryA(result, threshold=3)
```

![alt text][plot binaryA] 

[plot binaryA]: https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/binaryA.png "plot binary hm"

```python
# plotting clustered binary plot
binaryA(result, threshold=3, cluster=True)
```

![alt text][plot binaryA cluster] 

[plot binaryA cluster]: https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/binaryA_cluster.png "plot binary hm, cluster"

### 2.2.5 Calculate CoGAPS Statistics
`calcCoGAPSStat` calculates a statistic to determine if a pattern is enriched in a a particular set of measurements or samples.

```python
# sets is list of sets of measurements/samples
stats = calcCoGAPSStat(result, sets=['Hs.101174', 'Hs.1012'])
```

```
{'twoSidedPValue':               
Pattern1  0.496
Pattern2  0.353
Pattern3  0.289, 
'GSUpreg':                       
Pattern1  0.496
Pattern2  0.647
Pattern3  0.711, 
'GSDownreg':               
Pattern1  0.504
Pattern2  0.353
Pattern3  0.289, 
'GSActEst':               
Pattern1  0.008
Pattern2 -0.294
Pattern3 -0.422}
```

### 2.2.6 Calculate Gene GSS Statistic
`calcGeneGSStat` calculates the probability that a gene listed in a gene set behaves like other genes in the set within the given data set.

```python
stats = calcGeneGSStat(result, GStoGenes=['Hs.101174', 'Hs.1012'], numPerm=1000)
```

```
Hs.101174  0.422955
Hs.1012    0.391747
```

### 2.2.7 Compute Gene GS Probability

`computeGeneGSProb` computes the p-value for gene set membership using the CoGAPS-based statistics developed in Fertig et al. (2012). This statistic refines set membership for each candidate gene in a set specified in GSGenes by comparing the inferred activity of that gene to the average activity of the set.

```python
stats = computeGeneGSProb(result, GStoGenes=['Hs.101174', 'Hs.1012'])
```

```
Hs.101174  0.617193
Hs.1012    0.887583
```

# **3. Additional Features of CoGAPS**

## 3.1 Checkpoint System: Saving/Loading CoGAPS Runs
CoGAPS allows the user to save their progress throughout the run, and restart from the latest saved “checkpoint”. This is intended so that if the server crashes in the middle of a long run it doesn’t need to be restarted from the beginning. Set the `checkpointInterval` parameter to save checkpoints and pass a file name as `checkpointInFile` to load from a checkpoint.

A snippet of `params.yaml` is shown where we enable the checkpoint system, saving CoGAPS run.
```
## This file holds all parameters to be passed into PyCoGAPS.
...

run_params:
  checkpointOutFile: gaps_checkpoint.out
  # number of iterations between each checkpoint (set to 0 to disable checkpoints)
  checkpointInterval: 250
  # if this is provided, CoGAPS runs from the checkpoint contained in this file
  checkpointInFile: null
```

A snippet of `params.yaml` is shown where we now load the saved CoGAPS checkpoint file to continue the run.
```
## This file holds all parameters to be passed into PyCoGAPS.
...

run_params:
  checkpointOutFile: null
  # number of iterations between each checkpoint (set to 0 to disable checkpoints)
  checkpointInterval: null
  # if this is provided, CoGAPS runs from the checkpoint contained in this file
  checkpointInFile: gaps_checkpoint.out
```

## 3.2 Transposing Data
If your data is stored as samples x genes, CoGAPS allows you to pass `transposeData: True` and will automatically read the transpose of your data to get the required genes x samples configuration. 

A snippet of `params.yaml` is shown where we now load the saved CoGAPS checkpoint file to continue the run.
```
## This file holds all parameters to be passed into PyCoGAPS.
...

run_params:
  # T/F for transposing data while reading it in - useful for data that is stored as samples x genes since CoGAPS requires data to be genes x samples
  transposeData: True
```

## 3.3 Passing Uncertainty Matrix
In addition to providing the data, the user can also specify an uncertainty measurement - the standard deviation of each entry in the data matrix. By default, CoGAPS assumes that the standard deviation matrix is 10% of the data matrix. This is a reasonable heuristic to use, but for specific types of data you may be able to provide better information. Make sure to save your uncertainty file into the `data/` file.

A snippet of `params.yaml` is shown where we now load the saved CoGAPS checkpoint file to continue the run.
```
## This file holds all parameters to be passed into PyCoGAPS.
...

run_params:
  # uncertainty matrix - either a matrix or a supported file type
  uncertainty: data/GIST_uncertainty.csv
```

## 3.4 Distributed CoGAPS

3.4.1 [ Methods of Subsetting Data ](#341-methods-of-subsetting-data)    

### 3.4.1 Methods of Subsetting Data
The default method for subsetting the data is to uniformly break up the rows (cols) of the data. There is an alternative option where the user provides an annotation vector for the rownames (colnames) of the data and gives a weight to each category in the annotation vector. Equal sized subsets are then drawn by sampling all rows (cols) according to the weight of each category.

A snippet of `params.yaml` is shown below where `distributed_params` parameters are modified to subset the data.
```
## This file holds all parameters to be passed into PyCoGAPS.
...

distributed_params:
  # specify categories along the rows (cols) to use for weighted sampling
  samplingAnnotation: ['IM00', 'IM02', 'IM00']
  # weights associated with  samplingAnnotation
  samplingWeight: {'IM00': 2, 'IM02': 0.5}
```


Finally, the user can set `explicitSets` which is a list of character or numeric vectors indicating which names or indices of the data should be put into each set. Make sure to set nSets to the correct value before passing `explicitSets`.

A snippet of `params.yaml` is shown below where `distributed_params` parameters are modified to subset the data.
```
## This file holds all parameters to be passed into PyCoGAPS.
...

distributed_params:
  # number of sets to break data into
  nSets: 2
  # specify subsets by index or name
  explicitSets: ['IM00', 'IM02']
```

# **4. Citing CoGAPS**
If you use the CoGAPS package for your analysis, please cite Fertig et al. (2010)

If you use the gene set statistic, please cite Ochs et al. (2009)

# **References**
Fertig, Elana J., Jie Ding, Alexander V. Favorov, Giovanni Parmigiani, and Michael F. Ochs. 2010. “CoGAPS: An R/C++ Package to Identify Patterns and Biological Process Activity in Transcriptomic Data.” Bioinformatics 26 (21): 2792–3. https://doi.org/10.1093/bioinformatics/btq503.

Ochs, Michael F., Lori Rink, Chi Tarn, Sarah Mburu, Takahiro Taguchi, Burton Eisenberg, and Andrew K. Godwin. 2009. “Detection of Treatment-Induced Changes in Signaling Pathways in Gastrointestinal Stromal Tumors Using Transcriptomic Data.” Cancer Research 69 (23): 9125–32. https://doi.org/10.1158/0008-5472.CAN-09-1709.

Seung, Sebastian, and Daniel D. Lee. 1999. “Learning the Parts of Objects by Non-Negative Matrix Factorization.” Nature 401 (6755): 788–91. https://doi.org/10.1038/44565.

Stein-O’Brien, Genevieve L., Jacob L. Carey, Wai S. Lee, Michael Considine, Alexander V. Favorov, Emily Flam, Theresa Guo, et al. 2017. “PatternMarkers & Gwcogaps for Novel Data-Driven Biomarkers via Whole Transcriptome Nmf.” Bioinformatics 33 (12): 1892–4. https://doi.org/10.1093/bioinformatics/btx058.
