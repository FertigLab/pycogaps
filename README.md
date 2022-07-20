<img width="285" alt="image" src="https://user-images.githubusercontent.com/25310425/177400924-48b0c78a-16b5-4565-9de7-2a0f2b3d7ac6.png">

# **PyCoGAPS**

Coordinated Gene Activity in Pattern Sets (CoGAPS) implements a Bayesian MCMC matrix factorization algorithm, GAPS, and links it to gene set statistic methods to infer biological process activity. It can be used to perform sparse matrix factorization on any data, and when this data represents biomolecules, to do gene set analysis.

This package, PyCoGAPS, presents a unified Python interface, with a parallel, efficient underlying implementation in C++. The R implementation of CoGAPS can be found here: https://github.com/FertigLab/CoGAPS/

## **Table of Contents**

1. [ Using the PyCoGAPS Library ](#1-using-the-pycogaps-library)
2. [ Running PyCoGAPS Using Docker](#2-running-pycogaps-using-docker)   
3. [ Analyzing the PyCoGAPS Result ](#3-analyzing-the-pycogaps-result)
4. [ Additional Features of PyCoGAPS ](#4-additional-features-of-pycogaps)     
5. [ Citing PyCoGAPS ](#5-citing-pycogaps)



# **1. Using the PyCoGAPS library**

To install, please clone our GitHub repository as follows: 
```
git clone https://github.com/FertigLab/pycogaps.git --recursive
cd pycogaps 
python setup.py install
```
When PyCoGAPS has installed and built correctly, you should see this message:
```
Finished processing dependencies for pycogaps==0.0.1
```
Which means it is ready to use! You may need to install some Python dependencies before everything can build, so don’t be deterred if it takes a couple of tries to install.

We'll first begin with setting up your working environment and running CoGAPS on an example dataset with default parameters.

To use the PyCoGAPS python package, import dependencies as follows:

```
from parameters import *
from pycogaps_main import CoGAPS
import scanpy as sc
```
NOTE: if you wish to run distributed (parallel), please wrap all subsequent code in this check to avoid thread reentry issues:
```
if __name__ == "__main__":
```
Load input data (acceptable formats: h5ad, h5, csv, txt, mtx, tsv)
```
path = "/Users/jeanette/fertiglab/PDAC_Atlas_Pipeline/PDAC.h5ad"
adata = sc.read_h5ad(path)
```
We recommend log normalizing count data
```
sc.pp.log1p(adata)
```
Now, set run parameters by creating a CoParams object. 
```
params = CoParams(path)

setParams(params, {
    'nIterations': 50000,
    'seed': 42,
    'nPatterns': 8,
    'useSparseOptimization': True,
    'distributed': "genome-wide",
})
```
If you are running in parallel, distributed parameters can be modified like this:
```
params.setDistributedParams(nSets=15, minNS=8, maxNS=23, cut=8)
```
Now, start your CoGAPS run by passing your data object and parameter object. Since CoGAPS runs can take significant time to complete, we recommend keeping track of how run times scale with increasing patterns and iterations.
```
start = time.time()
result = CoGAPS(adata, params)
end = time.time()

print("TIME:", end - start)
```
While CoGAPS is running, you will see periodic status messages saying how many iterations have been completed, the current ChiSq value, and how much time has elapsed out of the estimated total runtime.
```
1000 of 50000, Atoms: 5424(A), 21232(P), ChiSq: 138364000, Time: 00:03:47 / 11:13:32
2000 of 50000, Atoms: 5394(A), 20568(P), ChiSq: 133824536, Time: 00:03:46 / 11:10:34
3000 of 50000, Atoms: 5393(A), 21161(P), ChiSq: 133621048, Time: 00:03:51 / 11:25:24
4000 of 50000, Atoms: 5527(A), 22198(P), ChiSq: 137671296, Time: 00:04:00 / 11:52:06
5000 of 50000, Atoms: 5900(A), 20628(P), ChiSq: 137228688, Time: 00:03:58 / 11:46:10
```
When the run is finished, CoGAPS will print a message like this:
```
GapsResult result object with 5900 features and 20628 samples
8 patterns were learned
```
We strongly recommend saving your result object as soon as it returns. One option to do so is using Python’s serialization library, pickle64.
```
print("Pickling...")
pickle.dump(result, open("./data/PDACresult_50kiterations.pkl", "wb"))
print("Pickling complete!")
```
Now you have successfully generated a CoGAPS result! To continue to visualization and analysis guides, please skip to the section below titled “Analyzing the PyCoGAPS Result” 

# **2. Running PyCoGAPS using Docker**

The second option for running PyCoGAPS is using a Docker image, which we will pull from the Docker repository, and this contains a set of instructions to build and run PyCoGAPS. With this Docker image, there's no need to install any dependencies, import packages, etc. as the environment is already set up and directly ready to run on your computer.

Please follow the steps below to run the PyCoGAPS vignette on Mac/Linux OS:
1. Install Docker at https://docs.docker.com/desktop/mac/install/ 
2. Open the Docker application or paste the following in terminal:
```
docker run -d -p 80:80 docker/getting-started
```
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
For MARCC users, we'll be building the pycogaps package and installing all dependencies in a conda environment.

Please follow the steps below to run the PyCoGAPS vignette on MARCC:
1. Copy the commands and paste in terminal   
**Note:** If you're encountering an error with `import pycogaps`, make sure that you have the intel/18.0 core module loaded instead of gcc.

```
git clone --recurse-submodules https://github.com/FertigLab/pycogaps.git
ml anaconda
conda create --name pycogaps python=3.8
conda activate pycogaps
cd pycogaps
pip install -r requirements.txt --user
python3 setup.py install
python3 vignette.py

```

This produces a CoGAPS run on a simple dataset with default parameters. You should then see the following output:
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

Pickling complete!
```

CoGAPS has successfully completed running and has saved the result file as `result.pkl` in a created `output/` folder.

Your working directory is the `PyCoGAPS` folder with the following structure and files:
```
PyCoGAPS
├── data
│   └── GIST.csv
├── params.yaml
├── output
│   └── result.pkl
```

Now, you're ready to run CoGAPS for analysis on your own data with custom parameters. 

In order to analyze your desired data, we'll need to input it and modify the default parameters before running CoGAPS. All parameter values can be modified directly in the `params.yaml` file already downloaded earlier. 

Please follow the steps below to run PyCoGAPS with custom parameters:
1. Open `params.yaml` with any text or code editor
2. Modify the `path` parameter value by replacing the default `data/GIST.csv` with `data/your-datafile-name`  
  **Note**: Make sure you have moved your data into the created `data/` folder
3. Modify any additional desired parameters and save
4. For Mac/Linux OS, run the following in terminal:
```
docker run -v $PWD:$PWD ashleyt2000/pycogaps:docker_pycogaps $PWD/params.yaml
```
4. For MARCC, run the following in terminal:
```
python3 vignette.py
```
## **Example Snippet of `params.yaml`**

A snippet of `params.yaml` is shown below, where we have changed some default parameter values to our own specified example values.

```
## This file holds all parameters to be passed into PyCoGAPS.
## To modify default parameters, simply replace parameter values below with user-specified values, and save file. 

# RELATIVE path to data -- make sure to move your data into the created data/ folder
path: data/liver_dataset.txt

# result output file name 
result_file: liver_result.pkl

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

## **Running CoGAPS in Parallel**
Non-Negative Matrix Factorization algorithms typically require long computation times and CoGAPS is no exception. In order to scale CoGAPS up to the size of data sets seen in practice we need to take advantage of modern hardware and parallelize the algorithm.

### **I. Multi-Threaded Parallelization**
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

### **II. Distributed CoGAPS**
For large datasets (greater than a few thousand genes or samples) the multi-threaded parallelization isn’t enough. It is more efficient to break up the data into subsets and perform CoGAPS on each subset in parallel, stitching the results back together at the end (Stein-O’Brien et al. (2017)).

In order to use these extensions, some additional parameters are required, specifically modifying the `distributed_params` in `params.yaml`. We first need to set `distributed`  to be `genome-wide.` Next, `nSets` specifies the number of subsets to break the data set into. `cut`, `minNS`, and `maxNS` control the process of matching patterns across subsets and in general should not be changed from defaults. More information about these parameters can be found in the original papers. 

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


# **3. Analyzing the PyCoGAPS Result**

## **Breaking Down the Result Object from CoGAPS**

<img src="https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/anndata-result.png" alt="anndata result obj" width="300" align="center">

A dictionary of the result as two representations is stored: an `anndata` object. CoGAPS stores the lower dimensional representation of the samples (P matrix) in the `.var` slot and the weight of the features (A matrix) in the `.obs` slot. The standard deviation across sample points for each matrix are stored in the `.uns` slots.

## **Analyzing the Result**
We provide two ways to analyze the result from PyCoGAPS. The first includes an interactive notebook interface using the web-based GenePattern Notebook (recommended for less experienced python/programming users), and the secoond includes running a python script from the command line (recommended for more experienced python/programming users). 

## **I. GenePattern Notebook** 
Here are the following steps to use the interactive GenePattern Notebook to analyze results:
1. Go to the PyCoGAPS Analysis Notebook found here: https://notebook.genepattern.org/hub/preview?id=440. 
3. Click 'Run' and open 'PyCoGAPS Analysis.ipynb'
4. Follow the instructions in the notebook to run your analysis.

## **II. Python Script via Terminal** 
In order to analyze the data, we'll need to make sure to install the necessary dependencies and import the built-in PyCoGAPS functions. 

Make sure you're in the `PyCoGAPS` folder, and copy the following commands in terminal, which will save plots generated from the example data in the `output/` folder:
```
cd output
curl -O https://raw.githubusercontent.com/FertigLab/pycogaps/master/PyCoGAPS/analysis_functions.py
curl -O https://raw.githubusercontent.com/FertigLab/pycogaps/master/PyCoGAPS/requirements_analysis.txt
pip install -r requirements_analysis.txt --user
python3 analysis_functions.py result.pkl
```
To analyze a different result, replace `result.pkl` with the path to your desired result file in the command line.

## **More on the Analysis Functions**
Below details each of the analysis functions included in the package. 

3.2.1 [ Default Plot ](#241-default-plot)  
3.2.2 [ Residuals Plot ](#242-residuals-plot)  
3.2.3 [ Pattern Markers Plot ](#243-pattern-markers-plot)  
3.2.4 [ Binary Plot ](#244-binarya-plot)  
3.2.5 [ Calculate CoGAPS Statistics ](#245-calculate-cogaps-statistics)  
3.2.6 [ Calculate Gene GSS Statistic ](#246-calculate-gene-gss-statistic)  
3.2.7 [ Calculate Gene GS Probability ](#247-calculate-gene-gs-probability)  



### **I. Default Plot**
By default, the `plot` function displays how the patterns vary across the samples.

```python
# plot result object returned from CoGAPS
plot(result)
```
<img src="https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/res_show.png" alt="show result function" width="400" align="center">





### **II. Residuals Plot**
`plotResiduals` calculates residuals and produces a heatmap.

```python
plotResiduals(result)
```

<img src="https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/plot_residuals.png" alt="plot residuals" width="400" align="center">


### **III. Pattern Markers Plot**
`plotPatternMarkers` plots a heatmap of the original data clustered by the pattern markers statistic, which computes the most associated pattern for each gene.

```python
plotPatternMarkers(result, legend_pos=None)
```

<img src="https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/plot_pm.png" alt="plot pattern markers" width="400" align="center">


### **IV. Binary Plot**
`binaryA` creates a binarized heatmap of the A matrix in which the value is 1 if the value in Amean is greater
than `threshold * Asd` and 0 otherwise.

```python
binaryA(result, threshold=3)
```

<img src="https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/binaryA.png" alt="plot binary hm" width="400" align="center">


```python
# plotting clustered binary plot
binaryA(result, threshold=3, cluster=True)
```

<img src="https://github.com/FertigLab/pycogaps/blob/update-setup-instructions/rm/binaryA_cluster.png" alt="plot binary hm, cluster" width="400" align="center">





### **V. Calculate CoGAPS Statistics**
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

### **VI. Calculate Gene GSS Statistic**
`calcGeneGSStat` calculates the probability that a gene listed in a gene set behaves like other genes in the set within the given data set.

```python
stats = calcGeneGSStat(result, GStoGenes=['Hs.101174', 'Hs.1012'], numPerm=1000)
```

```
Hs.101174  0.422955
Hs.1012    0.391747
```

### **VII. Compute Gene GS Probability**

`computeGeneGSProb` computes the p-value for gene set membership using the CoGAPS-based statistics developed in Fertig et al. (2012). This statistic refines set membership for each candidate gene in a set specified in GSGenes by comparing the inferred activity of that gene to the average activity of the set.

```python
stats = computeGeneGSProb(result, GStoGenes=['Hs.101174', 'Hs.1012'])
```

```
Hs.101174  0.617193
Hs.1012    0.887583
```

# **4. Additional Features of PyCoGAPS**

## **Checkpoint System: Saving/Loading CoGAPS Runs**
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

## **Transposing Data**
If your data is stored as samples x genes, CoGAPS allows you to pass `transposeData: True` and will automatically read the transpose of your data to get the required genes x samples configuration. 

A snippet of `params.yaml` is shown where we now load the saved CoGAPS checkpoint file to continue the run.
```
## This file holds all parameters to be passed into PyCoGAPS.
...

run_params:
  # T/F for transposing data while reading it in - useful for data that is stored as samples x genes since CoGAPS requires data to be genes x samples
  transposeData: True
```

## **Passing Uncertainty Matrix**
In addition to providing the data, the user can also specify an uncertainty measurement - the standard deviation of each entry in the data matrix. By default, CoGAPS assumes that the standard deviation matrix is 10% of the data matrix. This is a reasonable heuristic to use, but for specific types of data you may be able to provide better information. Make sure to save your uncertainty file into the `data/` file.

A snippet of `params.yaml` is shown where we now load the saved CoGAPS checkpoint file to continue the run.
```
## This file holds all parameters to be passed into PyCoGAPS.
...

run_params:
  # uncertainty matrix - either a matrix or a supported file type
  uncertainty: data/GIST_uncertainty.csv
```

## **Subsetting Data**
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

# **5. Citing PyCoGAPS**
If you use the CoGAPS package for your analysis, please cite Fertig et al. (2010)

If you use the gene set statistic, please cite Ochs et al. (2009)

# **References**
Fertig, Elana J., Jie Ding, Alexander V. Favorov, Giovanni Parmigiani, and Michael F. Ochs. 2010. “CoGAPS: An R/C++ Package to Identify Patterns and Biological Process Activity in Transcriptomic Data.” Bioinformatics 26 (21): 2792–3. https://doi.org/10.1093/bioinformatics/btq503.

Ochs, Michael F., Lori Rink, Chi Tarn, Sarah Mburu, Takahiro Taguchi, Burton Eisenberg, and Andrew K. Godwin. 2009. “Detection of Treatment-Induced Changes in Signaling Pathways in Gastrointestinal Stromal Tumors Using Transcriptomic Data.” Cancer Research 69 (23): 9125–32. https://doi.org/10.1158/0008-5472.CAN-09-1709.

Seung, Sebastian, and Daniel D. Lee. 1999. “Learning the Parts of Objects by Non-Negative Matrix Factorization.” Nature 401 (6755): 788–91. https://doi.org/10.1038/44565.

Stein-O’Brien, Genevieve L., Jacob L. Carey, Wai S. Lee, Michael Considine, Alexander V. Favorov, Emily Flam, Theresa Guo, et al. 2017. “PatternMarkers & Gwcogaps for Novel Data-Driven Biomarkers via Whole Transcriptome Nmf.” Bioinformatics 33 (12): 1892–4. https://doi.org/10.1093/bioinformatics/btx058.
