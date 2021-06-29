# Installation
```
git clone https://github.com/FertigLab/pycogaps --recursive
cd pycogaps
pip3 install .
```
# Usage

```
from PyCoGAPS import *

# replace with the path to your data, or use this provided example
path = "./data/GIST.csv" 

# run CoGAPS on your dataset
result = CoGAPS(path)

# create a GapsParameters object and set desired parameters
params = GapsParameters(path)
setParam(params, "maxThreads", 4) 
# and/or:
setParams(params, {
            'printMessages': True,
            'maxThreads': 4
        })

result = CoGAPS(path, params)
```
