from PyCoGAPS import *

def no_na_in_result(result):
    return (np.isnan(toNumpy(result['GapsResult'].Amean)).any()+
    np.isnan(toNumpy(result['GapsResult'].Asd)).any()+
    np.isnan(toNumpy(result['GapsResult'].Pmean)).any()+
    np.isnan(toNumpy(result['GapsResult'].Psd)).any())

csv_path = "./data/GIST.csv" 
mtx_path = "./data/GIST.mtx" 
tsv_path = "./data/GIST.tsv" 


csv_params = CoParams(csv_path)
mtx_params = CoParams(mtx_path)
tsv_params = CoParams(tsv_path)

setParams(csv_params, {"nIterations": 100,
                    "seed": 1,
                    "nPatterns": 7,
                    'hdfKey': 'counts'})
setParams(mtx_params, {"nIterations": 100,
                    "seed": 1})
setParams(tsv_params, {"nIterations": 100,
                    "seed": 1})

csv_params.print_all()

res = [None] * 3
res[0] = CoGAPS(csv_path, csv_params, outputFrequency=50, messages=False)
res[1] = CoGAPS(mtx_path, mtx_params, outputFrequency=50, messages=False)
res[2] = CoGAPS(tsv_path, tsv_params, outputFrequency=50, messages=False)

for r in res:
    assert(no_na_in_result(r) == False)

assert(toNumpy(res[0]['GapsResult'].Amean).shape[0] == 1363)
assert(toNumpy(res[0]['GapsResult'].Amean).shape[1] == 7)
assert(toNumpy(res[0]['GapsResult'].Pmean).shape[0] == 9)
assert(toNumpy(res[0]['GapsResult'].Pmean).shape[1] == 7)

