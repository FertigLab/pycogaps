from PyCoGAPS import *

path = "./data/GIST.csv"
params = CoParams(path)
adata = toAnndata(path)

if __name__ == '__main__':
    resultplaceholder = distributedCoGAPS(path, params, None)
    result = resultplaceholder.get()
    print(result)