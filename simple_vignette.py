if __name__ == '__main__':
    from PyCoGAPS import *
    import pickle
    print("This vignette was built using pycogaps version", getVersion())

    path = "data/GIST.csv"
    params = CoParams(path)

    setParams(params, {
        'nIterations': 5000,
        'seed': 42,
        'nPatterns': 10,
        'useSparseOptimization': True,
        'distributed': "genome-wide",
    })

    start = time.time()
    params.setDistributedParams(nSets=10, minNS=8, maxNS=23)
    result = CoGAPS(path, params)
    end = time.time()
    print("TIME:", end - start)
            
    print("Pickling...")
    pickle.dump(result, open("./data/simple_vig.pkl", "wb"))
    print("Pickling complete!")
 
