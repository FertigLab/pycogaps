if __name__ == '__main__':
    from PyCoGAPS import *
    import pickle
    print("This vignette was built using pycogaps version", getVersion())

    path = "data/GIST.csv"
    params = CoParams(path)

    setParams(params, {
        'nIterations': 1000,
        'seed': 0,
        'nPatterns': 3,
        # 'useSparseOptimization': True,
        'distributed': "genome-wide",
    })

    start = time.time()
    params.setDistributedParams(nSets=4, minNS=2, maxNS=6, cut=3)
    result = CoGAPS(path, params)
    end = time.time()
    print("TIME:", end - start)
            
    print("Pickling...")
    pickle.dump(result, open("./data/simple_vig.pkl", "wb"))
    print("Pickling complete!")
 
