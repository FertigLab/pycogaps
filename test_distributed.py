from PyCoGAPS import *

path = "./data/GIST.csv"
params = CoParams(path)
adata = toAnndata(path)

singlethreadres = CoGAPS(path, params)

if __name__ == '__main__':
    params.setDistributedParams(nSets=2, minNS=1, cut=3)
    params.coparams['subsetIndices'] = subset_data.createSets(adata, params)
    result = distributedCoGAPS(path, params, None)
    print(result)
    print("length: ", len(result))
    # print("Parallel chisqhistory:", result[0]["GapsResult"].chisqHistory, "\n")
    # print("Single-thread chisqhistory:", singlethreadres["GapsResult"].chisqHistory, "\n")
    # assert(singlethreadres["GapsResult"].chisqHistory == result[0]["GapsResult"].chisqHistory)
    # assert (singlethreadres["anndata"].shape == result[0]["anndata"].shape)



# test pickling / unpickling python objects so they can be handled by spawned processes
# leaving this commented out because it's a pain to test, but it's here if we need it
# import pickle
# # params.gaps.print()
# print("Now attempting to pickle GapsParameters object.....")
# pickle.dump(params.gaps, open("./data/save.p", "wb"))
# print("Done. Now attempting to unpickle...")
# unpickled = pickle.load(open("./data/save.p", "rb"))
# unpickled.print()
