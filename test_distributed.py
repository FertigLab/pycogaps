from PyCoGAPS import *

path = "./data/GIST.csv"
params = CoParams(path)
adata = toAnndata(path)

singlethreadres = CoGAPS(path, params)

if __name__ == '__main__':
    resultplaceholder = distributedCoGAPS(path, params, None)
    distresult = resultplaceholder.get()
    print(distresult)
    # single- and multi-threaded cogaps should yield slightly different results
    # because of different seeds, but we can compare some specific attributes
    assert(singlethreadres["GapsResult"].chisqHistory == distresult["GapsResult"].chisqHistory)
    assert (singlethreadres["anndata"].shape == distresult["anndata"].shape)


# test pickling / unpickling python objects so they can be handled by spawned processes
# leaving this commented out because it's a pain to test, but it's here if we need it
# import pickle
# # params.gaps.print()
# print("Now attempting to pickle GapsParameters object.....")
# pickle.dump(params.gaps, open("./data/save.p", "wb"))
# print("Done. Now attempting to unpickle...")
# unpickled = pickle.load(open("./data/save.p", "rb"))
# unpickled.print()
