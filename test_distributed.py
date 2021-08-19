from PyCoGAPS import *

path = "./data/GIST.csv"
params = CoParams(path)
adata = toAnndata(path)

# test pickling / unpickling python objects so they can be handled by spawned processes
# import pickle
# # params.gaps.print()
# print("Now attempting to pickle GapsParameters object.....")
# pickle.dump(params.gaps, open("./data/save.p", "wb"))
# print("Done. Now attempting to unpickle...")
# unpickled = pickle.load(open("./data/save.p", "rb"))
# unpickled.print()

if __name__ == '__main__':
    resultplaceholder = distributedCoGAPS(path, params, None)
    result = resultplaceholder.get()
    print(result)