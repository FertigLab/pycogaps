from distributed_cogaps import *


def calldistributed(path, params):
    print("name is ", __name__)
    if __name__ == "distcaller":
        return runDistributedCoGAPS(path, params)