import cppinterface

# python interface to CoGAPS
# this is what the outside world will call
# input: path to data file
# output: right now it just prints "hooty hoo"
def runCoGAPS(filepath):
    print("hooty hoo")
    print(filepath)
    print(cppinterface.cogaps_from_file_cpp(filepath))
    return 1