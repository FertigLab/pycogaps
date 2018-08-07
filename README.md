# Compilation

This makefile is pretty simple so you need to manually specify the path to your python C library. The target is `CogapsPy.so`. An Example compilation command looks like:

`make CogapsPy.so PYTHON_INC=/usr/include/python2.7/`.

To change the compiler set `CXX`.


Requires Python 2.7

To use, import Cogaps class from Cogaps module:

from Cogaps import Cogaps

Instantiate an instance of Cogaps by providing path for data file as the parameter:

result = Cogaps("data/GIST.tsv")

The Cogaps class currently has two graphing methods that can be called as follows:

result.graphPatternMarkerStats(pumpThreshold)

and

result.graphPmean()

The pumpThreshold parameter can either be an integer '0' for PUMP_UNIQUE or an integer '1' for PUMP_CUT.