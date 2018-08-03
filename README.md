# Compilation

This makefile is pretty simple so you need to manually specify the path to your python C library. The target is `CogapsPy.so`. An Example compilation command looks like:

`make CogapsPy.so PYTHON_INC=/usr/include/python2.7/`.

To change the compiler set `CXX`.