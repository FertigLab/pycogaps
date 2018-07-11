from distutils.core import setup, Extension

test = Extension('boosttest', ['testmodule.cpp'])

setup(name='boosttest', version='1.0', ext_modules=[test])
