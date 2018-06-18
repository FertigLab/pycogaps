from distutils.core import setup, Extension
setup(name='Cogaps', version='1.0',  \
      ext_modules=[Extension('Cogaps', ['Cogapsmodule.cpp'])])
