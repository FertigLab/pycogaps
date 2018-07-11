from distutils.core import setup, Extension

Cogaps = Extension('Cogaps', ['Cogapsmodule.c'])

setup(name='Cogaps', version='1.0', ext_modules=[Cogaps])
