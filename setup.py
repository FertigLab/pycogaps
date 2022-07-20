import setuptools
from setuptools import setup
from setuptools.command.build_ext import build_ext
import sys
import pybind11
from pybind11.setup_helpers import Pybind11Extension, build_ext

__version__ = '0.0.1'
import sys
print(sys.path)

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


ext_modules = [
    Pybind11Extension("pycogaps",
                      ['src/bindings.cpp',
                       'src/CoGAPS/src/GapsParameters.cpp',
                       'src/CoGAPS/src/GapsResult.cpp',
                       'src/CoGAPS/src/GapsRunner.cpp',
                       'src/CoGAPS/src/GapsStatistics.cpp',
                       'src/CoGAPS/src/atomic/Atom.cpp',
                       'src/CoGAPS/src/atomic/ConcurrentAtom.cpp',
                       'src/CoGAPS/src/atomic/AtomicDomain.cpp',
                       'src/CoGAPS/src/atomic/ConcurrentAtomicDomain.cpp',
                       'src/CoGAPS/src/atomic/ProposalQueue.cpp',
                       'src/CoGAPS/src/data_structures/HashSets.cpp',
                       'src/CoGAPS/src/data_structures/HybridMatrix.cpp',
                       'src/CoGAPS/src/data_structures/HybridVector.cpp',
                       'src/CoGAPS/src/data_structures/Matrix.cpp',
                       # 'src/CoGAPS/src/data_structures/Matrix.h',
                       'src/CoGAPS/src/data_structures/SparseIterator.cpp',
                       'src/CoGAPS/src/data_structures/SparseMatrix.cpp',
                       'src/CoGAPS/src/data_structures/SparseVector.cpp',
                       'src/CoGAPS/src/data_structures/Vector.cpp',
                       'src/CoGAPS/src/file_parser/CharacterDelimitedParser.cpp',
                       'src/CoGAPS/src/file_parser/FileParser.cpp',
                       'src/CoGAPS/src/file_parser/MtxParser.cpp',
                       'src/CoGAPS/src/file_parser/MatrixElement.cpp',
                       'src/CoGAPS/src/gibbs_sampler/DenseNormalModel.cpp',
                       'src/CoGAPS/src/gibbs_sampler/SparseNormalModel.cpp',
                       'src/CoGAPS/src/gibbs_sampler/AlphaParameters.cpp',
                       'src/CoGAPS/src/math/Math.cpp',
                       'src/CoGAPS/src/math/MatrixMath.cpp',
                       'src/CoGAPS/src/math/Random.cpp',
                       'src/CoGAPS/src/math/VectorMath.cpp',
                       'src/CoGAPS/src/test-runner.cpp'
                       ],
                      include_dirs=[
                          # Path to pybind11 headers
                          get_pybind_include(),
                          get_pybind_include(user=True),
                          'src/CoGAPS/src/include/',
                          'src/CoGAPS/src/',
                          'src/CoGAPS/src/data_structures'
                      ],
                      # depends=[
                      #   'src/CoGAPS/src/data_structures/Matrix.h',
                      #   'src/CoGAPS/src/data_structures/Vector.h'
                      # ],
                      language="c++"
                      ),
]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is required.
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    else:
        raise RuntimeError('Unsupported compiler -- C++14 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.13']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        opts.append("-I/src/CoGAPS/src/include/")
        opts.append("-I/src/CoGAPS/src/*")
        opts.append("-I/src/CoGAPS/src/data_structures/")
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)


setup(
    name='pycogaps',
    version=__version__,
    author='Jeanette Johnson',
    author_email='jjohn450@jhmi.edu',
    url='https://github.com/FertigLab/pycogaps',
    description='Python interface to the Non-Negative Matrix Factorization Algorithm CoGAPS',
    long_description='',
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.2'],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
    language="c++"
)
