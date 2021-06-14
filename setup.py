from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import sysconfig
import setuptools

from pybind11.setup_helpers import Pybind11Extension, build_ext
__version__ = '0.0.1'


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
    # Pybind11Extension("cogaps",
    #     ["src/bindings.cpp"],
    #     # Example: passing in the version to the compiled code
    #     define_macros = [('VERSION_INFO', __version__)],
    #     ),
    Pybind11Extension(
        name='cogaps',
        sources=[
            'src/bindings.cpp',
            'src/Rpackage/src/GapsParameters.cpp',
            'src/Rpackage/src/GapsResult.cpp',
            'src/Rpackage/src/GapsRunner.cpp',
            'src/Rpackage/src/GapsStatistics.cpp',
            'src/Rpackage/src/atomic/AtomicDomain.cpp',
            'src/Rpackage/src/atomic/ProposalQueue.cpp',
            'src/Rpackage/src/data_structures/HashSets.cpp',
            'src/Rpackage/src/data_structures/HybridMatrix.cpp',
            'src/Rpackage/src/data_structures/HybridVector.cpp',
            'src/Rpackage/src/data_structures/Matrix.cpp',
            'src/Rpackage/src/data_structures/SparseIterator.cpp',
            'src/Rpackage/src/data_structures/SparseMatrix.cpp',
            'src/Rpackage/src/data_structures/SparseVector.cpp',
            'src/Rpackage/src/data_structures/Vector.cpp',
            'src/Rpackage/src/file_parser/CsvParser.cpp',
            'src/Rpackage/src/file_parser/FileParser.cpp',
            'src/Rpackage/src/file_parser/GctParser.cpp',
            'src/Rpackage/src/file_parser/MtxParser.cpp',
            'src/Rpackage/src/file_parser/TsvParser.cpp',
            'src/Rpackage/src/gibbs_sampler/DenseGibbsSampler.cpp',
            'src/Rpackage/src/gibbs_sampler/SparseGibbsSampler.cpp',
            'src/Rpackage/src/math/Math.cpp',
            'src/Rpackage/src/math/MatrixMath.cpp',
            'src/Rpackage/src/math/Random.cpp',
            'src/Rpackage/src/math/VectorMath.cpp'            
        ],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            './src/Rpackage/src/include/'
        ],
        language='c++'
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

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.9']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        opts.append("-mmacosx-version-min=10.9")
        opts.append("-I src/Rpackage/src/include/*")
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
    name='cogaps',
    version=__version__,
    author='Thomas Sherman',
    author_email='tomsherman159@gmail.com',
    url='https://github.com/FertigLab/pycogaps',
    description='Non-Negative Matrix Factorization Algorithm',
    long_description='',
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.2'],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
    language="c++",
   # extra_compile_args=["-I src/Rpackage/src/include", "-nostdlib" "-undefined dynamic_lookup"],
)