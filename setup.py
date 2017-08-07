from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

basic_ext = Extension(
    'electradpy/sic_powerlaw',
    sources = ['csource/src/radiator.cxx','csource/src/synchrotron_impulsive_continuous.cxx'],
    language='c++',
    include_dirs=[numpy.get_include()],
    extra_compile_args=['-I./csource/incl', "-I/usr/local/include", "-Ofast", "-std=c++11"],
    extra_link_args=["-L/usr/local/lib", '-lgsl', '-lgslcblas']
)

ext_modules = [basic_ext]

setup(ext_modules=cythonize(ext_modules),
      name='elextradpy'

)


