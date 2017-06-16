from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext as build_pyx

setup(name = 'potts', ext_modules=[Extension('potts', ['potts.pyx'])], cmdclass = { 'build_ext': build_pyx })
