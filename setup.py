# Cython compile instructions

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

# Use python setup.py build --inplace
# to compile

setup(ext_modules = cythonize(Extension(
           "denoising",                                # the extension name
           sources=["denoising.pyx", "cpp/ising.cpp", "cpp/potts.cpp"], # the Cython source and
                                                  # additional C++ source files
           language="c++",                        # generate and compile C++ code
           extra_compile_args=["-std=c++11"],
           extra_link_args=["-std=c++11"]
      )))
