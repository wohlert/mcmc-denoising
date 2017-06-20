# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
# distutils: sources = ising.cpp potts.cpp

from libcpp.vector cimport vector

# c++ interface to cython
cdef extern from "denoising.hpp" namespace "denoising":
  cdef cppclass Ising:
    Ising(vector[vector[float]], float, float) except +
    vector[vector[int]] solve(int)

  cdef cppclass Potts:
    Potts(vector[vector[float]], float, float, int) except +
    vector[vector[float]] solve(int)

# creating a cython wrapper class
cdef class IsingMH:
  cdef Ising *thisptr      # hold a C++ instance which we're wrapping
  def __cinit__(self, y, beta, sigma):
      self.thisptr = new Ising(y, beta, sigma)
  def __dealloc__(self):
      del self.thisptr
  def solve(self, iterations):
      return self.thisptr.solve(iterations)

# creating a cython wrapper class
cdef class PottsMH:
  cdef Potts *thisptr
  def __cinit__(self, y, beta, sigma, bins):
      self.thisptr = new Potts(y, beta, sigma, bins)
  def __dealloc__(self):
      del self.thisptr
  def solve(self, iterations):
      return self.thisptr.solve(iterations)
