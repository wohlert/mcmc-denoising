# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
# distutils: sources = cpp/ising.cpp cpp/potts.cpp

from libcpp.vector cimport vector

# c++ interface to cython
cdef extern from "cpp/denoising.hpp" namespace "denoising":
  cdef cppclass Ising:
    Ising(vector[vector[float]], float, float) except +
    vector[vector[int]] metropolisHastings(int)

  cdef cppclass Potts:
    Potts(vector[vector[float]], float, float, int) except +
    vector[vector[float]] metropolisHastings(int)
    vector[vector[float]] metropolisHastings(int, vector[vector[float]])
    vector[vector[float]] MAP(int, float, float)
    vector[vector[vector[float]]] getHistory()

# creating a cython wrapper class
cdef class IsingMH:
  cdef Ising *thisptr      # hold a C++ instance which we're wrapping
  def __cinit__(self, y, beta, sigma):
      self.thisptr = new Ising(y, beta, sigma)
  def __dealloc__(self):
      del self.thisptr
  def metropolisHastings(self, iterations):
      return self.thisptr.metropolisHastings(iterations)

# creating a cython wrapper class
cdef class PottsMH:
  cdef Potts *thisptr
  def __cinit__(self, y, beta, sigma, bins):
      self.thisptr = new Potts(y, beta, sigma, bins)
  def __dealloc__(self):
      del self.thisptr
  def metropolisHastings(self, iterations):
      return self.thisptr.metropolisHastings(iterations)
  def metropolisHastings2(self, iterations, x):
      return self.thisptr.metropolisHastings(iterations, x)
  def MAP(self, iterations, init=4, diffusion=0.995):
      return self.thisptr.MAP(iterations, init, diffusion)
  def history(self):
      return self.thisptr.getHistory()
