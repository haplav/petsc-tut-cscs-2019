#!/usr/bin/python

# works on salomon cluster with
# module load git intel/2018a

if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--CC=mpicc',
    '--COPTFLAGS=-ipo -O3 -xCORE-AVX2',
    '--CXX=mpicxx',
    '--CXXOPTFLAGS=-ipo -O3 -xCORE-AVX2',
    '--F77=mpif77',
    '--F90=mpif90',
    '--FC=mpifc',
    '--FOPTFLAGS=-ipo -O3 -xCORE-AVX2',
    '--with-blaslapack-dir=' + os.getenv('MKLROOT'),
    '--with-debugging=0',
    '--with-make-np=8',
  ]
  configure.petsc_configure(configure_options)
