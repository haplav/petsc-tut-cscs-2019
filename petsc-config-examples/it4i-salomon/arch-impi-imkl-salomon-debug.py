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
    '--CXX=mpicxx',
    '--F77=mpif77',
    '--F90=mpif90',
    '--FC=mpifc',
    '--with-blaslapack-dir=' + os.getenv('MKLROOT'),
    '--with-make-np=8',
  ]
  configure.petsc_configure(configure_options)
