#!/usr/bin/python
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  MKL_DIR = os.getenv('MKLROOT')
  IMPI_DIR = os.getenv('I_MPI_ROOT')
  HDF5_DIR = os.getenv('H5ROOT')
  configure_options = [
    # Limit number of cpus used for compilation.
    # This is important for taito-shell.csc.fi (sinteractive) sessions.
    # For sbatch/salloc sessions with full node usage, this can be up to 24 or removed (will default to 13).
    '--with-make-np=4',

    '--with-debugging=1',
    # see https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
    '--with-blaslapack-lib=-L%s/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl' % MKL_DIR,
    '--with-blaslapack-include=%s/include' % MKL_DIR,
    '--with-mpidir=' + IMPI_DIR,
    '--with-mpiexec=srun',
    '--with-scalapack-lib=-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64', # we don't repeat what is already in --with-blaslapack-lib
    '--with-scalapack-include=%s/include' % MKL_DIR,
    '--download-hdf5',
    '--download-metis',
    '--download-mumps',
    '--download-parmetis',
    '--download-superlu',
    '--download-superlu_dist',
    '--with-cxx-dialect=C++11',
  ]
  configure.petsc_configure(configure_options)
