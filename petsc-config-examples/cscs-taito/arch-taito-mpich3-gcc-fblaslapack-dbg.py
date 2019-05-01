#!/usr/bin/python
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    # Limit number of cpus used for compilation.
    # This is important for taito-shell.csc.fi (sinteractive) sessions.
    # For sbatch/salloc sessions with full node usage, this can be up to 24 or removed (will default to 13).
    '--with-make-np=4',

    '--with-debugging=1',
    '--download-fblaslapack',
    '--download-hdf5',
    '--download-metis',
    '--download-mpich',
    # This option allows to use srun with mpich.
    # See https://wiki.mpich.org/mpich/index.php/Frequently_Asked_Questions#Q:_How_do_I_use_MPICH_with_slurm.3F
    #     https://slurm.schedmd.com/mpi_guide.html#mpich2
    # CPPFLAGS = workaround for missing pmi.h during mpich compilation
    '--download-mpich-configure-arguments=--with-pmi=slurm --with-pm=no CPPFLAGS=-I/usr/include/slurm', 
    '--with-mpiexec=srun',
    '--download-mumps',
    '--download-parmetis',
    '--download-scalapack',
    '--download-superlu',
    '--download-superlu_dist',
    '--with-cxx-dialect=C++11',
  ]
  configure.petsc_configure(configure_options)
