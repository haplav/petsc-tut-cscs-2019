#!/usr/bin/env python

if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--with-matlab=1',
    '--with-matlab-arch=glnxa64',
    '--with-matlab-engine=1',
    '--with-matlab-socket=0',
    '--with-blaslapack-dir=/apps/all/MATLAB/2015b-EDU',
    '--with-matlabengine-lib=-Wl,-rpath,/apps/all/MATLAB/2015b-EDU/sys/os/glnxa64:/apps/all/MATLAB/2015b-EDU/bin/glnxa64:/apps/all/MATLAB/2015b-EDU/extern/lib/glnxa64 -L/apps/all/MATLAB/2015b-EDU/bin/glnxa64 -L/apps/all/MATLAB/2015b-EDU/extern/lib/glnxa64 -leng -lmex -lmx -lmat -lut -lmwm_dispatcher -lmwopcmodel -lmwservices -lmwservices -lmwopcmodel -lmwopcmodel -lmwm_dispatcher -lmwmpath -lmwopcmodel -lmwservices -lmwopcmodel -lmwservices -lxerces-c',
    '--with-mpi=0',
    '--with-shared-libraries=1',
  ]
  configure.petsc_configure(configure_options)

