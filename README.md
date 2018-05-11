Hands-on examples for PETSc Basic &amp; Advanced Tutorial (PRACE Training Course), May 2018, Ostrava
==========

Welcome to the tutorial!

You can find slides at http://tinyurl.com/petsc-tut-2018.

Notes
-----

DMPlex + FEM example:
http://www.mcs.anl.gov/petsc/petsc-dev/src/snes/examples/tutorials/ex62.c.html

Convert a distributed vector to a sequential one using `VecScatter`:
```
Vec vin; /* distributed vector */
Vec vout; /* sequential vector (to be created) */
VecScatter ctx;
...
ierr = VecScatterCreateToZero(vin,&ctx,&vout);CHKERRQ(ierr);
ierr = VecScatterBegin(ctx,vin,vout,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
ierr = VecScatterEnd(ctx,vin,vout,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
...
ierr = VecDestroy(&vout);CHKERRQ(ierr);
```

If `git pull` fails because of local changes, do:
```
git stash
git pull
git stash pop
```
Try to compile the stashed example again. Sometimes you will have to resolve merge conflicts.

Filter the error output - only from rank 0:
```
mpirun -n 3 ./ex3 2>&1 | grep '^\[0\]'

```

Add this to the python configure script (`arch-impi-imkl-linux-opt.py`) to compile with optimizations:
```
    '--with-debugging=0',
    '--COPTFLAGS=-ipo -O3 -xCORE-AVX2',
    '--CXXOPTFLAGS=-ipo -O3 -xCORE-AVX2',
    '--FOPTFLAGS=-ipo -O3 -xCORE-AVX2',
```
