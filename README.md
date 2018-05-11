Hands-on examples for PETSc Basic &amp; Advanced Tutorial (PRACE Training Course), May 2018, Ostrava
==========

Welcome to the tutorial!

You can find slides at http://tinyurl.com/petsc-tut-2018.

Notes
-----

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
