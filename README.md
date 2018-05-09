Hands-on examples for PETSc Basic &amp; Advanced Tutorial (PRACE Training Course), May 2018, Ostrava
==========

Welcome to the tutorial!

You can find slides at http://tinyurl.com/petsc-tut-2018.

Notes
-----

Add this to the python configure script (`arch-impi-imkl-linux-opt.py`) to compile with optimizations:
```
    '--with-debugging=0',
    '--COPTFLAGS=-ipo -O3 -xCORE-AVX2',
    '--CXXOPTFLAGS=-ipo -O3 -xCORE-AVX2',
    '--FOPTFLAGS=-ipo -O3 -xCORE-AVX2',
```
