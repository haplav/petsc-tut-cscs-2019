1.  `cd ~/petsc-tut-cscs-2019/exercises && git pull`
2.  `make ex3 && mpirun -n 3 ./ex3  # what's wrong?`
3.  Set local size to `rank+1`, let the global size be computed by PETSc (`PETSC_DECIDE`).
4.  `make ex3 && mpirun -n 3 ./ex3  # what's wrong now?`
5.  Fix using `VecAssemblyBegin()`/`VecAssemblyBeginEnd()`.
6.  `make ex3 && mpirun -n 3 ./ex3`
7.  Look at the output of `VecView()`.
8.  All values of vector `x` are equal to `commsize+1`, why?
9.  Get the ownership range into `lo` and `hi` variables.
10. Set each process' values of `x` to the value `rank+1` just by changing the loop lower limit to `lo`.
11. Do the same, avoiding communication (loop from `lo` to `hi`, use `INSERT_VALUES` and `rank`).

### TIP: Filter the error output - only from rank 0 
```
mpirun -n 3 ./ex3 2>&1 | grep '\[0\]'
```
