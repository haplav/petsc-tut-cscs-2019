1. `cd ~/petsc-tut-cscs-2019/exercises`
2. `git pull`
3. `make ex1`
4. `srun -n 4 ./ex1`
5. notice direct calls to MPI are possible
6. compare the error output for PetscPrintfwith and without `CHKERRQ`
7. fix the error (change the communicator to `PETSC_COMM_WORLD`)
8. add similar call to `PetscPrintf`, but now with `PETSC_COMM_SELF` communicator – what‘s the problem?
9. add similar call to `PetscSynchronizedPrintf`/`PetscSynchronizedFlush` with `PETSC_COMM_WORLD`
10. set the value of variable myintfrom command line option `–myint` and print its value
11. try setting the same from the petscrcfileand/or environment variable `PETSC_OPTIONS` and overriding the value from command line
12. try out `-help` option
