1. `cd ~/petsc-tut-cscs-2019/exercises && git pull`
2. `make ex2`
3. `srun -n 4 ./ex2`
4. Write to `filetest.txt` (google `PetscViewerFileSetName`).
5. Set the file name dynamically from command line option `–f`, using the `filename` variable (already prepared) and `PetscOptionsGetString()` function.
6. Play around with the `–int` option. What’s the default? What's the effect of changing it?
7. `PetscViewerASCIISynchronizedPrintf()` output looks weird. What’s wrong? Fix it! (Use `PetscViewerFlush()`.)
8. Notice `PetscViewerASCIIPushSynchronized()`/`PetscViewerASCIIPopSynchronized()`. You can look at their manpages.
