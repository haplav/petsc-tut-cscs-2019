Hands-on examples for the PRACE Training Course "Scalable Parallel Computations with PETSc" held at CSC - IT Center for Science, Espoo, Finland on 2-3 May 2019
=======================================

Welcome to the tutorial!

You can find tutorial information and slides at https://events.prace-ri.eu/event/871/.

Preliminaries
-------------
Login to Taito:
```
ssh trainingXXX@taito.csc.fi
```
Load newer GIT (since the default one makes problems):
```
module add git
```
(Feel free to add this line to your ~/.bashrc file.)

Some GIT basics
---------------
To make a local copy of this repository:
```
# make sure we are at $HOME
cd

# clone into new petsc-tut-cscs-2019 directory
git clone https://github.com/haplav/petsc-tut-cscs-2019.git
```

To update your local copy:
```
cd ~/petsc-tut-cscs-2019       # path to local clone (contains .git subdirectory)
git pull
```

If `git pull` fails because of local changes, stash them first:
```
git stash
git pull
git stash pop
```
Sometimes you will have to resolve merge conflicts (let me know if you struggle with that).
Then try to compile the example again.

Tips & tricks
-------------

### DMPlex + FEM example
  http://www.mcs.anl.gov/petsc/petsc-dev/src/snes/examples/tutorials/ex62.c.html

### Convert a distributed vector to a sequential one using `VecScatter`
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

### Filter the error output - only from rank 0
```
mpirun -n 3 ./ex3 2>&1 | grep '^\[0\]'

```
