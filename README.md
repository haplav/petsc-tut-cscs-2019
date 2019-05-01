Hands-on examples for the PRACE Training Course "Scalable Parallel Computations with PETSc" held at CSC - IT Center for Science, Espoo, Finland on 2-3 May 2019
=======================================

Welcome to the tutorial!

You can find tutorial information and slides at https://events.prace-ri.eu/event/871/.

During the tutorial, please open the latest version of slides: https://bit.ly/2VCmLih.

Login to Taito
-------------
```
ssh trainingXXX@taito.csc.fi
```

Start interactive session at Taito
----------------------------------
Before running any example or compiling PETSc, start an interactive session on a compute node:
```
salloc -n 4 -t01:00:00
```
and at the end of each hands-on session
```
exit
```

Some GIT basics
---------------
First of all, load newer GIT (since the default one causes troubles):
```
module add git
```
(I suggest adding this line to the end of your `~/.bashrc` file.)

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

### Filter the error output - only from rank 0
```
srun -n 4 ./ex3 2>&1 | grep '^\[0\]'

```

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
