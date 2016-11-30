/* Adopted from http://www.mcs.anl.gov/petsc/petsc-3.7/src/vec/vec/examples/tutorials/ex2.c.html */

static char help[] = "Builds a parallel vector with 1 component on the first processor, 2 on the second, etc.\n\
  Then each processor adds one to all elements except the last rank.\n\n";

/*
  Include "petscvec.h" so that we can use vectors.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscis.h     - index sets
     petscviewer.h - viewers
*/
#include <petscvec.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  PetscInt       i,N,lo,hi;
  PetscScalar    integer = 1.0;
  Vec            x;

  PetscInitialize(&argc,&argv,(char*)0,help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  /*
     Create a parallel vector.
      - In this case, we specify the size of each processor's local
        portion, and PETSc computes the global size.  Alternatively,
        if we pass the global size and use PETSC_DECIDE for the
        local size PETSc will choose a reasonable partition trying
        to put nearly an equal number of elements on each processor.
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  //TODO task 3
  ierr = VecSetSizes(x,rank+1,rank+1);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecSet(x,1.0);CHKERRQ(ierr);

  /* get the computed global size */
  ierr = VecGetSize(x,&N);CHKERRQ(ierr);

  /*
     Set the vector elements.
      - Always specify global locations of vector entries.
      - Each processor can contribute any vector entries,
        regardless of which processor "owns" them; any nonlocal
        contributions will be transferred to the appropriate processor
        during the assembly process.
      - In this example, the flag ADD_VALUES indicates that all
        contributions will be added together.
  */
  //TODO task 8
  for (i=0; i<N; i++) {
    ierr = VecSetValues(x,1,&i,&integer,ADD_VALUES);CHKERRQ(ierr);
  }

  /*
     Assemble vector, using the 2-step process:
       VecAssemblyBegin(), VecAssemblyEnd()
     Computations can be done while messages are in transition
     by placing code between these two statements.
  */
  //TODO task 5

  /*
      View the vector.
  */
  ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  //TODO tasks 9 and 10
  ierr = VecZeroEntries(x);CHKERRQ(ierr);
  //ierr = VecGet...(x,&lo,&hi);CHKERRQ(ierr);
  for (i=0; i<N; i++) {
    ierr = VecSetValues(x,1,&i,&integer,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
  ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  //TODO task 11
  integer = rank+1;
  //for () {
  //
  //}
  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
  //ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /*
      Destroy the vector and finalize PETSc.
  */
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

