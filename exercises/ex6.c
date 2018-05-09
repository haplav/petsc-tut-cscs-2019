/* This does the same as ex5 but using DMDA. */

static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  DM             da;
  Vec            x, b;        /* approx solution, RHS */
  Vec            d;
  Mat            A;           /* linear system matrix */
  KSP            ksp;         /* linear solver context */
  PC             pc;          /* preconditioner context */
  PetscErrorCode ierr;
  PetscInt       i, n = 5, N;
  PetscInt       its;
  PetscInt       row[2];
  PetscInt       nel, nen;
  const PetscInt *e;
  PetscScalar    rho;
  PetscScalar    value[4], bvalue[2];
  PetscBool      nonzeroguess = PETSC_FALSE;

  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-nonzero_guess",&nonzeroguess,NULL);CHKERRQ(ierr);


  /*
     Create DMDA context for structured 1D problem.
  */
  ierr = MPI_Allreduce(&n, &N, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,N,1,1,NULL,&da);CHKERRQ(ierr);
  ierr = DMSetFromOptions(da);CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);

  /*
     Create vectors.  We get 1 vector from DM and
     then duplicate as needed.
  */
  ierr = DMCreateGlobalVector(da, &x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

  /*
     Set the same value to all vector entries.
  */
  if (nonzeroguess) {
    /* Set nonzero initial guess. Note we use ugly one here to affect number of iterations. */
    ierr = VecSet(x,-1e3);CHKERRQ(ierr);
  } else {
    ierr = VecSet(x,0.0);CHKERRQ(ierr);
  }
  ierr = VecSet(b,0.0);CHKERRQ(ierr);
    
  /*
     Get preallocated matrix from DM. 
  */
  ierr = DMCreateMatrix(da,&A);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);

  /*
    Gets an array containing the indices (in local coordinates) 
    of all the local elements.
    nel	- number of local elements
    nen	- number of one element's nodes
    e	- the local indices of the elements' vertices
  */
  ierr = DMDAGetElements(da, &nel, &nen, &e);CHKERRQ(ierr);

  /*
     Assemble matrix and RHS
  */
  value[0] = 1.0; value[1] = -1.0; value[2] = -1.0; value[3] = 1.0;
  bvalue[0] = 1.0; bvalue[1] = 1.0;
  for (i=0; i<nel; i++) {
    ierr = MatSetValuesLocal(A, 2, e+nen*i, 2, e+nen*i, value, ADD_VALUES);CHKERRQ(ierr);
    //TODO use VecSetValuesLocal() to set the values
  }
  ierr = DMDARestoreElements(da, &nel, &nen, &e);CHKERRQ(ierr);
  //TODO assemble matrix A
  //TODO assemble vector b
  ierr = VecDuplicate(b,&d);CHKERRQ(ierr);
  ierr = MatGetDiagonal(A, d);CHKERRQ(ierr);
  ierr = VecAbs(d);CHKERRQ(ierr);
  ierr = VecMax(d,NULL,&rho);CHKERRQ(ierr);

  row[0]=0; row[1]=N-1;
  bvalue[0] = 0.0; bvalue[1] = 0.0;
  ierr = VecSetValues(x,2,row,bvalue,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatZeroRowsColumns(A,2,row,rho,x,b); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

  /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following four statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions();
  */
  ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp,nonzeroguess);CHKERRQ(ierr);

  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Clean-up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&d);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary).
  */
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
