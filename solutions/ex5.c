
static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/
#include <petscksp.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            x, b;        /* approx solution, RHS */
  Vec            d;
  Mat            A;           /* linear system matrix */
  KSP            ksp;         /* linear solver context */
  PC             pc;          /* preconditioner context */
  PetscErrorCode ierr;
  PetscInt       i, n = 5, N, Istart, Iend;
  PetscInt       its;
  PetscInt       row[2], col[2], dbcidx[2];
  PetscScalar    value[4], bvalue[2];
  PetscReal      norm = PETSC_INFINITY;

  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create vectors.  Note that we form 1 vector from scratch and
     then duplicate as needed.
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x,n,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

  /* Set the same value to all vector entries.  */
  ierr = VecSet(b,0.0);CHKERRQ(ierr);
  //TODO task 4
  /* Set nonzero initial guess. It is used only if -ksp_initial_guess_nonzero 1, otherwise ignored. */
  ierr = VecSet(x,1000.0);CHKERRQ(ierr);
    
  /*
     Create matrix.  When using MatCreate(), the matrix format can
     be specified at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.
  */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,n,n,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetSize(A,&N,NULL);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

  /*
     Assemble matrix
  */
  value[0] = 1.0; value[1] = -1.0; value[2] = -1.0; value[3] = 1.0;
  bvalue[0] = 1.0; bvalue[1] = 1.0;
  if (Istart == 0) Istart = 1;
  for (i=Istart; i<Iend; i++) {
    row[0] = i-1; row[1] = i;
    col[0] = i-1; col[1] = i;
    ierr   = MatSetValues(A,2,row,2,col,value,ADD_VALUES);CHKERRQ(ierr);
    ierr   = VecSetValues(b,2,row,bvalue,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

  dbcidx[0]=0; dbcidx[1]=N-1;
  bvalue[0] = 0.0; bvalue[1] = 0.0;
  ierr = VecSetValues(x,2,dbcidx,bvalue,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
  /* task 5 */
  /* ierr = MatZeroRowsColumns(A,1,dbcidx,1.0,x,b);CHKERRQ(ierr); */
  /* task 8 */
  ierr = MatZeroRowsColumns(A,2,dbcidx,1.0,x,b);CHKERRQ(ierr);

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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Analyze the results
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* Get iteration count */
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of iterations: %D\n",its);CHKERRQ(ierr);
  //TODO task 6 - use MatGetDiagonal, VecAbs, VecMax; you can reuse PetscReal norm
  ierr = VecDuplicate(b,&d);CHKERRQ(ierr);
  ierr = MatGetDiagonal(A, d);CHKERRQ(ierr);
  ierr = VecAbs(d);CHKERRQ(ierr);
  ierr = VecMax(d,NULL,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"max(abs(diag(A))) = %.2e\n",norm);CHKERRQ(ierr);
  ierr = VecMin(d,NULL,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"min(abs(diag(A))) = %.2e\n",norm);CHKERRQ(ierr);
  /* Compute ||x|| */
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"||x|| = %.2e\n",norm);CHKERRQ(ierr);
  /* Compute ||b|| */
  ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"||b|| = %.2e\n",norm);CHKERRQ(ierr);
  /* Compute ||A*x-b|| */
  ierr = MatMult(A,x,d);CHKERRQ(ierr);
  ierr = VecAXPY(d,-1.0,b);CHKERRQ(ierr);
  ierr = VecNorm(d,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"||A*x-b|| = %.2e\n",norm);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Low-level access to direct solver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  {
    PC pc;
    Mat L;
    Vec x1;
    PetscReal norm;

    ierr = VecDuplicate(x,&x1);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCCHOLESKY);CHKERRQ(ierr);
    ierr = PCSetUp(pc);CHKERRQ(ierr);
    ierr = PCFactorGetMatrix(pc,&L);CHKERRQ(ierr);
    ierr = MatSolve(L,b,x1);CHKERRQ(ierr);

    ierr = VecAXPY(x1,-1.0,x);CHKERRQ(ierr);
    ierr = VecNorm(x1,NORM_2,&norm);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"||x-x1||=%g\n",norm);CHKERRQ(ierr);

    ierr = VecDestroy(&x1);CHKERRQ(ierr);
  }

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

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary).
  */
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
