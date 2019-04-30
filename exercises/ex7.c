/* This does the same as ex6 but using KSPSetCompute* API */

static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

static PetscErrorCode ComputeOperators(KSP ksp, Mat A, Mat B, void *ctx);
static PetscErrorCode ComputeRHS(KSP ksp, Vec b, void *ctx);
static PetscErrorCode ComputeInitialGuess(KSP ksp, Vec b, void *ctx);
static PetscErrorCode FormMatrix(DM da, Mat A);
static PetscErrorCode FormRHS(DM da, Vec b);
static PetscErrorCode FormInitialGuess(DM da, Vec x);

int main(int argc,char **args)
{
  DM             da;
  KSP            ksp;         /* linear solver context */
  PC             pc;          /* preconditioner context */
  PetscErrorCode ierr;
  PetscInt       n = 5, N;
  PetscInt       its;
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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /*
     Associate KSP with DM.
  */
  ierr = KSPSetDM(ksp,da);CHKERRQ(ierr);

  /*
     Set functions that compute the matrix, RHS and initial guess.
  */
  ierr = KSPSetComputeOperators(ksp,ComputeOperators,NULL);CHKERRQ(ierr);
  ierr = KSPSetComputeRHS(ksp,ComputeRHS,NULL);CHKERRQ(ierr);
  if (nonzeroguess) {
    ierr = KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,NULL);CHKERRQ(ierr);
  }

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
  ierr = KSPSolve(ksp,NULL,NULL);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Clean-up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary).
  */
  ierr = PetscFinalize();
  return ierr;
}

static PetscErrorCode FormMatrix(DM da, Mat A)
{
  PetscErrorCode ierr;
  Vec            d;
  PetscInt       i, nel, nen, N;
  const PetscInt *e;
  PetscInt       row[2];
  PetscScalar    rho;
  PetscScalar    value[4] = {1.0, -1.0, -1.0, 1.0};

  PetscFunctionBeginUser;
  ierr = MatGetSize(A,&N,NULL);CHKERRQ(ierr);

  ierr = DMDAGetElements(da, &nel, &nen, &e);CHKERRQ(ierr);
  for (i=0; i<nel; i++) {
    ierr = MatSetValuesLocal(A, 2, e+nen*i, 2, e+nen*i, value, ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = DMDARestoreElements(da, &nel, &nen, &e);CHKERRQ(ierr);

  ierr = MatCreateVecs(A,&d,NULL);CHKERRQ(ierr);
  ierr = MatGetDiagonal(A, d);CHKERRQ(ierr);
  ierr = VecAbs(d);CHKERRQ(ierr);
  ierr = VecMax(d,NULL,&rho);CHKERRQ(ierr);
  ierr = VecDestroy(&d);CHKERRQ(ierr);

  row[0]=0; row[1]=N-1;
  ierr = MatZeroRowsColumns(A,2,row,rho,NULL,NULL); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode FormRHS(DM da, Vec b)
{
  PetscErrorCode ierr;
  PetscInt       i, nel, nen, N;
  const PetscInt *e;
  PetscScalar    bvalue[2] = {1.0, 1.0};
  PetscInt       row[2];

  PetscFunctionBeginUser;
  ierr = VecGetSize(b,&N);CHKERRQ(ierr);

  ierr = DMDAGetElements(da, &nel, &nen, &e);CHKERRQ(ierr);
  for (i=0; i<nel; i++) {
    ierr = VecSetValuesLocal(b, 2, e+nen*i, bvalue, ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = DMDARestoreElements(da, &nel, &nen, &e);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

  /* account for Dirichlet b.c. */
  bvalue[0] = 0.0; bvalue[1] = 0.0;
  row[0]=0; row[1]=N-1;
  ierr = VecSetValues(b,2,row,bvalue,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode FormInitialGuess(DM da, Vec x)
{
  PetscErrorCode ierr;
  PetscScalar p = -1e3;

  PetscFunctionBeginUser;
  ierr = VecSet(x,p);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


static PetscErrorCode ComputeOperators(KSP ksp, Mat A, Mat B, void *ctx)
{
  PetscErrorCode ierr;
  DM dm;

  PetscFunctionBeginUser;
  ierr = KSPGetDM(ksp,&dm);CHKERRQ(ierr);
  ierr = FormMatrix(dm,A);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode ComputeRHS(KSP ksp, Vec b, void *ctx)
{
  PetscErrorCode ierr;
  DM dm;

  PetscFunctionBeginUser;
  ierr = KSPGetDM(ksp,&dm);CHKERRQ(ierr);
  ierr = FormRHS(dm,b);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode ComputeInitialGuess(KSP ksp, Vec x, void *ctx)
{
  PetscErrorCode ierr;
  DM dm;

  PetscFunctionBeginUser;
  ierr = KSPGetDM(ksp,&dm);CHKERRQ(ierr);
  ierr = FormInitialGuess(dm,x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
