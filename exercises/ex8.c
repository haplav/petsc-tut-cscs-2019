/* This does the same as ex7 but using SNES */

static char help[] = "Solves a tridiagonal linear system as a nonlinear problem using SNES.\n\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>

static PetscErrorCode ComputeJacobian(SNES snes, Vec x, Mat A, Mat B, void *ctx);
static PetscErrorCode ComputeFunction(SNES snes, Vec x, Vec f, void *ctx);
static PetscErrorCode FormMatrix(DM da, Mat A);
static PetscErrorCode FormRHS(DM da, Vec b);

int main(int argc,char **args)
{
  DM             da;
  SNES           snes;        /* nonlinear solver context */
  KSP            ksp;         /* linear solver context */
  PC             pc;          /* preconditioner context */
  PetscErrorCode ierr;
  PetscInt       n = 5, N;

  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

  /*
     Create DMDA context for structured 1D problem.
  */
  ierr = MPI_Allreduce(&n, &N, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,N,1,1,NULL,&da);CHKERRQ(ierr);
  ierr = DMSetFromOptions(da);CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the nonlinear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create nonlinear solver context
  */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /*
     Associate SNES with DM.
  */
  ierr = SNESSetDM(snes,da);CHKERRQ(ierr);

  /*
     Set functions that compute the matrix and RHS.
  */
  ierr = SNESSetJacobian(snes,NULL,NULL,ComputeJacobian,NULL);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,NULL,ComputeFunction,NULL);CHKERRQ(ierr);

  /*
     Set linear and nonlinear solver defaults for this problem (optional).
  */
  ierr = SNESSetType(snes,SNESNEWTONLS);CHKERRQ(ierr);
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);

  /*
    Set runtime options.
  */
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESSolve(snes,NULL,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Clean-up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
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

static PetscErrorCode ComputeJacobian(SNES snes, Vec x, Mat A, Mat B, void *ctx)
{
  PetscErrorCode ierr;
  DM dm;

  PetscFunctionBeginUser;
  ierr = SNESGetDM(snes,&dm);CHKERRQ(ierr);
  ierr = FormMatrix(dm,A);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode ComputeFunction(SNES snes, Vec x, Vec f, void *ctx)
{
  PetscErrorCode ierr;
  DM dm;
  Mat A;

  PetscFunctionBeginUser;
  ierr = SNESGetDM(snes,&dm);CHKERRQ(ierr);

  /* f = -b */
  ierr = FormRHS(dm,f);CHKERRQ(ierr);
  ierr = VecScale(f,-1.0);CHKERRQ(ierr);

  /* create and form the matrix if it has not been created yet */
  ierr = DMGetApplicationContext(dm,&A);CHKERRQ(ierr);
  if (!A) {
    ierr = DMSetMatType(dm,MATAIJ);CHKERRQ(ierr);
    ierr = DMCreateMatrix(dm,&A);CHKERRQ(ierr);
    ierr = FormMatrix(dm,A);CHKERRQ(ierr);
    ierr = DMSetApplicationContext(dm,A);CHKERRQ(ierr);
    ierr = DMSetApplicationContextDestroy(dm,(PetscErrorCode (*)(void**))MatDestroy);CHKERRQ(ierr);
  }

  /* f = A*x - b */
  ierr = MatMultAdd(A,x,f,f);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
