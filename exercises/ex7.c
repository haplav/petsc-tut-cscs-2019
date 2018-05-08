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

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  DM             da;
  KSP            ksp;         /* linear solver context */
  PC             pc;          /* preconditioner context */
  PetscErrorCode ierr;
  PetscInt       n = 5, N;
  PetscBool      nonzeroguess = PETSC_TRUE;

  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-nonzero_guess",&nonzeroguess,NULL);CHKERRQ(ierr);


  /*
     Create DMDA context for structured 1D problem.
  */
  ierr = MPI_Allreduce(&n, &N, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,N,1,1,NULL,&da);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /*
     Associate KSP with DM.
  */
  ierr = KSPSetDM(ksp,da);CHKERRQ(ierr);
  ierr = KSPSetComputeOperators(ksp,ComputeOperators,NULL);CHKERRQ(ierr);
  ierr = KSPSetComputeRHS(ksp,ComputeRHS,NULL);CHKERRQ(ierr);
  if (nonzeroguess) {
    ierr = KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,NULL);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  }


  ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);

  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = KSPSolve(ksp,NULL,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Finalize
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "FormMatrix"
static PetscErrorCode FormMatrix(DM da, Mat A)
{
  PetscErrorCode ierr;
  Vec            d;
  PetscInt       i, nel, nen, N;
  const PetscInt *e;
  PetscInt       row[2];
  PetscScalar    rho;
  PetscScalar    value[4];

  PetscFunctionBeginUser;
  ierr = MatGetSize(A,&N,NULL);CHKERRQ(ierr);

  ierr = DMDAGetElements(da, &nel, &nen, &e);CHKERRQ(ierr);
 
  value[0] = 1.0; value[1] = -1.0; value[2] = -1.0; value[3] = 1.0;
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
  //ierr = MatZeroRowsColumns(A,2,row,rho,x,b); CHKERRQ(ierr);
  ierr = MatZeroRowsColumns(A,2,row,rho,NULL,NULL); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormRHS"
static PetscErrorCode FormRHS(DM da, Vec b)
{
  PetscErrorCode ierr;
  PetscInt N;

  PetscFunctionBeginUser;
  ierr = VecGetSize(b,&N);CHKERRQ(ierr);
  ierr = VecSet(b,1.0);CHKERRQ(ierr);

  /* take into account Dirichlet b.c. */
  {
    PetscInt row[2] = {0, N-1};
    PetscScalar val[2] = {0.0, 0.0};

    ierr = VecSetValues(b,2,row,val,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormInitialGuess"
static PetscErrorCode FormInitialGuess(DM da, Vec x)
{
  PetscErrorCode ierr;
  PetscScalar p = -1e3;

  PetscFunctionBeginUser;
  ierr = VecSet(x,p);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeOperators"
static PetscErrorCode ComputeOperators(KSP ksp, Mat A, Mat B, void *ctx)
{
  PetscErrorCode ierr;
  DM dm;

  PetscFunctionBeginUser;
  ierr = KSPGetDM(ksp,&dm);CHKERRQ(ierr);
  ierr = FormMatrix(dm,A);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
static PetscErrorCode ComputeRHS(KSP ksp, Vec b, void *ctx)
{
  PetscErrorCode ierr;
  DM dm;

  PetscFunctionBeginUser;
  ierr = KSPGetDM(ksp,&dm);CHKERRQ(ierr);
  ierr = FormRHS(dm,b);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeInitialGuess"
static PetscErrorCode ComputeInitialGuess(KSP ksp, Vec x, void *ctx)
{
  PetscErrorCode ierr;
  DM dm;

  PetscFunctionBeginUser;
  ierr = KSPGetDM(ksp,&dm);CHKERRQ(ierr);
  ierr = FormInitialGuess(dm,x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
