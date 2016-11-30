/* This does the same as ex7 but using SNES */

static char help[] = "Solves a tridiagonal linear system as a nonlinear problem using SNES.\n\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>

static PetscErrorCode ComputeJacobian(SNES snes, Vec x, Mat A, Mat B, void *ctx);
static PetscErrorCode ComputeFunction(SNES snes, Vec x, Vec f, void *ctx);
static PetscErrorCode FormMatrix(DM da, Mat A);
static PetscErrorCode FormRHS(DM da, Vec b);

#undef __FUNCT__
#define __FUNCT__ "main"
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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the nonlinear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /*
     Associate SNES with DM.
  */
  ierr = SNESSetDM(snes,da);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,NULL,NULL,ComputeJacobian,NULL);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,NULL,ComputeFunction,NULL);CHKERRQ(ierr);

  ierr = SNESSetType(snes,SNESNEWTONLS);CHKERRQ(ierr);
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESSolve(snes,NULL,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Finalize
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
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
#define __FUNCT__ "ComputeJacobian"
static PetscErrorCode ComputeJacobian(SNES snes, Vec x, Mat A, Mat B, void *ctx)
{
  PetscErrorCode ierr;
  DM dm;

  PetscFunctionBeginUser;
  ierr = SNESGetDM(snes,&dm);CHKERRQ(ierr);
  ierr = FormMatrix(dm,A);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeFunction"
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
