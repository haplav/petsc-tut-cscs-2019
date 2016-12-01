/* This does the same as ex7 but using SNES */

static char help[] = "Solves a tridiagonal linear system as a nonlinear problem using SNES.\n\n";

#include <fllopqps.h>
#include <petscdm.h>
#include <petscdmda.h>

static PetscErrorCode FormMatrix(DM da, Mat A);
static PetscErrorCode FormRHS(DM da, Vec b);

typedef struct {

  DM dm; /* distributed array data structure */
  Mat H; /* Quadratic Objective term */
  Vec B; /* Linear Objective term */

} AppCtx;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Tao		 tao;
  PetscErrorCode ierr;
  PetscInt       n = 5, N;
  AppCtx         user;
  PetscScalar    ub=1.75;
  Vec		 x,xl,xu,xqp;
  QP             qp;
  QPS            qps;
  

  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

  /*
     Create DMDA context for structured 1D problem.
  */
  ierr = MPI_Allreduce(&n, &N, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,N,1,1,NULL,&user.dm);CHKERRQ(ierr);

  /*
     Create vectors.  We get 1 vector from DM and
     then duplicate as needed.
  */
  ierr = DMCreateGlobalVector(user.dm, &x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.B);CHKERRQ(ierr);
  ierr = FormRHS(user.dm, user.B);CHKERRQ(ierr);
  /*
     Set the same value to all vector entries.
  */
  ierr = VecSet(x,0.0);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&xl);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&xu);CHKERRQ(ierr);  
  ierr = VecSet(xl, 0);CHKERRQ(ierr);
  ierr = VecSet(xu, ub);CHKERRQ(ierr);

  /*
     Get preallocated matrix from DM. 

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.
  */
  ierr = DMCreateMatrix(user.dm,&user.H);CHKERRQ(ierr);
  ierr = MatSetFromOptions(user.H);CHKERRQ(ierr);
  ierr = FormMatrix(user.dm, user.H);CHKERRQ(ierr);
  
  /* Prescribe the QP problem. */
  ierr = QPCreate(PETSC_COMM_WORLD, &qp);CHKERRQ(ierr);
  //ierr = QPSetOperator(qp, user.H, QP_SYM_SYMMETRIC);CHKERRQ(ierr);
  ierr = QPSetOperator(qp, user.H);CHKERRQ(ierr);
  ierr = QPSetRhs(qp, user.B);CHKERRQ(ierr); 
  ierr = QPSetBox(qp, xl, xu);CHKERRQ(ierr);
  
  /* Create the QP solver (QPS). */
  ierr = QPSCreate(PETSC_COMM_WORLD, &qps);CHKERRQ(ierr);
  
  /* Create the TAO solver*/
  //ierr = QPSSetType(qps,QPSTAO);CHKERRQ(ierr);

  /* Insert the QP problem into the solver. */
  ierr = QPSSetQP(qps, qp);CHKERRQ(ierr);
  
  
  /* Set the QPS monitor */
  ierr = QPSMonitorSet(qps,QPSMonitorDefault,NULL,0);CHKERRQ(ierr);
  
  /* Set QPS options from the options database (overriding the defaults). */
  ierr = QPSSetFromOptions(qps);CHKERRQ(ierr);
  
  /* Solve the QP */
  ierr = QPSSolve(qps);CHKERRQ(ierr);
  
  /* Get the solution vector */
  ierr = QPGetSolutionVector(qp, &xqp);CHKERRQ(ierr);
  
  ierr = VecView(xqp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /* Free PETSc data structures */
  ierr = QPSDestroy(&qps);CHKERRQ(ierr);
  ierr = QPDestroy(&qp);CHKERRQ(ierr);
  ierr = VecDestroy(&xl);CHKERRQ(ierr);
  ierr = VecDestroy(&xu);CHKERRQ(ierr);
  ierr = MatDestroy(&user.H);CHKERRQ(ierr);
  ierr = VecDestroy(&user.B);CHKERRQ(ierr);
  /* Free TAO data structures */

  ierr = DMDestroy(&user.dm);CHKERRQ(ierr);
  
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
  PetscBool      assembled;

  PetscFunctionBeginUser;
  ierr = MatAssembled(A,&assembled);CHKERRQ(ierr);
  if (assembled){ierr = MatZeroEntries(A);CHKERRQ(ierr);}

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

