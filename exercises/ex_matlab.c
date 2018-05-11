#include <petscmat.h>
#include <petscmatlab.h>

int main(int argc,char **args)
{
  Mat            A;
  Vec            x,x_matlab,b;
  PetscMatlabEngine meng;
  MPI_Comm       comm;
  PetscReal      norm;
  PetscScalar    n;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc,&args,(char*)0,(char*)0);if (ierr) return ierr;
  A = NULL; b = NULL; x = NULL; x_matlab = NULL;
  comm = PETSC_COMM_WORLD;
  meng = PETSC_MATLAB_ENGINE_(comm);

  ierr = PetscMatlabEngineEvaluate(meng,"mscript;");CHKERRQ(ierr);
  ierr = PetscMatlabEnginePrintOutput(meng,PETSC_STDOUT);CHKERRQ(ierr);
  ierr = PetscMatlabEngineEvaluate(meng,"A = mfunction(A)");CHKERRQ(ierr);
  ierr = PetscMatlabEnginePrintOutput(meng,PETSC_STDOUT);CHKERRQ(ierr);
  ierr = PetscMatlabEngineEvaluate(meng,"b = [1; 2; 3]");CHKERRQ(ierr);
  ierr = PetscMatlabEnginePrintOutput(meng,PETSC_STDOUT);CHKERRQ(ierr);
  ierr = PetscMatlabEngineEvaluate(meng,"x = A*b;");CHKERRQ(ierr);
  ierr = PetscMatlabEnginePrintOutput(meng,PETSC_STDOUT);CHKERRQ(ierr);

  ierr = PetscMatlabEngineEvaluate(meng,"n = nnz(A)");CHKERRQ(ierr);
  ierr = PetscMatlabEnginePrintOutput(meng,PETSC_STDOUT);CHKERRQ(ierr);
  ierr = PetscMatlabEngineGetArray(meng,1,1,&n,"n");CHKERRQ(ierr);
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,3,3,(PetscInt)n,NULL,&A);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)A,"A");CHKERRQ(ierr);
  ierr = PetscMatlabEngineGet(meng,(PetscObject)A);CHKERRQ(ierr);

  ierr = VecCreateSeq(PETSC_COMM_SELF,3,&b);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)b,"b");CHKERRQ(ierr);
  ierr = PetscMatlabEngineGet(meng,(PetscObject)b);CHKERRQ(ierr);

  ierr = VecCreateSeq(PETSC_COMM_SELF,3,&x_matlab);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)x_matlab,"x");CHKERRQ(ierr);
  ierr = PetscMatlabEngineGet(meng,(PetscObject)x_matlab);CHKERRQ(ierr);

  ierr = VecDuplicate(x_matlab,&x);CHKERRQ(ierr);
  ierr = MatMult(A,b,x);CHKERRQ(ierr);

  ierr = VecAXPY(x_matlab,-1.0,x);CHKERRQ(ierr);
  ierr = VecNorm(x_matlab,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"\n||x-x_matlab||=%.8f\n",norm);CHKERRQ(ierr);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&x_matlab);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
