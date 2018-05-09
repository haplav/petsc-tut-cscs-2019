/* Adopted from $PETSC_DIR/src/sys/classes/viewer/examples/tutorials/ex1.c */

/*
   Build with
   make ex2
*/

static char help[] = "Appends to an ASCII file.\n\n";

#include <petscviewer.h>

int main(int argc,char **args)
{
  PetscViewer    viewer;
  PetscInt       i, n=5;
  PetscErrorCode ierr;
  char           filename[PETSC_MAX_PATH_LEN] = "";
  PetscBool      flg;
  PetscMPIInt    rank;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  ierr = PetscOptionsGetString(NULL,NULL,"-f",filename,sizeof(filename),&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Option -f set to value \"%s\".\n",filename);CHKERRQ(ierr);
  }
  //TODO task 6
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Option -n set to value \"%d\".\n",n);CHKERRQ(ierr);
  }

  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND);CHKERRQ(ierr);

  //TODO tasks 4, 5

  for (i = 0; i < n; ++i) {
    ierr = PetscViewerASCIIPrintf(viewer, "test line %d\n", i);CHKERRQ(ierr);
  }

  ierr = PetscViewerASCIIPrintf(viewer, "========================\n\n");CHKERRQ(ierr);

  ierr = PetscViewerASCIIPushSynchronized(viewer);CHKERRQ(ierr);
  for (i = rank*n; i < rank*n+n; ++i) {
    ierr = PetscViewerASCIISynchronizedPrintf(viewer, "test line %d\n", i);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIISynchronizedPrintf(viewer, "---\n");CHKERRQ(ierr);
  //TODO task 7
  ierr = PetscViewerASCIIPopSynchronized(viewer);CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
