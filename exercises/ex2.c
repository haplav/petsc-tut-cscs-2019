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
  PetscInt       i, myint=5;
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
  /* TODO task 6 */
  ierr = PetscOptionsGetInt(NULL,NULL,"-myint",&myint,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Option -myint set to value \"%D\".\n",myint);CHKERRQ(ierr);
  }

  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND);CHKERRQ(ierr);

  /* TODO tasks 4, 5 - insert below */

  for (i = 0; i < myint; ++i) {
    ierr = PetscViewerASCIIPrintf(viewer, "test line %d\n", i);CHKERRQ(ierr);
  }

  ierr = PetscViewerASCIIPrintf(viewer, "========================\n\n");CHKERRQ(ierr);

  ierr = PetscViewerASCIIPushSynchronized(viewer);CHKERRQ(ierr);
  for (i = rank*myint; i < rank*myint+myint; ++i) {
    ierr = PetscViewerASCIISynchronizedPrintf(viewer, "test line %d\n", i);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIISynchronizedPrintf(viewer, "---\n");CHKERRQ(ierr);
  /* TODO task 7 - add PetscViewerFlush() call */
  ierr = PetscViewerASCIIPopSynchronized(viewer);CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
