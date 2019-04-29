
static char help[] = "Empty program.\n\n";

#include <petscksp.h>

int main(int argc,char **args)
{
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscFinalize();
  return ierr;
}

