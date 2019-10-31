
static char help[] = "VecScatter example. Invocation: ./ex_vecscatter -view_{u,v}0 -view_{iu,iv,u,v,sc}1 -view_{iu,iv,u,v,sc}2\n\n";

#include <petscvec.h>

int main(int argc,char **args)
{
  Vec u,v;
  IS iu,iv;
  PetscInt *iua, *iva;
  PetscInt i;
  VecScatter sc;
  PetscInt ulo,vhi,n,N,l=2;
  MPI_Comm comm;
  PetscRandom rnd;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  comm = PETSC_COMM_WORLD;
  ierr = PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL);CHKERRQ(ierr);

  ierr = PetscRandomCreate(comm,&rnd);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(rnd, 0.0, 100.0);CHKERRQ(ierr);

  ierr = VecCreateMPI(comm, 3, PETSC_DECIDE, &u);CHKERRQ(ierr);
  ierr = VecCreateMPI(comm, 4, PETSC_DECIDE, &v);CHKERRQ(ierr);
  ierr = VecSetRandom(u, rnd);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(u,&ulo,NULL);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(v,NULL,&vhi);CHKERRQ(ierr);
  ierr = VecGetLocalSize(v,&n);CHKERRQ(ierr);
  ierr = VecGetSize(v,&N);CHKERRQ(ierr);
  ierr = PetscMalloc2(l,&iua,l,&iva);CHKERRQ(ierr);
  ierr = VecSet(v, 0.0);CHKERRQ(ierr);
  ierr = VecViewFromOptions(u,NULL,"-view_u0");CHKERRQ(ierr);
  ierr = VecViewFromOptions(v,NULL,"-view_v0");CHKERRQ(ierr);

  for (i=0; i<l; i++) {
    iua[i]=ulo+i;
    iva[i]=vhi-1-i;
  }
  ierr = ISCreateGeneral(comm,l,iua,PETSC_USE_POINTER,&iu);CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm,l,iva,PETSC_USE_POINTER,&iv);CHKERRQ(ierr);
  ierr = ISViewFromOptions(iu,NULL,"-view_iu1");CHKERRQ(ierr);
  ierr = ISViewFromOptions(iv,NULL,"-view_iv1");CHKERRQ(ierr);
  ierr = VecScatterCreate(u,iu,v,iv,&sc);CHKERRQ(ierr);
  ierr = VecScatterViewFromOptions(sc,NULL,"-view_sc1");CHKERRQ(ierr);
  ierr = VecScatterBegin(sc,u,v,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(  sc,u,v,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&sc);CHKERRQ(ierr);
  ierr = ISDestroy(&iu);CHKERRQ(ierr);
  ierr = ISDestroy(&iv);CHKERRQ(ierr);
  ierr = VecViewFromOptions(u,NULL,"-view_u1");CHKERRQ(ierr);
  ierr = VecViewFromOptions(v,NULL,"-view_v1");CHKERRQ(ierr);

  ierr = VecSet(v, 0.0);CHKERRQ(ierr);
  for (i=0; i<l; i++) {
    iua[i]=i;
    iva[i]=N-1-i;
  }
  ierr = ISCreateGeneral(comm,l,iua,PETSC_USE_POINTER,&iu);CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm,l,iva,PETSC_USE_POINTER,&iv);CHKERRQ(ierr);
  ierr = ISViewFromOptions(iu,NULL,"-view_iu2");CHKERRQ(ierr);
  ierr = ISViewFromOptions(iv,NULL,"-view_iv2");CHKERRQ(ierr);
  ierr = VecScatterCreate(u,iu,v,iv,&sc);CHKERRQ(ierr);
  ierr = VecScatterViewFromOptions(sc,NULL,"-view_sc2");CHKERRQ(ierr);
  ierr = VecScatterBegin(sc,u,v,ADD_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(  sc,u,v,ADD_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&sc);CHKERRQ(ierr);
  ierr = ISDestroy(&iu);CHKERRQ(ierr);
  ierr = ISDestroy(&iv);CHKERRQ(ierr);
  ierr = VecViewFromOptions(u,NULL,"-view_u2");CHKERRQ(ierr);
  ierr = VecViewFromOptions(v,NULL,"-view_v2");CHKERRQ(ierr);

  ierr = PetscRandomDestroy(&rnd);CHKERRQ(ierr);
  ierr = PetscFree2(iua,iva);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

