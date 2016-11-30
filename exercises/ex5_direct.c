
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Low-level access to direct solver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  {
    PC pc;
    Mat L;
    Vec x1;
    PetscReal norm;

    ierr = VecDuplicate(x,&x1);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCCHOLESKY);CHKERRQ(ierr);
    ierr = PCSetUp(pc);CHKERRQ(ierr);
    ierr = PCFactorGetMatrix(pc,&L);CHKERRQ(ierr);
    ierr = MatSolve(L,b,x1);CHKERRQ(ierr);

    ierr = VecAXPY(x1,-1.0,x);CHKERRQ(ierr);
    ierr = VecNorm(x1,NORM_2,&norm);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"||x-x1||=%g\n",norm);CHKERRQ(ierr);

    ierr = VecDestroy(&x1);CHKERRQ(ierr);
  }
