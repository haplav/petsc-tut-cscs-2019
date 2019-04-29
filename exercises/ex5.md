# Hands-on 5
```
# usual steps:
cd ~/petsc-tut-cscs-2019/exercises && git pull
make ex5 && mpirun -n 3 ./ex5
```
1. Monitor convergence: `-ksp_monitor`.
2. View solver details: `-ksp_view`.
3. View reason of iteration termination: `-ksp_converged_reason`
   - System matrix is singular.
   - Iterative solver diverges.
   - Direct solver (`-ksp_type preonly –pc_type cholesky`) fails due to a zero pivot.
4. Notice `KSPSetInitialGuessNonzero()` call. Look at its manual page.
5. Enforce Dirichlet boundary conditions using `MatZeroRowsColumns()` (see [below](#dirichlet-boundary-conditions)).
6. Avoid coefficient jumps when enforcing Dirichlet boundary conditions (see [below](#avoid-coefficient-jumps)), compare number of iterations.
7. View solution using `-ksp_view_solution`.
8. The string is fixed on one side – modify code to get it fixed on both sides. (Just one number!)
9. Try various solvers and preconditioners (`-ksp_type`, `-pc_type`),   e.g. `gmres + jacobi`.
10. Bonus: low-level access to direct solvers (see [below](#low-level-access-to-direct-solvers)).


**See also manual pages of `KSPSetFromOptions()`, `KSPSolve()`, `KSPSetType`. Also try `-help` option.**

## Dirichlet Boundary Conditions
* Enforce Dirichlet boundary conditions for DOFs with indices in the index set `dbcidx`. In MATLAB it would be something like
    ```matlab
    diag = 1.0
    for i = dbcidx
      x(i) = 0.0     % or other value of prescribed displacement
      A(i,:) = 0;  A(:,i) = 0;
      A(i,i) = diag;  b(i) = diag*x(i);
    end
    ```
* Then zero the rows using
    ```c
    MatZeroRowsIS(A, dbcidx, diag, x, b);
    MatZeroRowsColumnsIS(A, dbcidx, diag, x, b); // zero also columns, preserving symmetry
    ```

## Avoid coefficient jumps
* We can replace
    ```matlab
    diag = 1.0
    ```
    with the maximum of absolute values of diagonal entries:
    ```matlab
    diag = max(abs(diag(A)));
    ```
* Use `MatGetDiagonal()`, `VecAbs()`, `VecMax()` to get `diag` in PETSc.

## Low-level access to direct solvers
* Paste code snippet `ex5_direct.c` to the end of `ex5.c` (before Clean-up section).
* Look at the output – what does it say?


