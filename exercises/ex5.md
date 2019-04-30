# Hands-on 5
```
# usual steps:
cd ~/petsc-tut-cscs-2019/exercises && git pull
make ex5 && mpirun -n 3 ./ex5
```
1. Monitor convergence: `-ksp_monitor`.
2. View solver details: `-ksp_view`.
3. View reason of iteration termination with `-ksp_converged_reason` and final residual with `-ksp_final_residual`.
   - System matrix is singular indefinite.
   - GMRES iterative solver (`-ksp_type gmres`; this is a default) diverges with `DIVERGED_ITS` reason (maximum iteration count reached).
   - CG iterative solver (`-ksp_type cg`) diverges with `DIVERGED_INDEFINITE_MAT` reason.
   - PETSc direct solver (`-ksp_type preonly –pc_type cholesky -pc_factor_mat_solver_type petsc`; sequential only, so run with `mpirun -n 1`) fails due to a zero pivot.
   - MUMPS direct solver (`-ksp_type preonly –pc_type cholesky -pc_factor_mat_solver_type mumps`) works because MUMPS supports singular matrix factorizations. However, norm of solution is very large which makes us scared.
   - LSQR iterative solver (`-ksp_type lsqr`) converges to a good-looking least-square solution within a few iterations.
4. Enforce Dirichlet boundary conditions using `MatZeroRowsColumns()` (see [below](#dirichlet-boundary-conditions)). For now, fix just one side of the string (`ndbc=1`). Run the solvers above again.
5. Look the code section `Analyze the results` and its output. To have some information about the matrix as well, print `max(abs(diag(A)))` and `min(abs(diag(A)))`.
   - Hint: Use `VecDuplicate()`, `MatGetDiagonal()`, `VecAbs()`, `VecMax()` in this order.
6. View solution using `-ksp_view_solution`. Try different problem sizes using `-n`.
7. The string is fixed on one side – modify code to get it fixed on both sides. (Just one number!)
8. Play around with tolerances (`-ksp_rtol`, `-ksp_atol`, `-ksp_max_it`). Try various solvers and preconditioners (`-ksp_type`, `-pc_type`),   e.g. `gmres + jacobi`.
9. Notice that all values of `x` are initially set to 1000. However, initial residual norm displayed in the monitor output (`-ksp_monitor`) does not correspond to that.
   - This is because KSP ignores the initial guess by default and uses zero vector.
   - This can be changed with `-ksp_initial_guess_nonzero` option or `KSPSetInitialGuessNonzero()` function.
   - Use `-ksp_type gmres -n 16` and compare `-ksp_initial_guess_nonzero {0,1}` - monitor output, number of iterations, and `-ksp_view`.
10. Bonus: low-level access to direct solvers (see [below](#low-level-access-to-direct-solvers)).


**See also manual pages of `KSPSetFromOptions()`, `KSPSolve()`, `KSPSetType`. Also try `-help` option.**

## Dirichlet Boundary Conditions
* Enforce Dirichlet boundary conditions for DOFs with indices in the index set `dbcidx`. In MATLAB it would be something like
    ```matlab
    rho = 1.0  % or some other value to avoid jumps, e.g. max(abs(diag(A))) or maxeig(A)
    for i = dbcidx
      x(i) = 0.0     % or other value of prescribed displacement
      A(i,:) = 0;  A(:,i) = 0;
      A(i,i) = rho;  b(i) = rho*x(i);
    end
    ```
* Then zero the rows using
    ```c
    MatZeroRows(A, ndbc, dbcidx, rho, x, b); /* ndbc is number of zeroed rows */
    MatZeroRowsColumns(A, ndbc, dbcidx, rho, x, b); /* zero also columns, preserving symmetry */
    ```
    or, if you have `dbcidx` as `IS`,
    ```c
    MatZeroRowsIS(A, dbcidx, rho, x, b);
    MatZeroRowsColumnsIS(A, dbcidx, rho, x, b); /* zero also columns, preserving symmetry */
    ```

## Low-level access to direct solvers
* Paste code snippet `ex5_direct.c` to the end of `ex5.c` (before Clean-up section).
* Look at the output – what does it say?


