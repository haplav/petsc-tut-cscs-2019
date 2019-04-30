1.  `cd ~/petsc-tut-cscs-2019/exercises && git pull`
2.  `make ex4 && srun -n 4 ./ex4`
3.  Look at the code for matrix assembly.
4.  Add "element contributions" to vector `b` using `VecSetValues()`.
    - Add values `[1;1]` to positions `[0;1], ..., [n-2;n-1]`.
    - You can reuse the loop with `MatSetValues()`.  
5.  Try various views of the matrix and vectors
    - programmatically - see documentation of `PetscViewerPushFormat()`,
    - by option (`-vec_view`, `-mat_view {,::ascii_dense, ::ascii_info, ::ascii_info_detail, ::ascii_matlab}`).
6.  Change type of matrix to dense
    - programmatically (see `MatSetType`, `MATDENSE`),
    - by option (`-mat_type dense`).
7.  Compute `x = A * b`  (`MatMult()`).
8.  Compute `x = 2 * A * b` (`MatScale` or `VecScale`) and check the result `[-2;2;0;2;-2;]`.
