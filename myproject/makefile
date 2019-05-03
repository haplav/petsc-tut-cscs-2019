CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
CLEANFILES       = rhs.vtk solution.vtk
APPNAME          = myproj

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Your application-specific targets can go here.
# Note that programs starting with ex% (like all PETSc provided examples)
# don't need a specific target if you include ${PETSC_DIR}/lib/petsc/conf/test

${APPNAME}: ${APPNAME}.o | chkopts
	-${CLINKER} -o $@ $^  ${PETSC_LIB}
	${DSYMUTIL} $@

include ${PETSC_DIR}/lib/petsc/conf/test
