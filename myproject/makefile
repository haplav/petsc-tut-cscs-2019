CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
CLEANFILES       = rhs.vtk solution.vtk
APPNAME          = myproj

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Your application-specific targets can go here.

${APPNAME}: ${APPNAME}.o | chkopts
	-${CLINKER} -o $@ $^  ${PETSC_LIB}
	${DSYMUTIL} $@

include ${PETSC_DIR}/lib/petsc/conf/test
