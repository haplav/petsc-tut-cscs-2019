
ALL:

CFLAGS	        =
FFLAGS	        = -I${PETSC_DIR}/include/petsc/finclude
CPPFLAGS        =
FPPFLAGS        =
EXAMPLESC       = ex6.c
EXAMPLESF       = 

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

ex6: ex6.o   chkopts
	-${CLINKER} -o ex6 ex6.o  ${PETSC_KSP_LIB} ${SLEPC_EPS_LIB}
	${RM} -f ex6.o
