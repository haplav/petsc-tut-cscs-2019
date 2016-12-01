
ALL:

CFLAGS	        =
FFLAGS	        = -I${PETSC_DIR}/include/petsc/finclude
CPPFLAGS        =
FPPFLAGS        =
EXAMPLESC       = ex1.c
EXAMPLESF       = 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

ex1: ex1.o   chkopts
	-${CLINKER} -o ex1 ex1.o  ${PETSC_SYS_LIB}
	${RM} -f ex1.o

ex2: ex2.o   chkopts
	-${CLINKER} -o ex2 ex2.o  ${PETSC_SYS_LIB}
	${RM} -f ex2.o

ex3: ex3.o   chkopts
	-${CLINKER} -o ex3 ex3.o  ${PETSC_VEC_LIB}
	${RM} -f ex3.o

ex4: ex4.o   chkopts
	-${CLINKER} -o ex4 ex4.o  ${PETSC_KSP_LIB}
	${RM} -f ex4.o

ex5: ex5.o   chkopts
	-${CLINKER} -o ex5 ex5.o  ${PETSC_KSP_LIB}
	${RM} -f ex5.o

ex6: ex6.o   chkopts
	-${CLINKER} -o ex6 ex6.o  ${PETSC_SNES_LIB}
	${RM} -f ex6.o

ex7: ex7.o   chkopts
	-${CLINKER} -o ex7 ex7.o  ${PETSC_SNES_LIB}
	${RM} -f ex7.o

ex8: ex8.o   chkopts
	-${CLINKER} -o ex8 ex8.o  ${PETSC_SNES_LIB}
	${RM} -f ex8.o

ex9: ex9.o   chkopts
	-${CLINKER} -o ex9 ex9.o  ${PETSC_SNES_LIB}
	${RM} -f ex9.o

ex10: ex10.o   chkopts
	-${CLINKER} -o ex10 ex10.o  ${PETSC_TAO_LIB}
	${RM} -f ex10.o
