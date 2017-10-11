all: main

COMPILER = mpifort
COMPFLAG = -Wall -Ofast

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

objects = properties.o petsc_lib.o particle_props.o particle_lib.o laplace_lib.o

#--------------------------------------------------------------------------
debug: COMPFLAG += -Wall -Wextra -pedantic -g -Og -fimplicit-none -fbacktrace -fcheck=all
debug: main
#--------------------------------------------------------------------------

main: $(objects) main.o chkopts
	-${FLINKER} $(COMPFLAG) -o main $(objects) main.o ${PETSC_KSP_LIB}
