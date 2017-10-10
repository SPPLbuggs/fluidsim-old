all: main

COMPILER = mpifort
COMPFLAG = -Wall -Ofast
PETSC_INCL = -I$(PETSC_DIR)/include
PETSC_INCL2 = -I$(PETSC_DIR)/arch-linux2-c-debug/include
LIB_INCL = -L$(PETSC_DIR)/arch-linux2-c-debug/lib

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

objects = properties.o particle_props.o particle_lib.o laplace_lib.o sfc_lib.o

#--------------------------------------------------------------------------
PETSC = $(PETSC_INCL) $(PETSC_INCL2) $(LIB_INCL) -lpetsc
#--------------------------------------------------------------------------
debug: COMPFLAG += -Wall -Wextra -pedantic -g -Og -fimplicit-none -fbacktrace -fcheck=all
debug: main
#--------------------------------------------------------------------------

main: $(objects) main.o
	-${FLINKER} $(COMPFLAG) -o main $(objects) main.o $(PETSC)

%.o: %.F90
	$(COMPILER) $(COMPFLAG) -c $< $(PETSC)

clearscr :
	clear

fresh : | clean clearscr main
