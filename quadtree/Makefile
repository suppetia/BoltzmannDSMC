#!/bin/bash

MAIN=test_quadtree

OBJECTS=types.o quadtree.o quadtree_io.o matrix_io.o $(MAIN).o

# Fortran compiler
FC=gfortran
# Fortran compiler flags
CFLAGS=

# Linker
LD=gfortran
# Linker flags
LDFLAGS=


build: CFLAGS += -O2
build: LDFLAGS += -O2
build: $(MAIN).out


debug: CFLAGS += -Wall -g
debug: LDFLAGS += -Wall -g
debug: $(MAIN).out


# build the final executable
$(MAIN).out: $(OBJECTS)
	$(LD) $(OBJECTS) -o$(MAIN).out $(LDFLAGS)
	
%.o: %.f90
	$(FC) -c $(CFLAGS) $< -o $@


clean:
	rm -f *.o
	rm -f *.out
	rm -f *.mod
