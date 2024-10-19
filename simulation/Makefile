#!/bin/bash

MAIN=test_simulation

OBJECTS=simulation.o

MODDIR=$(shell pwd)/include
LIBDIR=$(shell pwd)/lib
HDF5_LIBDIR="/usr/lib/x86_64-linux-gnu/hdf5/serial"
HDF5_INCLUDEDIR="/usr/include/hdf5/serial"

# Fortran compiler
FC=gfortran
# Fortran compiler flags
CFLAGS:=-I$(MODDIR)

# Linker
LD=gfortran
# Linker flags
LDFLAGS= -L$(LIBDIR) -L$(HDF5_LIBDIR) -lquadtree -lutil -lhdf5 -lhdf5_fortran

# archive utility
AR=ar
ARFLAGS =

# export all variables to sub-Makefiles
export

build: CFLAGS += -O3 -ftree-vectorize
build: LDFLAGS += -O3
build: $(MAIN).out

debug: CFLAGS += -Wall -g -fbounds-check
debug: LDFLAGS += -Wall -g -fbounds-check
debug: $(MAIN).out


$(MAIN).out: submodules $(OBJECTS) $(MAIN).o
	$(LD) $(OBJECTS) $(MAIN).o -o$(MAIN).out $(LDFLAGS)


submodules: 
	$(MAKE) -C util
	$(MAKE) -C quadtree

%.o: %.f90
	$(FC) -c $(CFLAGS) $< -o $@ -J$(MODDIR)

clean:
	rm -f $(MAIN).out *.o *.mod
	rm -f $(MODDIR)/*.mod
	$(MAKE) -C util clean
	$(MAKE) -C quadtree clean

