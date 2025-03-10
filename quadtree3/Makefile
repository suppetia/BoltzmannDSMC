#!/bin/bash

MAIN=run_simulation

OBJECTS=

MODDIR=$(shell pwd)/include
LIBDIR=$(shell pwd)/lib
HDF5_LIBDIR="/usr/lib/x86_64-linux-gnu/hdf5/serial"
HDF5_INCLUDEDIR="/usr/include/hdf5/serial"

# Fortran compiler
FC=gfortran
# Fortran compiler flags
CFLAGS:=-I$(MODDIR) -fopenmp

# Linker
LD=gfortran
# Linker flags
LDFLAGS= -L$(LIBDIR) -L$(HDF5_LIBDIR) -lsimulation -lquadtree -lutil -lhdf5 -lhdf5_fortran -lcfgio -fopenmp

# archive utility
AR=ar
ARFLAGS =

# export all variables to sub-Makefiles
export

build: CFLAGS += -O3 -ftree-vectorize -fopt-info-all=optinfo.all
build: LDFLAGS += -O3
build: $(MAIN).out

debug: CFLAGS += -Wall -g -pg -fbounds-check
debug: LDFLAGS += -Wall -pg -fbounds-check
debug: $(MAIN).out


$(MAIN).out: submodules $(OBJECTS) $(MAIN).o
	$(LD) $(OBJECTS) $(MAIN).o -o$(MAIN).out $(LDFLAGS)


submodules: 
	$(MAKE) -C util
	$(MAKE) -C quadtree
	$(MAKE) -C simulation

%.o: %.f90
	$(FC) -c $(CFLAGS) $< -o $@ -J$(MODDIR)

clean:
	rm -f $(MAIN).out *.o *.mod
	# cd $(MODDIR) && rm $(patsubst %.o,m_%.mod, $(OBJECTS))
	$(MAKE) -C util clean
	$(MAKE) -C quadtree clean
	$(MAKE) -C simulation clean

