#!/bin/bash

MAIN=test_quadtree
PYTHON_INTERFACE=quadtree_io_py
JULIA_INTERFACE=quadtree_io_jl

OBJECTS=types.o quadtree.o quadtree_io.o matrix_io.o $(MAIN).o 

HDF5_LIBDIR="/usr/lib/x86_64-linux-gnu/hdf5/serial"
HDF5_INCLUDEDIR="/usr/include/hdf5/serial"


# Fortran compiler
FC=gfortran
# Fortran compiler flags
CFLAGS= -fPIC -I$(HDF5_INCLUDEDIR)

# Linker
LD=gfortran
# Linker flags
LDFLAGS= -L$(HDF5_LIBDIR) -lhdf5 -lhdf5_fortran

# Python installation
PY=$(shell which python3)


build: CFLAGS += -O2
build: LDFLAGS += -O2
build: $(MAIN).out


debug: CFLAGS += -Wall -g -fbounds-check
debug: LDFLAGS += -Wall -g -fbounds-check
debug: $(MAIN).out


julia: $(OBJECTS) $(JULIA_INTERFACE).o
	# create a shared library
	$(FC) -shared -o $(JULIA_INTERFACE).so $(OBJECTS) $(JULIA_INTERFACE).o

python: $(OBJECTS) $(PYTHON_INTERFACE).f90
	$(PY) -m numpy.f2py -c -m $(PYTHON_INTERFACE) $(PYTHON_INTERFACE).f90 $(OBJECTS) --f2cmap typemap.dat
# python: $(PYTHON_INTERFACE).pyf $(OBJECTS)
# 	$(PY) -m numpy.f2py -c $(PYTHON_INTERFACE).pyf $(PYTHON_INTERFACE).f90 $(OBJECTS)
#
# $(PYTHON_INTERFACE).pyf: $(PYTHON_INTERFACE).f90
# 	$(PY) -m numpy.f2py $(PYTHON_INTERFACE).f90 -m $(PYTHON_INTERFACE) -h $(PYTHON_INTERFACE).pyf

# build the final executable
$(MAIN).out: $(OBJECTS)
	$(LD) $(OBJECTS) -o$(MAIN).out $(LDFLAGS)
	
%.o: %.f90
	$(FC) -c $(CFLAGS) $< -o $@


clean:
	rm -f *.o
	rm -f *.so
	rm -f *.out
	rm -f *.mod
	rm -f *.pyf
