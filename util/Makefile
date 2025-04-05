#!/bin/bash

# directory where the .mod files are placed
MODDIR ?= ../include
# directory where the libraries are placed
LIBDIR ?= ../lib


ARCHIVE=libutil.a

#OBJECTS=types.o matrix_io.o data_structures.o util.o
OBJECTS=types.o datastructures.o util.o matrix_io.o


# paths to the HDF5 installation
HDF5_LIBDIR ?="/usr/lib/x86_64-linux-gnu/hdf5/serial"
HDF5_INCLUDEDIR ?="/usr/include/hdf5/serial"

# Fortran compiler
FC ?= gfortran
# Fortran compiler flags
CFLAGS += -I$(HDF5_INCLUDEDIR) -I$(LIBDIR) -J$(MODDIR)

# Linker
LD ?= gfortran
# Linker flags
LDFLAGS ?=

# archive utility
AR ?= ar
ARFLAGS ?=

all: $(ARCHIVE)

$(ARCHIVE): $(OBJECTS) $(LIBDIR)
	$(AR) r $@ $(OBJECTS)
	cp $(ARCHIVE) $(LIBDIR)

%.o: %.f90 $(MODDIR)
	$(FC) $(CFLAGS) -c $< -o $@

$(MODDIR):
	mkdir -p $(MODDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

clean:
	rm -f $(OBJECTS) $(ARCHIVE)
	cd $(MODDIR) && rm $(patsubst %.o,m_%.mod, $(OBJECTS))
	rm -f $(LIBDIR)/$(ARCHIVE)

