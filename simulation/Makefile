#!/bin/bash

# directory where the .mod files are placed
MODDIR ?=../include
# directory where the libraries are placed
LIBDIR ?=../lib


ARCHIVE=libsimulation.a

OBJECTS=simulation.o simulation_io.o

# Fortran compiler
FC ?=gfortran
# Fortran compiler flags
CFLAGS += -I$(LIBDIR) -J$(MODDIR) -fopenmp

# Linker
LD ?=gfortran
# Linker flags
LDFLAGS ?= -fopenmp

# archive utility
AR ?= ar
ARFLAGS ?=

all: $(ARCHIVE)

$(ARCHIVE): $(OBJECTS) $(LIBDIR)
	$(AR) r $@ $(OBJECTS)
	cp $(ARCHIVE) $(LIBDIR)

%.o: %.f90 $(MODDIR)
	$(FC) -c $(CFLAGS) $< -o $@

$(MODDIR):
	mkdir -p $(MODDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

clean:
	rm -f $(OBJECTS) $(ARCHIVE)
	cd $(MODDIR) && rm $(patsubst %.o,m_%.mod, $(OBJECTS))
	rm -f $(LIBDIR)/$(ARCHIVE)

