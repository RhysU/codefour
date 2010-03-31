# $HeadURL$
# $Id$
CFLAGS= -O3 -I.
#FFLAGS= -g -W -pedantic -ggdb -gstabs+ -g3
FFLAGS= -O

programs=weno5.x weno3.x
objects=assorted.o doublePrecision.o flux.o rhside.o \
		reconstruct3.o reconstruct5.o

all: $(programs)

main3.o: FFLAGS += -cpp -DRECONSTRUCT_FUNCTION=reconstruct3
main3.o: main.f90
	$(FC) -o $@ -c $(FFLAGS) $<

weno3.x: main3.o $(objects)
	$(LD) -o $@ $^

main5.o: FFLAGS += -cpp -DRECONSTRUCT_FUNCTION=reconstruct5
main5.o: main.f90
	$(FC) -o $@ -c $(FFLAGS) $<

weno5.x: main5.o $(objects)
	$(LD) -o $@ $^

# Module dependencies
assorted.f90:     doublePrecision.mod
flux.f90:         doublePrecision.mod
main.f90:         doublePrecision.mod
reconstruct3.f90: doublePrecision.mod
reconstruct5.f90: doublePrecision.mod
rhside.f90:       doublePrecision.mod

clean:
	@rm -fv *.mod *.o *.x

# Use GNU compilers if choice not present in environment
ifndef CXX
CXX=g++
endif
ifndef CC
CC=gcc
endif
ifndef FC
FC=gfortran
endif
ifeq ($(FC),f77)
FC=gfortran
endif

LD=${FC}
RANLIB=touch
AR=ar r

.SUFFIXES:
.SUFFIXES: .C .o
.SUFFIXES: .c .o
.SUFFIXES: .f .o
.SUFFIXES: .f90 .o
.SUFFIXES: .f90 .mod

.C.o:
	$(CXX) -c $(CFLAGS) $<

.c.o:
	$(CC) -c $(cFLAGS) $<

.f.o:
	$(FC) -c $(FFLAGS) $<

.f90.o:
	$(FC) -c $(FFLAGS) $<

.f90.mod:
	$(FC) -c $(FFLAGS) $<
