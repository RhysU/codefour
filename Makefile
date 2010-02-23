CFLAGS= -O3 -I.
#FFLAGS= -g -W -pedantic -ggdb -gstabs+ -g3
FFLAGS= -O

objects=assorted.o doublePrecision.o flux.o main.o reconstruct.o rhside.o

weno5.x: $(objects)
	$(LD) -o $@ $(OBJS)

# Module dependencies
assorted.f90:    doublePrecision.mod
flux.f90:        doublePrecision.mod
main.f90:        doublePrecision.mod
reconstruct.f90: doublePrecision.mod
rhside.f90:      doublePrecision.mod

clean:
	@rm -fv *.o *.a *.mod

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
