.SUFFIXES:
.SUFFIXES: .C .o 
.SUFFIXES: .c .o
.SUFFIXES: .f .o
.SUFFIXES: .f90 .o
.SUFFIXES: .f90 .mod

CC= g++
cc = gcc
fc = gfortran

CFLAGS= -O3 -I. 
#FFLAGS= -g -W -pedantic -ggdb -gstabs+ -g3
FFLAGS= -O 

LD=gfortran
RANLIB=touch
AR= ar r

.C.o:
	$(CC) -c $(CFLAGS) $<

.c.o:	
	$(cc) -c $(cFLAGS) $<

.f.o:	
	$(fc) -c $(FFLAGS) $<

.f90.o:	
	$(fc) -c $(FFLAGS) $<

.f90.mod:	
	$(fc) -c $(FFLAGS) $<

#
# object files 
#
OBJS= main.o doublePrecision.o reconstruct.o assorted.o flux.o rhside.f90

# the dependencies between modules 
main.f90: doublePrecision.mod
reconstruct.f90: doublePrecision.mod
assorted.f90: doublePrecision.mod
flux.f90: doublePrecision.mod
rhside.f90: doublePrecision.mod

weno5.x: $(OBJS)
	$(LD) -o $@ $(OBJS) 

clean:
	rm -f *.o *.a *.mod 
