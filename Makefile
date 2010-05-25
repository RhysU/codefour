# $HeadURL$
# $Id$
FFLAGS= -g -O0 -warn -I.

programs=weno5.x weno3.x weno54.x weno32.x viscouscheck.x
common=assorted.o doublePrecision.o flux.o rhside.o

all: $(programs)

main3.o: FFLAGS += -DWENOORDER=3
main3.o: main.F90
	$(FC) -o $@ -c $(FFLAGS) $<

weno3.x: main3.o reconstruct3.o viscousnop.o $(common)
	$(LD) -o $@ $^

main32.o: FFLAGS += -DWENOORDER=3 -DVISCOUSORDER=2
main32.o: main.F90
	$(FC) -o $@ -c $(FFLAGS) $<

weno32.x: main32.o reconstruct3.o viscous2.o $(common)
	$(LD) -o $@ $^

main5.o: FFLAGS += -DWENOORDER=5
main5.o: main.F90
	$(FC) -o $@ -c $(FFLAGS) $<

weno5.x: main5.o reconstruct5.o viscousnop.o $(common)
	$(LD) -o $@ $^

main54.o: FFLAGS += -DWENOORDER=5 -DVISCOUSORDER=4
main54.o: main.F90
	$(FC) -o $@ -c $(FFLAGS) $<

weno54.x: main54.o reconstruct5.o viscous4.o $(common)
	$(LD) -o $@ $^

viscouscheck.x: viscouscheck.o viscousnop.o viscous2.o viscous4.o assorted.o
	$(LD) -o $@ $^

# Module dependencies
assorted.F90:     doublePrecision.mod
flux.F90:         doublePrecision.mod
main.F90:         doublePrecision.mod
reconstruct3.F90: doublePrecision.mod
reconstruct5.F90: doublePrecision.mod
rhside.F90:       doublePrecision.mod
viscous2.F90:     doublePrecision.mod
viscous4.F90:     doublePrecision.mod
viscousnop.F90:   doublePrecision.mod

clean:
	@rm -fv *.mod *.o *.x *__genmod.f90 *__genmod.mod

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
.SUFFIXES: .F90 .o
.SUFFIXES: .F90 .mod

.C.o:
	$(CXX) -c $(CFLAGS) $<

.c.o:
	$(CC) -c $(cFLAGS) $<

.f.o:
	$(FC) -c $(FFLAGS) $<

.F90.o:
	$(FC) -c $(FFLAGS) $<

.F90.mod:
	$(FC) -c $(FFLAGS) $<
