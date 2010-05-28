# $HeadURL$
# $Id$
FFLAGS= -g -O3 -warn -I.

programs=weno5.x weno3.x weno54.x weno32.x viscouscheck.x
common=doublePrecision.o assorted.o flux.o rhside.o

all: $(programs)

main3.o: FFLAGS += -DWENOORDER=3
main3.o: main.F90
	$(FC) $(FFLAGS) -c -o $@ $<

weno3.x: main3.o reconstruct3.o viscousnop.o $(common)
	$(LD) -o $@ $^

main32.o: FFLAGS += -DWENOORDER=3 -DVISCOUSORDER=2
main32.o: main.F90
	$(FC) $(FFLAGS) -c -o $@ $<

weno32.x: main32.o reconstruct3.o viscous2.o $(common)
	$(LD) -o $@ $^

main5.o: FFLAGS += -DWENOORDER=5
main5.o: main.F90
	$(FC) $(FFLAGS) -c -o $@ $<

weno5.x: main5.o reconstruct5.o viscousnop.o $(common)
	$(LD) -o $@ $^

main54.o: FFLAGS += -DWENOORDER=5 -DVISCOUSORDER=4
main54.o: main.F90
	$(FC) $(FFLAGS) -c -o $@ $<

weno54.x: main54.o reconstruct5.o viscous4.o $(common)
	$(LD) -o $@ $^

viscouscheck.x: viscouscheck.o viscousnop.o viscous2.o viscous4.o assorted.o
	$(LD) -o $@ $^

# Module dependencies
assorted.F90:         doublePrecision.mod
flux.F90:             doublePrecision.mod
main.F90:             doublePrecision.mod
reconstruct3.F90:     doublePrecision.mod
reconstruct5.F90:     doublePrecision.mod
rhside.F90:           doublePrecision.mod
viscous2.F90:         doublePrecision.mod
viscous4.F90:         doublePrecision.mod
viscousnop.F90:       doublePrecision.mod
viscouscheck.F90:     doublePrecision.mod

clean:
	@rm -fv  *.mod *.o *.x *__genmod.f90 *__genmod.mod
	@rm -rfv docs

docs:
	@doxygen

# Use HDF5-enabled toolchain
FC=h5fc
LD=${FC}
RANLIB=touch
AR=ar r

.PHONY: clean docs

.SUFFIXES:
.SUFFIXES: .f .o
.SUFFIXES: .F90 .o
.SUFFIXES: .F90 .mod

.f.o:
	$(FC) $(FFLAGS) -c $<

.F90.o:
	$(FC) $(FFLAGS) -c $<

.F90.mod:
	$(FC) $(FFLAGS) -c $<
