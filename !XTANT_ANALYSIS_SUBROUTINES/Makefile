# Build the analysis programs listed in sources.txt.
SOURCES := $(shell cat sources.txt)
TARGETS := $(SOURCES:.f90=.x)
FC = ifort
FFLAGS = -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics

.PHONY: default clean

default: $(TARGETS)

%.x: %.f90
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.obj
	rm -f *.x
