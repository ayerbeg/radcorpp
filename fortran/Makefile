# Build the two Fortran programs translated by radcor (gfortran >= 10).
# The flags are the minimum needed by the 2005 source; -O0 keeps the arithmetic of the reference outputs.
FC     ?= gfortran
FFLAGS ?= -std=legacy -ffixed-line-length-none -fdec-char-conversions -O0 -g

all: unpolarized/nradcor polarized/nradcor

unpolarized/nradcor: unpolarized/nradcor.f
	$(FC) $(FFLAGS) -o $@ $<
polarized/nradcor: polarized/nradcor.f
	$(FC) $(FFLAGS) -o $@ $<

clean:
	rm -f unpolarized/nradcor polarized/nradcor
.PHONY: all clean
