// Drivers that reproduce the two Fortran programs end to end, with their input
// formats and output files:
//   runUnpolarized : unpolarized_RC/nradcor.f  (ele-input -> output/E94010_xs_unfolded.dat)
//   runPolarized   : polarized_RC/nradcor.f    (input/radcor_input.dat + input/ds.inp -> output/ds.rc{1,2,3})
// Extra outputs (not in the Fortran): a per-point "details" table with the two tail integrals,
// the soft factor, the radiated model, the Born result, RC and errors; optionally a ROOT file
// with the same content.
#pragma once
#include <string>

#include "radcor/Options.h"

namespace radcor {

// Float reproduces the REAL*4 Fortran (parity tests); Double is the production mode.
enum class Precision { Float, Double };

// File names of a run. Relative names are taken relative to dir.
struct RunFiles {
  std::string dir = ".";       // working directory (like the Fortran cwd)
  std::string input;           // main input file (default: ele-input / input/radcor_input.dat)
  std::string data;            // polarized data file (default: input/ds.inp)
  std::string output;          // main output (default as in Fortran)
  std::string details;         // optional details table ("" = none)
  std::string root;            // optional ROOT file ("" = none)
  bool verbose = true;         // print "Es iteration" progress like the Fortran
};

// Run the unpolarized / polarized program; return 0 on success, throw std::runtime_error on bad input.
int runUnpolarized(const Options& opt, const RunFiles& files, Precision prec);
int runPolarized(const Options& opt, const RunFiles& files, Precision prec);

// Fortran edit descriptors Ew.d and Fw.d (gfortran style).
std::string fortranE(double v, int w, int d);
std::string fortranF(double v, int w, int d);

}  // namespace radcor
