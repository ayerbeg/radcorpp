// Drivers for a complete run: take a RunInput (from any format, see Input.h), fill the
// Engine, unfold, and write the outputs:
//   runUnpolarized : as the unpolarized nradcor.f (main output in its format, E94010_xs_unfolded.dat)
//   runPolarized   : as the polarized nradcor.f   (main output ds.rc{1,2,3})
// plus the optional per-point "details" table (text) and ROOT tree: the two tail integrals,
// the soft factor, the radiated model, the Born result, RC, errors and acceptance factor.
#pragma once
#include <string>

#include "radcor/Input.h"
#include "radcor/Options.h"

namespace radcor {

// Float reproduces the REAL*4 Fortran (parity tests); Double is the production mode.
enum class Precision { Float, Double };

// Output file names of a run. Relative names are taken relative to dir, which is also where
// the E94-010 acceptance tables input/<cell>_qfs_finite_acceptance.out are looked for.
struct RunFiles {
  std::string dir = ".";       // working directory (like the Fortran cwd)
  std::string output;          // main output ("" = the Fortran default name)
  std::string details;         // optional details table ("" = none)
  std::string root;            // optional ROOT file ("" = none)
  bool verbose = true;         // print "Es iteration" progress like the Fortran
};

// Run the unpolarized / polarized calculation on `in` (checked with checkInput first).
// Return 0 on success; throw std::runtime_error on invalid input or unwritable files.
int runUnpolarized(const RunInput& in, const Options& opt, const RunFiles& files, Precision prec);
int runPolarized(const RunInput& in, const Options& opt, const RunFiles& files, Precision prec);

// Path of file f relative to dir (absolute paths and dir "." are left as they are).
std::string joinPath(const std::string& dir, const std::string& f);

// Fortran edit descriptors Ew.d and Fw.d (gfortran style).
std::string fortranE(double v, int w, int d);
std::string fortranF(double v, int w, int d);

}  // namespace radcor
