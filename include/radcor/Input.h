// Input of a radcor run, independent of the file format.
//
// RunInput holds everything a run needs (target, angle, spectra, points). It is filled by an
// InputFormat (one class per file format) and checked by validate(); the drivers in Run.h
// take a RunInput and never look at files. Adding a format = writing one InputFormat class
// and registering it in inputFormat().
//
// Formats provided:
//   "ele"  : the unpolarized Fortran file ele-input (fixed line order, list-directed numbers)
//   "pol"  : the polarized Fortran pair radcor_input.dat (header) + ds.inp (data table)
//   "csv"  : a table with named columns, one row per point, run settings in "#!" lines
//   "root" : a TTree "points" with the CSV columns as branches, run settings in a TNamed "run"
//            (only if the library was built with ROOT)
#pragma once
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace radcor {

// The points of one measured quantity on one spectrum (the polarized case has two per spectrum).
struct ConfigData {
  std::vector<double> nu;    // energy transfer E0 - E' [MeV], strictly increasing
  std::vector<double> xs;    // radiated cross section, any unit (<= 0 marks a hole to fill)
  std::vector<double> stat;  // statistical error, same unit (scaled to the Born value in the output)
  std::vector<double> syst;  // systematic error, same unit (unpolarized only; 0 if absent)
};

// One spectrum: one beam energy at the run's angle.
struct SpectrumInput {
  static constexpr double unset = std::numeric_limits<double>::quiet_NaN();
  std::string title;           // a name (only the first 7 characters exist in ele-input)
  double e0 = 0;               // nominal beam energy [MeV]
  double tTarget = 0;          // full target thickness [r.l.]; half of it is placed before the vertex
  double tBefore = 0;          // other material before the vertex [r.l.]
  double tAfter = 0;           // other material after the vertex [r.l.]
  double rcSys = 0;            // fractional uncertainty of the correction (enters SYST only)
  double norm = 1;             // multiplies xs and stat
  double lossS = unset;        // collisional loss of the beam [MeV]; unset: option "deltas"
  double lossP = unset;        // collisional loss of the scattered electron [MeV]; unset: option "deltap"
  std::string accTable;        // finite-acceptance table file; empty: option "acc_cells"
  std::vector<ConfigData> configs;  // 1 (unpolarized) or 2 (polarized)
};

// Everything a run reads.
struct RunInput {
  std::string kind = "unpol";  // "unpol" or "pol"
  int nitt = 10;               // number of unfolding passes
  double Z = 0, A = 0;         // target charge and mass number
  double theta = 0;            // scattering angle [deg]
  double ZW = 0, AW = 0, tTheta = 0;  // wall Z, A and target angle: in the Fortran format, not used
  std::vector<SpectrumInput> spectra;  // in increasing e0
  std::string source;          // where it was read from (for the logs)
  int nconfig() const { return kind == "pol" ? 2 : 1; }
};

// Names of the files that hold one input. Most formats use only `main`.
struct InputFiles {
  std::string main;  // ele-input, radcor_input.dat, the CSV file or the ROOT file
  std::string data;  // "pol" format only: the ds.inp table
};

// One file format: reads a RunInput and writes it back.
class InputFormat {
 public:
  virtual ~InputFormat() = default;
  virtual std::string name() const = 0;
  // Read the files; `kind` ("unpol"/"pol") is what the caller expects.
  virtual RunInput read(const InputFiles& f, const std::string& kind) const = 0;
  // Write `in` to the files (throws if the format cannot hold it).
  virtual void write(const RunInput& in, const InputFiles& f) const = 0;
};

// The format registered under `name` ("ele", "pol", "csv", "root"); throws if unknown.
const InputFormat& inputFormat(const std::string& name);

// Format chosen from the file name: ".csv" -> csv, ".root" -> root, otherwise the Fortran
// format of the run kind (ele for unpol, pol for pol).
std::string guessFormat(const std::string& file, const std::string& kind);

// List of problems that would make the run wrong or impossible (empty = fine):
// target and angle in range, spectra in increasing e0, nu strictly increasing within a spectrum,
// finite numbers, thicknesses >= 0, norm > 0, the polarized configurations on the same nu list.
std::vector<std::string> validate(const RunInput& in);

// validate() and throw std::runtime_error listing every problem.
void checkInput(const RunInput& in);

// Short human-readable summary (one line per spectrum).
std::string describeInput(const RunInput& in);

}  // namespace radcor
