// radcor — command-line front end.
//
//   radcor unpol [options]   reproduce unpolarized_RC/nradcor.f
//   radcor pol   [options]   reproduce polarized_RC/nradcor.f
//
// Options:
//   --dir D            working directory containing the Fortran-style inputs (default .)
//   --input F          main input (unpol: ele-input; pol: input/radcor_input.dat)
//   --data F           polarized data (default input/ds.inp)
//   --out F            main output (default as the Fortran)
//   --details F        per-point table: tails, soft factor, radiated model, Born, RC, errors
//   --root F           same content as a ROOT TTree (if built with ROOT)
//   --precision P      float (REAL*4 emulation, default) | double
//   --set key=value    any switch of radcor::Options (repeatable), e.g. --set legacy=0
//   -q                 quiet
#include <cstring>
#include <iostream>
#include <string>

#include "radcor/Legacy.h"
#include "radcor/Options.h"

// Short help printed on -h or on a usage error.
static void usage() {
  std::cerr << "usage: radcor unpol|pol [--dir D] [--input F] [--data F] [--out F] [--details F] [--root F]\n"
               "                       [--precision float|double] [--set key=value ...] [-q]\n"
               "keys: internal external collisional_loss finite_acc radiate delta(=cute) simpson_min(=num)\n"
               "      simpson_step(=del) nitt amass deltas deltap\n"
               "      b del_overwrite spence_legacy ejsz_cap spline_bug spline_y2zero tbta_lag radiate_overwrite\n"
               "      acc_cells legacy (legacy=0 turns off all known legacy bugs)\n";
}

// Parse the command line into Options and RunFiles, then run the selected program.
int main(int argc, char** argv) {
  if (argc < 2) { usage(); return 1; }
  const std::string mode = argv[1];
  if (mode != "unpol" && mode != "pol") { usage(); return 1; }
  radcor::Options opt;
  if (mode == "pol") {  // polarized Fortran defaults (P:123-135)
    opt.amass = 938.272;
    opt.collisional_loss = false;
    opt.finite_acc = false;
    opt.deltas.clear();
    opt.deltap.clear();
  }
  radcor::RunFiles files;
  radcor::Precision prec = radcor::Precision::Float;
  try {
    for (int a = 2; a < argc; ++a) {
      std::string s = argv[a];
      auto next = [&]() -> std::string {
        if (a + 1 >= argc) throw std::runtime_error("missing value after " + s);
        return argv[++a];
      };
      if (s == "--dir") files.dir = next();
      else if (s == "--input") files.input = next();
      else if (s == "--data") files.data = next();
      else if (s == "--out") files.output = next();
      else if (s == "--details") files.details = next();
      else if (s == "--root") files.root = next();
      else if (s == "--precision") {
        std::string p = next();
        if (p == "float") prec = radcor::Precision::Float;
        else if (p == "double") prec = radcor::Precision::Double;
        else throw std::runtime_error("precision must be float or double");
      } else if (s == "--set") radcor::setOption(opt, next());
      else if (s == "-q") files.verbose = false;
      else if (s == "-h" || s == "--help") { usage(); return 0; }
      else throw std::runtime_error("unknown argument " + s);
    }
    if (files.verbose)
      std::cout << "radcor " << mode << " (" << (prec == radcor::Precision::Float ? "float" : "double") << ")\n  "
                << radcor::describeOptions(opt) << "\n";
    return mode == "unpol" ? radcor::runUnpolarized(opt, files, prec) : radcor::runPolarized(opt, files, prec);
  } catch (const std::exception& e) {
    std::cerr << "radcor: error: " << e.what() << "\n";
    return 2;
  }
}
