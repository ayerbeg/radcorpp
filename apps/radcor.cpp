// radcor -- command-line front end.
//
//   radcor unpol   [options]   unfold (or radiate) unpolarized spectra, as the unpolarized nradcor.f
//   radcor pol     [options]   the same for two polarized cross-section differences
//   radcor check   unpol|pol [input options]          read and validate an input, print a summary
//   radcor convert unpol|pol [input options] --to F [--to-data F] [--to-format X]
//                                                     write the input in another format
// Input options:
//   --dir D            directory for relative file names (default .)
//   --input F          main input (unpol: ele-input; pol: input/radcor_input.dat)
//   --data F           polarized data table, format pol only (default input/ds.inp)
//   --format X         ele | pol | csv | root (default: from the file name, see guessFormat)
// Run options:
//   --out F            main output (default as the Fortran)
//   --details F        per-point table: tails, soft factor, radiated model, Born, RC, errors
//   --root F           same content as a ROOT TTree (if built with ROOT)
//   --precision P      float (REAL*4 emulation, default) | double
//   --set key=value    any switch of radcor::Options (repeatable), e.g. --set legacy=0
//   -q                 quiet
#include <cstring>
#include <iostream>
#include <string>

#include "radcor/Input.h"
#include "radcor/Options.h"
#include "radcor/Run.h"

// Short help printed on -h or on a usage error.
static void usage() {
  std::cerr << "usage: radcor unpol|pol [--dir D] [--input F] [--data F] [--format ele|pol|csv|root]\n"
               "                       [--out F] [--details F] [--root F] [--precision float|double]\n"
               "                       [--set key=value ...] [-q]\n"
               "       radcor check   unpol|pol [--dir D] [--input F] [--data F] [--format X]\n"
               "       radcor convert unpol|pol [--dir D] [--input F] [--data F] [--format X]\n"
               "                       --to F [--to-data F] [--to-format ele|pol|csv|root]\n"
               "keys: internal external collisional_loss finite_acc radiate delta(=cute) simpson_min(=num)\n"
               "      simpson_step(=del) nitt amass deltas deltap\n"
               "      b del_overwrite spence_legacy ejsz_cap spline_bug spline_y2zero tbta_lag radiate_overwrite\n"
               "      acc_cells legacy (legacy=0 turns off all known legacy bugs)\n";
}

// Parse the command line, read the input in the requested format, then run, check or convert.
int main(int argc, char** argv) {
  if (argc < 2) { usage(); return 1; }
  std::string cmd = argv[1];
  if (cmd == "-h" || cmd == "--help") { usage(); return 0; }
  int a0 = 2;
  std::string kind = cmd;
  if (cmd == "check" || cmd == "convert") {
    if (argc < 3) { usage(); return 1; }
    kind = argv[2];
    a0 = 3;
  }
  if (kind != "unpol" && kind != "pol") { usage(); return 1; }

  radcor::Options opt;
  if (kind == "pol") {  // polarized Fortran defaults (P:123-135)
    opt.amass = 938.272;
    opt.collisional_loss = false;
    opt.finite_acc = false;
    opt.deltas.clear();
    opt.deltap.clear();
  }
  radcor::RunFiles files;
  radcor::Precision prec = radcor::Precision::Float;
  std::string input, data, format, to, toData, toFormat;
  try {
    for (int a = a0; a < argc; ++a) {
      std::string s = argv[a];
      auto next = [&]() -> std::string {
        if (a + 1 >= argc) throw std::runtime_error("missing value after " + s);
        return argv[++a];
      };
      if (s == "--dir") files.dir = next();
      else if (s == "--input") input = next();
      else if (s == "--data") data = next();
      else if (s == "--format") format = next();
      else if (s == "--out") files.output = next();
      else if (s == "--details") files.details = next();
      else if (s == "--root") files.root = next();
      else if (s == "--to") to = next();
      else if (s == "--to-data") toData = next();
      else if (s == "--to-format") toFormat = next();
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

    // --- read the input
    if (input.empty()) input = kind == "unpol" ? "ele-input" : "input/radcor_input.dat";
    if (format.empty()) format = radcor::guessFormat(input, kind);
    if (format == "pol" && data.empty()) data = "input/ds.inp";
    const radcor::InputFiles inFiles{radcor::joinPath(files.dir, input),
                                     data.empty() ? "" : radcor::joinPath(files.dir, data)};
    const radcor::RunInput in = radcor::inputFormat(format).read(inFiles, kind);

    if (cmd == "check") {
      std::cout << radcor::describeInput(in);
      const auto problems = radcor::validate(in);
      for (const auto& p : problems) std::cout << "  problem: " << p << "\n";
      std::cout << (problems.empty() ? "input OK\n" : "input NOT usable\n");
      return problems.empty() ? 0 : 3;
    }
    if (cmd == "convert") {
      if (to.empty()) throw std::runtime_error("convert needs --to FILE");
      if (toFormat.empty()) toFormat = radcor::guessFormat(to, kind);
      if (toFormat == "pol" && toData.empty()) throw std::runtime_error("format pol needs --to-data FILE (ds.inp)");
      radcor::checkInput(in);
      radcor::inputFormat(toFormat).write(in, {radcor::joinPath(files.dir, to),
                                                toData.empty() ? "" : radcor::joinPath(files.dir, toData)});
      if (files.verbose) std::cout << "wrote " << to << (toData.empty() ? "" : " + " + toData) << " (" << toFormat << ")\n";
      return 0;
    }

    if (files.verbose)
      std::cout << "radcor " << kind << " (" << (prec == radcor::Precision::Float ? "float" : "double") << ")\n  "
                << radcor::describeOptions(opt) << "\n";
    return kind == "unpol" ? radcor::runUnpolarized(in, opt, files, prec) : radcor::runPolarized(in, opt, files, prec);
  } catch (const std::exception& e) {
    std::cerr << "radcor: error: " << e.what() << "\n";
    return 2;
  }
}
