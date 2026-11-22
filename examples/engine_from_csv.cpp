// engine_from_csv.cpp -- use the radcor library from your own program: read an input in any
// format into a RunInput, then drive the engine directly (no Fortran-style output stage).
//
// usage: engine_from_csv INPUT [nitt]      INPUT: .csv, .root or an ele-input file
// Prints E0, nu, measured xs, Born xs and RC = measured/Born for every point.
// Settings: double precision, all legacy bugs off, no E94-010 specifics.
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "radcor/Engine.h"
#include "radcor/Input.h"

int main(int argc, char** argv) {
  if (argc < 2) { std::cerr << "usage: engine_from_csv INPUT [nitt]\n"; return 1; }
  try {
    // 1. read and check the input (format from the file name)
    const std::string file = argv[1];
    const radcor::RunInput in = radcor::inputFormat(radcor::guessFormat(file, "unpol")).read({file, ""}, "unpol");
    radcor::checkInput(in);
    const int nitt = argc > 2 ? std::stoi(argv[2]) : in.nitt;

    // 2. settings: defaults, then the recommended ones
    radcor::Options opt;
    radcor::setOption(opt, "legacy=0");
    opt.collisional_loss = false;
    opt.finite_acc = false;

    // 3. the engine: target, angle, then one addSpectrum per beam energy
    radcor::Engine<double> eng(in.Z, in.A, in.theta, opt);
    for (const auto& s : in.spectra) {
      std::vector<double> xs = s.configs[0].xs;
      for (auto& x : xs) x *= s.norm;  // normalisation; no hole filling here
      eng.addSpectrum(s.e0, 0.0, 0.0, s.configs[0].nu, xs, 0.5 * s.tTarget, s.tBefore, s.tAfter);
    }
    eng.run(nitt);

    // 4. results: one PointResult per point (1-based within each spectrum)
    std::printf("#      E0        nu        xs_in      xs_born        RC\n");
    for (size_t i = 0; i < in.spectra.size(); ++i)
      for (int j = 1; j <= eng.spectra()[i].np; ++j) {
        const auto& r = eng.results()[i][j];
        std::printf("%9.2f %9.3f %12.5e %12.5e %9.5f\n", in.spectra[i].e0, r.nu, r.xs_in, r.xs_born, r.xs_in / r.xs_born);
      }
  } catch (const std::exception& e) {
    std::cerr << "engine_from_csv: " << e.what() << "\n";
    return 2;
  }
  return 0;
}
