// engine_from_csv.cpp -- drive the radcor engine directly from a CSV table, without the
// Fortran-style input files. Shows the minimum a program must provide to the engine.
//
// CSV columns (a header line first; the columns must be in exactly this order):
//   E0,nu,xs,stat,syst,t_target,t_before,t_after
// one row per point; rows of one spectrum are consecutive, nu increasing within a spectrum,
// spectra in increasing E0. Energies in MeV, thicknesses in radiation lengths, xs in any unit.
//
// usage: engine_from_csv points.csv Z A theta_deg [nitt]
// Prints E0, nu, measured xs, Born xs and RC = measured/Born for every point.
// Settings: double precision, all legacy bugs off, no E94-010 specifics.
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "radcor/Engine.h"

struct Block {                // one spectrum as read from the CSV
  double e0 = 0, tTarget = 0, tBefore = 0, tAfter = 0;
  std::vector<double> nu, xs;
};

// Read the CSV into spectra, checking the ordering the engine relies on.
static std::vector<Block> readCsv(const std::string& file) {
  std::ifstream in(file);
  if (!in) throw std::runtime_error("cannot open " + file);
  std::string line;
  std::getline(in, line);  // header
  std::vector<Block> out;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    for (auto& c : line) if (c == ',') c = ' ';
    std::istringstream ss(line);
    double e0, nu, xs, stat, syst, tt, tb, ta;
    if (!(ss >> e0 >> nu >> xs >> stat >> syst >> tt >> tb >> ta)) throw std::runtime_error("bad row: " + line);
    if (out.empty() || e0 != out.back().e0) {
      if (!out.empty() && e0 < out.back().e0) throw std::runtime_error("spectra must be in increasing E0");
      out.push_back({e0, tt, tb, ta, {}, {}});
    }
    Block& b = out.back();
    if (!b.nu.empty() && nu <= b.nu.back()) throw std::runtime_error("nu must increase within a spectrum");
    b.nu.push_back(nu);
    b.xs.push_back(xs);
  }
  return out;
}

int main(int argc, char** argv) {
  if (argc < 5) { std::cerr << "usage: engine_from_csv points.csv Z A theta_deg [nitt]\n"; return 1; }
  try {
    const auto blocks = readCsv(argv[1]);
    const double Z = std::stod(argv[2]), A = std::stod(argv[3]), theta = std::stod(argv[4]);
    const int nitt = argc > 5 ? std::stoi(argv[5]) : 10;

    radcor::Options opt;                       // defaults, then the recommended settings
    radcor::setOption(opt, "legacy=0");
    opt.collisional_loss = false;
    opt.finite_acc = false;

    radcor::Engine<double> eng(Z, A, theta, opt);
    for (const auto& b : blocks)               // no collisional loss: lossS = lossP = 0
      eng.addSpectrum(b.e0, 0.0, 0.0, b.nu, b.xs, 0.5 * b.tTarget, b.tBefore, b.tAfter);
    eng.run(nitt);

    std::printf("#      E0        nu        xs_in      xs_born        RC\n");
    for (size_t i = 0; i < blocks.size(); ++i)
      for (int j = 1; j <= eng.spectra()[i].np; ++j) {
        const auto& r = eng.results()[i][j];
        std::printf("%9.2f %9.3f %12.5e %12.5e %9.5f\n", blocks[i].e0, r.nu, r.xs_in, r.xs_born, r.xs_in / r.xs_born);
      }
  } catch (const std::exception& e) {
    std::cerr << "engine_from_csv: " << e.what() << "\n";
    return 2;
  }
  return 0;
}
