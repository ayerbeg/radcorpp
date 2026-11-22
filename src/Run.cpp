// Drivers: take a RunInput, fill the Engine, unfold, and write the outputs in the Fortran
// formats, plus the optional details/ROOT tables. No file parsing here (see Input.cpp).
#include "radcor/Run.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "Details.h"
#include "radcor/Engine.h"
#include "radcor/Interp.h"

namespace radcor {

// ---------------------------------------------------------------------------
// Fortran-style formatting
// Ew.d: mantissa 0.ddd, exponent E+xx, right-aligned in w characters ('*' fill on overflow).
std::string fortranE(double v, int w, int d) {
  std::string s;
  if (std::isnan(v)) s = "NaN";
  else if (std::isinf(v)) s = v > 0 ? "Infinity" : "-Infinity";
  else if (v == 0.0) s = "0." + std::string(d, '0') + "E+00";
  else {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*e", d - 1, std::fabs(v));  // d.ddde+xx
    std::string b(buf);
    auto epos = b.find('e');
    std::string mant = b.substr(0, epos);
    int ex = std::stoi(b.substr(epos + 1)) + 1;
    std::string digits;
    for (char c : mant) if (c != '.') digits += c;
    char eb[16];
    if (std::abs(ex) <= 99) std::snprintf(eb, sizeof eb, "E%c%02d", ex < 0 ? '-' : '+', std::abs(ex));
    else std::snprintf(eb, sizeof eb, "%c%03d", ex < 0 ? '-' : '+', std::abs(ex));
    s = std::string(v < 0 ? "-" : "") + "0." + digits + eb;
  }
  if (static_cast<int>(s.size()) > w) return std::string(w, '*');
  return std::string(w - s.size(), ' ') + s;
}

// Fw.d: fixed point, d decimals, right-aligned in w characters ('*' fill on overflow).
std::string fortranF(double v, int w, int d) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%*.*f", w, d, v);
  std::string s(buf);
  if (static_cast<int>(s.size()) > w) return std::string(w, '*');
  return s;
}

// ---------------------------------------------------------------------------
std::string joinPath(const std::string& dir, const std::string& f) {
  if (f.empty() || f[0] == '/' || dir.empty() || dir == ".") return f;
  return dir + "/" + f;
}

// Invariant mass W of the hadronic final state, proton mass (GET_W, U:1920).
template <class T>
static T invariantMassW(T E0, T TH, T XNU) {
  const T XM = T(938.272);
  const T PI = std::acos(T(-1));
  const T THR = TH * PI / T(180);
  const T S2 = fm::sq(std::sin(THR / T(2)));
  const T EP = E0 - XNU;
  const T Q2 = T(4) * E0 * EP * S2;
  const T W2 = fm::sq(XM) - Q2 + T(2) * XM * XNU;
  return std::sqrt(W2);
}

// Header of the input file (shared by the unpolarized and polarized formats): NITT, NE,
// target, wall, and one thickness/normalisation line per spectrum.
// Collisional loss of spectrum i: the value given in the input if any, else the option list
// (v[i], or 0 if the list is shorter).
static double lossOf(double given, const std::vector<double>& v, int i) {
  if (!std::isnan(given)) return given;
  return i < static_cast<int>(v.size()) ? v[i] : 0.0;
}

// Complete description of a run, written into the details/ROOT outputs.
template <class T>
static std::string runSettings(const char* mode, const Options& opt, const RunInput& in, int nitt) {
  std::ostringstream s;
  s << "radcor " << mode << " (" << (sizeof(T) == 4 ? "float" : "double") << ")\n"
    << "  input=" << in.source << "\n"
    << "  Z=" << in.Z << " A=" << in.A << " theta=" << in.theta << " NE=" << in.spectra.size() << " NITT=" << nitt << "\n"
    << "  " << describeOptions(opt);
  return s.str();
}

// Copy the engine's per-point result into an output row (energies, terms, Born, RC).
template <class T>
static DetailRow detailOf(const PointResult<T>& p) {
  DetailRow d;
  d.es = double(p.es); d.ep = double(p.ep); d.nu = double(p.nu); d.xs_in = double(p.xs_in);
  d.tail_es = double(p.tailEs); d.tail_ep = double(p.tailEp); d.soft_factor = double(p.softFactor);
  d.rad_model = double(p.radModel); d.xs_born = double(p.xs_born); d.xs_born_prev = double(p.xs_born_prev);
  d.rc = d.xs_born != 0 ? d.xs_in / d.xs_born : 0.0;
  return d;
}

// ===========================================================================
// Unpolarized
// Whole unpolarized run in precision T: apply the normalisation and hole filling, unfold, then
// write the main output (Born on the nominal nu grid) or, in RADIATE mode, print the radiated
// model; finally the details/ROOT tables.
template <class T>
static int runUnpolT(const RunInput& in, const Options& opt, const RunFiles& files) {
  checkInput(in);
  if (in.kind != "unpol") throw std::runtime_error("unpolarized run needs an unpolarized input");
  const int ne = static_cast<int>(in.spectra.size());
  const int nitt = opt.nitt_override > 0 ? opt.nitt_override : in.nitt;
  const std::string settings = runSettings<T>("unpol", opt, in, nitt);

  Engine<T> eng(T(in.Z), T(in.A), T(in.theta), opt);
  if (files.verbose) eng.progress = [](T e, int n) { std::cout << "   " << double(e) << "   " << n << "\n"; };

  struct Raw { std::vector<T> nu, xs, stat, syst; T e0; };
  std::vector<Raw> raw(ne);
  for (int i = 0; i < ne; ++i) {
    const SpectrumInput& S = in.spectra[i];
    const ConfigData& c = S.configs[0];
    Raw& R = raw[i];
    R.e0 = T(S.e0);
    const int np = static_cast<int>(c.nu.size());
    for (int j = 0; j < np; ++j) {
      R.nu.push_back(T(c.nu[j])); R.xs.push_back(T(c.xs[j])); R.stat.push_back(T(c.stat[j]));
      R.syst.push_back(c.syst.empty() ? T(0) : T(c.syst[j]));
    }
    // RNORME and hole filling, in the Fortran order (U:345-361): a point with xs <= 0 takes the
    // mean of its neighbours (the right one not yet normalised, as in the Fortran).
    const T rn = T(S.norm);
    for (int j = 1; j <= np; ++j) {
      T& x = R.xs[j - 1];
      x = x * rn;
      R.stat[j - 1] = R.stat[j - 1] * rn;
      if (j >= 2 && j <= np - 1 && x <= T(0)) {
        if (R.xs[j - 2] != T(0) && R.xs[j] != T(0)) {
          x = (R.xs[j - 2] + R.xs[j]) / T(2);
          R.stat[j - 1] = R.stat[j - 2] + R.stat[j];
        }
      }
    }
    const T lossS = opt.collisional_loss ? T(lossOf(S.lossS, opt.deltas, i)) : T(0);
    const T lossP = opt.collisional_loss ? T(lossOf(S.lossP, opt.deltap, i)) : T(0);
    eng.addSpectrum(R.e0, lossS, lossP, R.nu, R.xs, T(0.5) * T(S.tTarget), T(S.tBefore), T(S.tAfter));
  }

  // finite acceptance tables (U:192-194, GET_FINITE_COR): per-spectrum file if given, else the
  // E94-010 cell names of option acc_cells
  std::vector<std::vector<T>> accNu(ne), accF(ne);
  if (opt.finite_acc) {
    for (int i = 0; i < ne; ++i) {
      std::string name = in.spectra[i].accTable;
      if (name.empty()) {
        if (i >= static_cast<int>(opt.acc_cells.size()))
          throw std::runtime_error("finite_acc: no acceptance table for spectrum " + std::to_string(i + 1));
        name = "input/" + opt.acc_cells[i] + "_qfs_finite_acceptance.out";
      }
      const std::string f = joinPath(files.dir, name);
      std::ifstream is(f);
      if (!is) throw std::runtime_error("cannot open " + f);
      std::string line;
      accNu[i].push_back(T(0)); accF[i].push_back(T(0));  // 1-based
      while (std::getline(is, line)) {
        if (!line.empty() && line[0] == '#') continue;
        std::istringstream ss(line);
        double a, b, c, d;
        if (ss >> a >> b >> c >> d) { accNu[i].push_back(T(a)); accF[i].push_back(T(d)); }
      }
    }
  }

  // Output header (U:171-190); written before the run as in the Fortran.
  std::string outName = joinPath(files.dir, files.output.empty() ? "output/E94010_xs_unfolded.dat" : files.output);
  std::ofstream out(outName);
  if (!out) throw std::runtime_error("cannot open output " + outName);
  out << "# Energy in MeV.  cross sections in pb/MeV-sr.    \n# XSi : before RC\n# XSf : after  RC\n";
  auto ra = [](const std::string& s, int w) { return std::string(w > (int)s.size() ? w - s.size() : 0, ' ') + s; };
  if (opt.finite_acc) {
    out << "# XSfc: after  RC, with finite acceptance correction\n#\n";
    out << "#" << ra("E0", 5) << " " << ra("E-E'", 6) << ra("XSi", 12) << ra("XSf", 12) << ra("XSfc", 12)
        << ra("STAT", 12) << ra("SYST", 12) << "\n";
  } else {
    out << "#\n#" << ra("E0", 5) << " " << ra("E-E'", 6) << ra("XSi", 12) << ra("XSf", 12) << ra("STAT", 12)
        << ra("SYST", 12) << "\n";
  }
  out << "#" << std::string(72, '-') << "\n";

  eng.run(nitt);
  const auto& res = eng.results();
  const auto& sp = eng.spectra();

  std::vector<DetailRow> rows;
  for (int i = 0; i < ne; ++i)
    for (int j = 1; j <= sp[i].np; ++j) {
      const auto& p = res[i][j];
      DetailRow d = detailOf<T>(p);
      d.spectrum = i + 1; d.point = j; d.e0 = double(sp[i].e0);
      // errors of the Born value, same formula as the main output but at these kinematics
      const double scale = d.xs_in != 0 ? d.xs_born / d.xs_in : 0.0;
      d.stat = double(raw[i].stat[j - 1]) * scale;
      d.syst = std::sqrt(std::pow(double(raw[i].syst[j - 1]) * scale, 2) +
                         std::pow(std::abs(d.xs_in - d.xs_born) * in.spectra[i].rcSys, 2));
      if (opt.finite_acc) {
        T f = 0;
        interpPoly(static_cast<int>(accNu[i].size()) - 1, accNu[i], accF[i], T(d.nu), f, 15, 0);
        d.acc_factor = double(f);
      }
      rows.push_back(d);
    }

  if (opt.radiate) {  // U:584-588 then STOP (U:634)
    const int i = ne - 1;
    for (int j = 1; j <= sp[i].np; ++j) {
      const auto& p = res[i][j];
      const T nuR = sp[i].es - sp[i].ep[j];
      const T wR = invariantMassW(sp[i].es, T(in.theta), nuR);
      std::cout << fortranF(double(sp[i].es), 7, 1) << fortranF(double(nuR), 7, 1) << fortranF(double(wR), 7, 1)
                << fortranE(double(p.xs_in), 12, 4) << fortranE(double(p.radModel), 12, 4) << "\n";
    }
  } else {
    for (int i = 0; i < ne; ++i) {
      const int np = sp[i].np;
      // Born and its errors on the grid actually used (nu = E_s - E_p after the shift) ...
      std::vector<T> nuUsed(np + 1), born(np + 1), statBorn(np + 1), systBorn(np + 1), meas(np + 1);
      for (int j = 1; j <= np; ++j) {
        const auto& p = res[i][j];
        meas[j] = p.xs_in;
        born[j] = p.xs_born;
        const T scale = born[j] / meas[j];
        statBorn[j] = raw[i].stat[j - 1] * scale;
        const T systRC = std::abs(meas[j] - born[j]) * T(in.spectra[i].rcSys);
        T s2 = fm::sq(raw[i].syst[j - 1] * scale) + fm::sq(systRC);
        systBorn[j] = std::sqrt(s2);
        nuUsed[j] = sp[i].es - sp[i].ep[j];
      }
      // ... then read back by cubic spline at the nominal nu of each input row (SPLINE/SPLINT)
      const T y2z = T(opt.spline_y2zero);
      auto curvBorn = splineSecondDerivs(nuUsed, born, np, T(0), T(0), opt.spline_bug, y2z);
      auto curvStat = splineSecondDerivs(nuUsed, statBorn, np, T(0), T(0), opt.spline_bug, y2z);
      auto curvSyst = splineSecondDerivs(nuUsed, systBorn, np, T(0), T(0), opt.spline_bug, y2z);
      const T shift = opt.collisional_loss ? T(lossOf(in.spectra[i].lossS, opt.deltas, i)) + T(lossOf(in.spectra[i].lossP, opt.deltap, i)) : T(0);
      T accFactor = 0;
      for (int j = 1; j <= np - 1; ++j) {  // last point not written (U:690)
        const T nuOut = opt.collisional_loss ? nuUsed[j] + shift : nuUsed[j];
        const T bornOut = splineEval(nuUsed, born, curvBorn, np, nuOut);
        const T statOut = splineEval(nuUsed, statBorn, curvStat, np, nuOut);
        const T systOut = splineEval(nuUsed, systBorn, curvSyst, np, nuOut);
        out << fortranF(double(sp[i].e0), 6, 1) << " " << fortranF(double(nuOut), 6, 1)
            << fortranE(double(meas[j]), 12, 4) << fortranE(double(bornOut), 12, 4);
        if (opt.finite_acc) {
          interpPoly(static_cast<int>(accNu[i].size()) - 1, accNu[i], accF[i], nuOut, accFactor, 15, 0);
          out << fortranE(double(bornOut * accFactor), 12, 4);
        }
        out << fortranE(double(statOut), 12, 4) << fortranE(double(systOut), 12, 4) << "\n";
      }
    }
  }
  if (!files.details.empty()) writeDetailsText(joinPath(files.dir, files.details), rows, settings);
  if (!files.root.empty() && !writeDetailsRoot(joinPath(files.dir, files.root), rows, settings))
    std::cerr << "warning: built without ROOT, no ROOT file written\n";
  return 0;
}

// ===========================================================================
// Polarized
// Whole polarized run in precision T: unfold each configuration (cross-section difference)
// with its own Engine, write ds.rc{1,2,3}.
template <class T>
static int runPolT(const RunInput& in, const Options& opt, const RunFiles& files) {
  checkInput(in);
  if (in.kind != "pol") throw std::runtime_error("polarized run needs a polarized input");
  const int ne = static_cast<int>(in.spectra.size());
  const int nitt = opt.nitt_override > 0 ? opt.nitt_override : in.nitt;
  const std::string settings = runSettings<T>("pol", opt, in, nitt);

  struct Raw { T e0; std::vector<T> nu[2], xs[2], err[2]; };
  std::vector<Raw> raw(ne);
  std::vector<std::unique_ptr<Engine<T>>> eng;
  for (int k = 0; k < 2; ++k) {
    eng.emplace_back(new Engine<T>(T(in.Z), T(in.A), T(in.theta), opt));
    for (int i = 0; i < ne; ++i) {
      const SpectrumInput& S = in.spectra[i];
      const ConfigData& c = S.configs[k];
      Raw& R = raw[i];
      R.e0 = T(S.e0);
      for (size_t j = 0; j < c.nu.size(); ++j) {
        R.nu[k].push_back(T(c.nu[j])); R.xs[k].push_back(T(c.xs[j])); R.err[k].push_back(T(c.stat[j]));
      }
      std::vector<T> xs = R.xs[k];
      const T rn = T(S.norm);
      for (auto& x : xs) x = x * rn;
      for (auto& e : R.err[k]) e = e * rn;
      const T lossS = opt.collisional_loss ? T(lossOf(S.lossS, opt.deltas, i)) : T(0);
      const T lossP = opt.collisional_loss ? T(lossOf(S.lossP, opt.deltap, i)) : T(0);
      eng[k]->addSpectrum(R.e0, lossS, lossP, R.nu[k], xs, T(0.5) * T(S.tTarget), T(S.tBefore), T(S.tAfter));
    }
    if (files.verbose) std::cout << "  configuration " << k + 1 << "\n";
    eng[k]->run(nitt);
  }

  std::string outDefault = opt.internal && opt.external ? "output/ds.rc3" : (opt.internal ? "output/ds.rc2" : "output/ds.rc1");
  std::string outName = joinPath(files.dir, files.output.empty() ? outDefault : files.output);
  std::ofstream out(outName);
  if (!out) throw std::runtime_error("cannot open output " + outName);

  std::vector<DetailRow> rows;
  for (int i = 0; i < ne; ++i) {
    const auto& s1 = eng[0]->spectra()[i];
    const auto& s2 = eng[1]->spectra()[i];
    const int np = s1.np;
    std::vector<T> nuUsed(np + 1), born1(np + 1), born2(np + 1);
    for (int j = 1; j <= np; ++j) {
      nuUsed[j] = s1.es - s1.ep[j];
      born1[j] = s1.mott * s1.bornRatio[j];
      born2[j] = s2.mott * s2.bornRatio[j];
    }
    const T shift = opt.collisional_loss ? T(lossOf(in.spectra[i].lossS, opt.deltas, i)) + T(lossOf(in.spectra[i].lossP, opt.deltap, i)) : T(0);
    T DS1F = 0, DS2F = 0;  // Born of configurations 1 and 2 at the output nu
    for (int j = 1; j <= np; ++j) {
      const T DS1I = s1.mott * s1.measRatio[j];
      const T DS2I = s2.mott * s2.measRatio[j];
      const T nuOut = opt.collisional_loss ? nuUsed[j] + shift : nuUsed[j];
      interpPoly(np, nuUsed, born1, nuOut, DS1F, 1, 1);
      interpPoly(np, nuUsed, born2, nuOut, DS2F, 1, 1);
      if (j == 1) { if (DS1I == T(0)) DS1F = 0; if (DS2I == T(0)) DS2F = 0; }
      const T wOut = invariantMassW(s1.es, T(in.theta), nuOut);
      out << fortranF(double(raw[i].e0), 10, 1) << fortranF(double(nuOut), 10, 1) << fortranF(double(wOut), 10, 1)
          << fortranE(double(DS1F), 15, 5) << fortranE(double(raw[i].err[0][j - 1]), 15, 5)
          << fortranE(double(DS2F), 15, 5) << fortranE(double(raw[i].err[1][j - 1]), 15, 5) << "\n";
    }
    for (int k = 0; k < 2; ++k) {
      const auto& res = eng[k]->results()[i];
      for (int j = 1; j <= np; ++j) {
        DetailRow d = detailOf<T>(res[j]);
        d.spectrum = i + 1; d.config = k + 1; d.point = j; d.e0 = double(raw[i].e0);
        d.stat = double(raw[i].err[k][j - 1]);  // copied unchanged, as the Fortran output does
        d.syst = 0;
        rows.push_back(d);
      }
    }
  }
  if (!files.details.empty()) writeDetailsText(joinPath(files.dir, files.details), rows, settings);
  if (!files.root.empty() && !writeDetailsRoot(joinPath(files.dir, files.root), rows, settings))
    std::cerr << "warning: built without ROOT, no ROOT file written\n";
  return 0;
}

// Public entry points: dispatch to the float (Fortran parity) or double instantiation.
int runUnpolarized(const RunInput& in, const Options& opt, const RunFiles& files, Precision prec) {
  return prec == Precision::Float ? runUnpolT<float>(in, opt, files) : runUnpolT<double>(in, opt, files);
}
int runPolarized(const RunInput& in, const Options& opt, const RunFiles& files, Precision prec) {
  return prec == Precision::Float ? runPolT<float>(in, opt, files) : runPolT<double>(in, opt, files);
}

// ---------------------------------------------------------------------------
// One row per point and configuration, fixed-width columns, settings in '#!' header lines.
void writeDetailsText(const std::string& file, const std::vector<DetailRow>& rows, const std::string& settings) {
  std::ofstream o(file);
  if (!o) throw std::runtime_error("cannot open " + file);
  o << "# radcor details (last pass). Energies MeV; cross sections in input units.\n"
       "# xs_in = measured; tail_Es / tail_Ep = E_s' / E_p' tail integrals (Fortran BBB / CCC);\n"
       "# soft = soft+virtual factor (CFACT); rad_model = soft*Born + tails; xs_born = unfolded Born;\n"
       "# prev = Born at pass NITT-1; RC = xs_in/xs_born; stat, syst = errors of xs_born; acc = acceptance factor\n";
  std::istringstream ss(settings);
  for (std::string l; std::getline(ss, l);) o << "#! " << l << "\n";
  o << "#spec conf  pt        E0        Es        Ep        nu        xs_in      tail_Es      tail_Ep"
       "         soft    rad_model      xs_born         prev           RC         stat         syst       acc\n";
  char buf[512];
  for (const auto& d : rows) {
    std::snprintf(buf, sizeof buf,
                  "%5d %4d %4d %9.2f %9.3f %9.3f %9.3f %12.5e %12.5e %12.5e %12.6f %12.5e %12.5e %12.5e %12.6f"
                  " %12.5e %12.5e %9.6f\n",
                  d.spectrum, d.config, d.point, d.e0, d.es, d.ep, d.nu, d.xs_in, d.tail_es, d.tail_ep,
                  d.soft_factor, d.rad_model, d.xs_born, d.xs_born_prev, d.rc, d.stat, d.syst, d.acc_factor);
    o << buf;
  }
}

}  // namespace radcor
