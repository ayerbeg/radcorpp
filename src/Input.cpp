// Input formats (ele-input, polarized pair, CSV, ROOT), validation and summary of a RunInput.
#include "radcor/Input.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "RootIO.h"

namespace radcor {

namespace {

// ---------------------------------------------------------------------------
// Small helpers

// Shortest decimal text that reads back as exactly the same double.
std::string num(double v) {
  if (std::isnan(v)) return "";
  char buf[64];
  auto r = std::to_chars(buf, buf + sizeof buf, v);
  return std::string(buf, r.ptr);
}

std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  const auto b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

// Parse a whole field as a double (Fortran 'D' exponents accepted); throws with context.
double toDouble(std::string t, const std::string& where) {
  for (auto& c : t) if (c == 'd' || c == 'D') c = 'E';
  t = trim(t);
  char* end = nullptr;
  const double v = std::strtod(t.c_str(), &end);
  if (t.empty() || end != t.c_str() + t.size()) throw std::runtime_error("bad number '" + t + "' " + where);
  return v;
}

// Fortran list-directed reader: numbers separated by blanks/commas, spanning records; the
// remainder of a record after the last needed value is skipped, as in Fortran.
class ListReader {
 public:
  explicit ListReader(const std::string& file) : in_(file), file_(file) {
    if (!in_) throw std::runtime_error("cannot open " + file);
  }
  // Next n numbers, reading as many lines as needed; the rest of the last line is discarded.
  std::vector<double> read(int n) {
    std::vector<double> v;
    std::string line;
    while (static_cast<int>(v.size()) < n) {
      if (!std::getline(in_, line)) throw std::runtime_error("unexpected end of " + file_);
      ++lineNo_;
      for (auto& c : line) if (c == ',') c = ' ';
      std::istringstream ss(line);
      std::string tok;
      while (static_cast<int>(v.size()) < n && ss >> tok) {
        if (tok[0] == '!' || tok[0] == '/') break;
        for (auto& c : tok) if (c == 'd' || c == 'D') c = 'E';
        try { v.push_back(std::stod(tok)); }
        catch (...) { throw std::runtime_error("bad number '" + tok + "' at " + file_ + ":" + std::to_string(lineNo_)); }
      }
    }
    return v;
  }
  // Next whole line (the spectrum title, which is not numeric).
  std::string readLine() {
    std::string line;
    if (!std::getline(in_, line)) throw std::runtime_error("unexpected end of " + file_);
    ++lineNo_;
    return line;
  }

 private:
  std::ifstream in_;
  std::string file_;
  int lineNo_ = 0;
};

// Header shared by both Fortran formats (U:242-262): NITT, NE, Z A THETA, ZW AW TTHETA, and
// NE lines "TX TIW TFW RC_SYS RNORME". Returns NE; thickness lines go to `thick`.
int readFortranHeader(ListReader& r, RunInput& in, std::vector<std::vector<double>>& thick) {
  in.nitt = static_cast<int>(r.read(1)[0]);
  const int ne = static_cast<int>(r.read(1)[0]);
  auto v = r.read(3); in.Z = v[0]; in.A = v[1]; in.theta = v[2];
  v = r.read(3); in.ZW = v[0]; in.AW = v[1]; in.tTheta = v[2];
  if (ne < 1) throw std::runtime_error("NE must be >= 1");
  for (int k = 0; k < ne; ++k) thick.push_back(r.read(5));
  return ne;
}

void setThickness(SpectrumInput& s, const std::vector<double>& t) {
  s.tTarget = t[0]; s.tBefore = t[1]; s.tAfter = t[2]; s.rcSys = t[3]; s.norm = t[4];
}

void writeFortranHeader(std::ostream& o, const RunInput& in) {
  o << num(in.nitt) << "                ! NITT   number of unfolding passes\n"
    << in.spectra.size() << "                 ! NE     number of spectra\n"
    << num(in.Z) << " " << num(in.A) << " " << num(in.theta) << "          ! Z A THETA[deg]\n"
    << num(in.ZW) << " " << num(in.AW) << " " << num(in.tTheta) << "          ! ZW AW TTHETA (not used)\n";
  for (const auto& s : in.spectra)
    o << num(s.tTarget) << " " << num(s.tBefore) << " " << num(s.tAfter) << " " << num(s.rcSys) << " "
      << num(s.norm) << "   ! TX TIW TFW RC_SYS RNORME  E0=" << num(s.e0) << "\n";
}

// Per-spectrum settings that the Fortran formats cannot hold are dropped with a warning.
void warnDropped(const RunInput& in, const char* fmt) {
  for (const auto& s : in.spectra)
    if (!std::isnan(s.lossS) || !std::isnan(s.lossP) || !s.accTable.empty()) {
      std::cerr << "warning: format " << fmt << " cannot store per-spectrum loss_s/loss_p/acc_table; "
                   "they are dropped (give them as options instead)\n";
      return;
    }
}

// Proton-mass invariant mass W [MeV], used only to fill the informative W column of ds.inp.
double invariantMass(double e0, double thetaDeg, double nu) {
  const double M = 938.272, s2 = std::pow(std::sin(thetaDeg * M_PI / 360.0), 2);
  const double w2 = M * M - 4 * e0 * (e0 - nu) * s2 + 2 * M * nu;
  return w2 > 0 ? std::sqrt(w2) : 0.0;
}

// ---------------------------------------------------------------------------
// "ele": unpolarized Fortran input
class EleFormat : public InputFormat {
 public:
  std::string name() const override { return "ele"; }

  RunInput read(const InputFiles& f, const std::string& kind) const override {
    if (kind == "pol") throw std::runtime_error("format ele is unpolarized; use format pol, csv or root");
    RunInput in;
    in.kind = "unpol";
    in.source = f.main;
    ListReader r(f.main);
    std::vector<std::vector<double>> thick;
    const int ne = readFortranHeader(r, in, thick);
    for (int i = 0; i < ne; ++i) {
      SpectrumInput s;
      s.title = trim(r.readLine());  // TITLE (A7)
      auto hd = r.read(4);           // ES0 EMX DE NP
      s.e0 = hd[0];
      const int np = static_cast<int>(hd[3]);
      setThickness(s, thick[i]);
      ConfigData c;
      for (int j = 0; j < np; ++j) {
        auto p = r.read(4);
        c.nu.push_back(p[0]); c.xs.push_back(p[1]); c.stat.push_back(p[2]); c.syst.push_back(p[3]);
      }
      s.configs.push_back(std::move(c));
      in.spectra.push_back(std::move(s));
    }
    return in;
  }

  void write(const RunInput& in, const InputFiles& f) const override {
    if (in.kind != "unpol") throw std::runtime_error("format ele holds unpolarized input only");
    warnDropped(in, "ele");
    if (in.spectra.size() > 13) std::cerr << "warning: NE > 13: the Fortran program cannot read this file\n";
    std::ofstream o(f.main);
    if (!o) throw std::runtime_error("cannot write " + f.main);
    writeFortranHeader(o, in);
    for (size_t i = 0; i < in.spectra.size(); ++i) {
      const auto& s = in.spectra[i];
      const auto& c = s.configs.at(0);
      if (c.nu.size() > 1499) std::cerr << "warning: NP > 1499: the Fortran program cannot read this file\n";
      o << (s.title.empty() ? "spec" + std::to_string(i + 1) : s.title) << "\n";
      o << num(s.e0) << ", 0, 0, " << c.nu.size() << "      ! ES0 EMX DE NP\n";
      for (size_t j = 0; j < c.nu.size(); ++j)
        o << num(c.nu[j]) << " " << num(c.xs[j]) << " " << num(c.stat.empty() ? 0.0 : c.stat[j]) << " "
          << num(c.syst.empty() ? 0.0 : c.syst[j]) << "\n";
    }
  }
};

// ---------------------------------------------------------------------------
// "pol": polarized Fortran input, header file + ds.inp table
class PolFormat : public InputFormat {
 public:
  std::string name() const override { return "pol"; }

  RunInput read(const InputFiles& f, const std::string& kind) const override {
    if (kind == "unpol") throw std::runtime_error("format pol is polarized; use format ele, csv or root");
    if (f.data.empty()) throw std::runtime_error("format pol needs the data table (ds.inp)");
    RunInput in;
    in.kind = "pol";
    in.source = f.main + " + " + f.data;
    ListReader r(f.main);
    std::vector<std::vector<double>> thick;
    const int ne = readFortranHeader(r, in, thick);
    // ds.inp: E0 NU W DS1 ES1 DS2 ES2; a new spectrum when E0 changes (P:253-317).
    // A point belongs to configuration k only if its error is >= 0.
    std::ifstream din(f.data);
    if (!din) throw std::runtime_error("cannot open " + f.data);
    std::string line;
    double e0Prev = 0;
    while (std::getline(din, line)) {
      if (!line.empty() && line[0] == '#') continue;
      std::istringstream ss(line);
      double e0, nu, w, d1, e1, d2, e2;
      if (!(ss >> e0 >> nu >> w >> d1 >> e1 >> d2 >> e2)) break;  // Fortran ERR=/END= exit
      if (in.spectra.empty() || static_cast<float>(e0) != static_cast<float>(e0Prev)) {
        in.spectra.emplace_back();
        in.spectra.back().e0 = e0;
        in.spectra.back().configs.resize(2);
        e0Prev = e0;
      }
      auto& s = in.spectra.back();
      const double d[2] = {d1, d2}, e[2] = {e1, e2};
      for (int k = 0; k < 2; ++k)
        if (e[k] >= 0) {
          s.configs[k].nu.push_back(nu); s.configs[k].xs.push_back(d[k]); s.configs[k].stat.push_back(e[k]);
        }
    }
    const int nsp = static_cast<int>(in.spectra.size());
    if (nsp != ne)
      std::cerr << "warning: NE=" << ne << " in " << f.main << " but " << nsp << " spectra in " << f.data
                << " (the Fortran does not check)\n";
    if (nsp > ne) throw std::runtime_error("more spectra in the data table than thickness lines in the header");
    for (int i = 0; i < nsp; ++i) setThickness(in.spectra[i], thick[i]);
    return in;
  }

  void write(const RunInput& in, const InputFiles& f) const override {
    if (in.kind != "pol") throw std::runtime_error("format pol holds polarized input only");
    if (f.data.empty()) throw std::runtime_error("format pol needs a second file name for the data table");
    warnDropped(in, "pol");
    std::ofstream h(f.main);
    if (!h) throw std::runtime_error("cannot write " + f.main);
    writeFortranHeader(h, in);
    std::ofstream o(f.data);
    if (!o) throw std::runtime_error("cannot write " + f.data);
    o << "#E0 NU W DeltaSig1 Error1 DeltaSig2 Error2   (error < 0: point absent in that configuration)\n";
    for (const auto& s : in.spectra) {
      std::set<double> all;  // union of the nu lists of the two configurations
      for (const auto& c : s.configs) all.insert(c.nu.begin(), c.nu.end());
      for (double nu : all) {
        o << num(s.e0) << " " << num(nu) << " " << num(invariantMass(s.e0, in.theta, nu));
        for (const auto& c : s.configs) {
          auto it = std::find(c.nu.begin(), c.nu.end(), nu);
          if (it == c.nu.end()) { o << " 0 -1"; continue; }
          const size_t j = it - c.nu.begin();
          o << " " << num(c.xs[j]) << " " << num(c.stat[j]);
        }
        o << "\n";
      }
    }
  }
};

// ---------------------------------------------------------------------------
// "csv": named columns, one row per point.
//   #! kind=unpol nitt=10 Z=2 A=3 theta=15.5        run settings (any number of "#!" lines)
//   E0,nu,xs,stat,syst,t_before,t_after             header: the columns present, any order
//   862,40.5,1.2e3,...                              one row per point
// Required columns: E0, nu, xs. Optional: config (1|2, polarized), stat, syst, t_target, t_before,
// t_after, rc_sys, norm, loss_s, loss_p, acc_table, title. Spectrum columns (t_*, rc_sys, norm,
// loss_*, acc_table, title) must be the same on every row of a spectrum. An empty field = default.
class CsvFormat : public InputFormat {
 public:
  std::string name() const override { return "csv"; }

  RunInput read(const InputFiles& f, const std::string& kind) const override {
    std::ifstream is(f.main);
    if (!is) throw std::runtime_error("cannot open " + f.main);
    return parse(is, f.main, kind);
  }

  // The CSV parser proper, on any stream (also used by the ROOT reader).
  static RunInput parse(std::istream& is, const std::string& source, const std::string& kind) {
    InputFiles f{source, ""};
    RunInput in;
    in.source = source;
    in.kind = kind.empty() ? "unpol" : kind;
    bool kindGiven = false;
    std::vector<std::string> cols;
    std::map<double, size_t> index;  // E0 -> spectrum, in order of first appearance
    std::string line;
    int lineNo = 0;
    static const std::set<std::string> known = {"E0", "nu", "xs", "config", "stat", "syst", "t_target",
                                                "t_before", "t_after", "rc_sys", "norm", "loss_s",
                                                "loss_p", "acc_table", "title"};
    while (std::getline(is, line)) {
      ++lineNo;
      const std::string where = "at " + f.main + ":" + std::to_string(lineNo);
      if (line.rfind("#!", 0) == 0) {  // run settings
        std::istringstream ss(line.substr(2));
        std::string kv;
        while (ss >> kv) {
          const auto eq = kv.find('=');
          if (eq == std::string::npos) throw std::runtime_error("expected key=value " + where);
          const std::string k = kv.substr(0, eq), v = kv.substr(eq + 1);
          if (k == "kind") { in.kind = v; kindGiven = true; }
          else if (k == "nitt") in.nitt = static_cast<int>(toDouble(v, where));
          else if (k == "Z") in.Z = toDouble(v, where);
          else if (k == "A") in.A = toDouble(v, where);
          else if (k == "theta") in.theta = toDouble(v, where);
          else if (k == "ZW") in.ZW = toDouble(v, where);
          else if (k == "AW") in.AW = toDouble(v, where);
          else if (k == "ttheta") in.tTheta = toDouble(v, where);
          else throw std::runtime_error("unknown run setting '" + k + "' " + where);
        }
        continue;
      }
      if (trim(line).empty() || line[0] == '#') continue;
      std::vector<std::string> fld;
      { std::stringstream ss(line); std::string t; while (std::getline(ss, t, ',')) fld.push_back(trim(t)); }
      if (line.back() == ',') fld.push_back("");
      if (cols.empty()) {  // header line
        cols = fld;
        for (const auto& c : cols)
          if (!known.count(c)) throw std::runtime_error("unknown column '" + c + "' " + where);
        for (const char* c : {"E0", "nu", "xs"})
          if (std::find(cols.begin(), cols.end(), c) == cols.end())
            throw std::runtime_error(std::string("missing column ") + c + " in " + f.main);
        if (!kindGiven && !kind.empty()) in.kind = kind;
        continue;
      }
      if (fld.size() != cols.size()) throw std::runtime_error("wrong number of fields " + where);
      std::map<std::string, std::string> row;
      for (size_t k = 0; k < cols.size(); ++k) row[cols[k]] = fld[k];
      auto get = [&](const char* k, double def) {
        auto it = row.find(k);
        return it == row.end() || it->second.empty() ? def : toDouble(it->second, where);
      };
      const double e0 = get("E0", 0);
      auto it = index.find(e0);
      const bool fresh = it == index.end();
      if (fresh) {
        it = index.emplace(e0, in.spectra.size()).first;
        in.spectra.emplace_back();
        auto& s = in.spectra.back();
        s.e0 = e0;
        s.configs.resize(in.nconfig());
      }
      SpectrumInput& s = in.spectra[it->second];
      // spectrum-level columns: set from the first row, then required to be equal
      auto spec = [&](const char* k, double& field, double def) {
        const double v = get(k, def);
        if (fresh) field = v;
        else if (!(v == field || (std::isnan(v) && std::isnan(field))))
          throw std::runtime_error(std::string("column ") + k + " changes within spectrum E0=" + num(e0) + " " + where);
      };
      spec("t_target", s.tTarget, 0); spec("t_before", s.tBefore, 0); spec("t_after", s.tAfter, 0);
      spec("rc_sys", s.rcSys, 0); spec("norm", s.norm, 1);
      spec("loss_s", s.lossS, SpectrumInput::unset); spec("loss_p", s.lossP, SpectrumInput::unset);
      for (const char* k : {"acc_table", "title"}) {
        auto r = row.find(k);
        const std::string v = r == row.end() ? "" : r->second;
        std::string& field = std::string(k) == "title" ? s.title : s.accTable;
        if (fresh) field = v;
        else if (v != field) throw std::runtime_error(std::string("column ") + k + " changes within spectrum " + where);
      }
      const int cfg = static_cast<int>(get("config", 1));
      if (cfg < 1 || cfg > in.nconfig())
        throw std::runtime_error("config must be 1" + std::string(in.nconfig() == 2 ? " or 2 " : " ") + where);
      auto& c = s.configs[cfg - 1];
      c.nu.push_back(get("nu", 0)); c.xs.push_back(get("xs", 0));
      c.stat.push_back(get("stat", 0)); c.syst.push_back(get("syst", 0));
    }
    if (cols.empty()) throw std::runtime_error("no header line in " + f.main);
    if (!kind.empty() && in.kind != kind)
      throw std::runtime_error(f.main + " is kind=" + in.kind + " but a " + kind + " run was requested");
    return in;
  }

  void write(const RunInput& in, const InputFiles& f) const override {
    std::ofstream o(f.main);
    if (!o) throw std::runtime_error("cannot write " + f.main);
    o << "# radcor input table: one row per point. Energies MeV, thicknesses in radiation lengths.\n"
      << "#! kind=" << in.kind << " nitt=" << in.nitt << " Z=" << num(in.Z) << " A=" << num(in.A)
      << " theta=" << num(in.theta) << " ZW=" << num(in.ZW) << " AW=" << num(in.AW) << " ttheta=" << num(in.tTheta)
      << "\nE0,config,nu,xs,stat,syst,t_target,t_before,t_after,rc_sys,norm,loss_s,loss_p,acc_table,title\n";
    auto clean = [](std::string s) { std::replace(s.begin(), s.end(), ',', ' '); return s; };
    for (const auto& s : in.spectra)
      for (size_t k = 0; k < s.configs.size(); ++k) {
        const auto& c = s.configs[k];
        for (size_t j = 0; j < c.nu.size(); ++j)
          o << num(s.e0) << "," << k + 1 << "," << num(c.nu[j]) << "," << num(c.xs[j]) << ","
            << num(c.stat.empty() ? 0.0 : c.stat[j]) << "," << num(c.syst.empty() ? 0.0 : c.syst[j]) << ","
            << num(s.tTarget) << "," << num(s.tBefore) << "," << num(s.tAfter) << "," << num(s.rcSys) << ","
            << num(s.norm) << "," << num(s.lossS) << "," << num(s.lossP) << "," << clean(s.accTable) << ","
            << clean(s.title) << "\n";
      }
  }
};

// ---------------------------------------------------------------------------
// "root": TTree "points" (the CSV columns as branches) + TNamed "run" (the "#!" settings).
class RootFormat : public InputFormat {
 public:
  std::string name() const override { return "root"; }
  RunInput read(const InputFiles& f, const std::string& kind) const override {
    RunInput in = readRootInput(f.main, kind);
    in.source = f.main;
    return in;
  }
  void write(const RunInput& in, const InputFiles& f) const override { writeRootInput(in, f.main); }
};

}  // namespace

// ---------------------------------------------------------------------------
RunInput parseCsv(std::istream& is, const std::string& source, const std::string& kind) {
  return CsvFormat::parse(is, source, kind);
}

const InputFormat& inputFormat(const std::string& name) {
  static const EleFormat ele;
  static const PolFormat pol;
  static const CsvFormat csv;
  static const RootFormat root;
  if (name == "ele") return ele;
  if (name == "pol") return pol;
  if (name == "csv") return csv;
  if (name == "root") return root;
  throw std::runtime_error("unknown input format '" + name + "' (ele, pol, csv, root)");
}

std::string guessFormat(const std::string& file, const std::string& kind) {
  auto ends = [&](const std::string& e) {
    return file.size() >= e.size() && file.compare(file.size() - e.size(), e.size(), e) == 0;
  };
  if (ends(".csv") || ends(".CSV")) return "csv";
  if (ends(".root")) return "root";
  return kind == "pol" ? "pol" : "ele";
}

std::vector<std::string> validate(const RunInput& in) {
  std::vector<std::string> p;
  auto bad = [](double x) { return !std::isfinite(x); };
  if (in.kind != "unpol" && in.kind != "pol") p.push_back("kind must be unpol or pol, not '" + in.kind + "'");
  if (in.nitt < 1) p.push_back("nitt must be >= 1");
  if (bad(in.Z) || in.Z <= 0) p.push_back("target Z must be > 0");
  if (bad(in.A) || in.A <= 0) p.push_back("target A must be > 0");
  if (bad(in.theta) || in.theta <= 0 || in.theta >= 180) p.push_back("theta must be in (0, 180) deg");
  if (in.spectra.empty()) p.push_back("no spectra");
  for (size_t i = 0; i < in.spectra.size(); ++i) {
    const auto& s = in.spectra[i];
    const std::string id = "spectrum " + std::to_string(i + 1) + " (E0=" + num(s.e0) + "): ";
    if (bad(s.e0) || s.e0 <= 0) p.push_back(id + "E0 must be > 0");
    if (i > 0 && !(s.e0 > in.spectra[i - 1].e0)) p.push_back(id + "spectra must be in strictly increasing E0");
    for (double t : {s.tTarget, s.tBefore, s.tAfter})
      if (bad(t) || t < 0) { p.push_back(id + "thicknesses must be >= 0"); break; }
    if (bad(s.norm) || s.norm <= 0) p.push_back(id + "norm must be > 0");
    if (bad(s.rcSys) || s.rcSys < 0) p.push_back(id + "rc_sys must be >= 0");
    for (double l : {s.lossS, s.lossP})
      if (std::isinf(l)) p.push_back(id + "collisional loss must be finite");
    if (static_cast<int>(s.configs.size()) != in.nconfig()) {
      p.push_back(id + "expected " + std::to_string(in.nconfig()) + " configuration(s)");
      continue;
    }
    for (size_t k = 0; k < s.configs.size(); ++k) {
      const auto& c = s.configs[k];
      const std::string cid = id + (in.nconfig() == 2 ? "config " + std::to_string(k + 1) + ": " : "");
      const size_t n = c.nu.size();
      if (n < 2) { p.push_back(cid + "at least 2 points needed"); continue; }
      if (c.xs.size() != n || c.stat.size() != n || (!c.syst.empty() && c.syst.size() != n)) {
        p.push_back(cid + "columns of different length");
        continue;
      }
      for (size_t j = 0; j < n; ++j) {
        if (bad(c.nu[j]) || bad(c.xs[j]) || bad(c.stat[j]) || (!c.syst.empty() && bad(c.syst[j]))) {
          p.push_back(cid + "non-finite value at point " + std::to_string(j + 1));
          break;
        }
        if (j > 0 && !(c.nu[j] > c.nu[j - 1])) {
          p.push_back(cid + "nu must increase strictly (point " + std::to_string(j + 1) + ", nu=" + num(c.nu[j]) + ")");
          break;
        }
      }
      if (c.nu.back() >= s.e0) p.push_back(cid + "nu must be < E0");
    }
    if (in.nconfig() == 2 && s.configs[0].nu != s.configs[1].nu)
      p.push_back(id + "the two configurations must have the same nu points (Fortran limitation B-7)");
  }
  return p;
}

void checkInput(const RunInput& in) {
  const auto p = validate(in);
  if (p.empty()) return;
  std::string msg = "invalid input " + in.source + ":";
  for (const auto& s : p) msg += "\n  " + s;
  throw std::runtime_error(msg);
}

std::string describeInput(const RunInput& in) {
  std::ostringstream o;
  o << in.kind << " input " << in.source << "\n  Z=" << num(in.Z) << " A=" << num(in.A) << " theta=" << num(in.theta)
    << " nitt=" << in.nitt << " spectra=" << in.spectra.size() << "\n";
  char buf[256];
  for (const auto& s : in.spectra) {
    const auto& c = s.configs.at(0);
    std::snprintf(buf, sizeof buf, "  E0=%8.1f  points=%5zu  nu=[%.1f, %.1f]  t_target=%g t_before=%g t_after=%g norm=%g",
                  s.e0, c.nu.size(), c.nu.empty() ? 0.0 : c.nu.front(), c.nu.empty() ? 0.0 : c.nu.back(), s.tTarget,
                  s.tBefore, s.tAfter, s.norm);
    o << buf << (s.title.empty() ? "" : "  [" + s.title + "]") << "\n";
  }
  return o.str();
}

}  // namespace radcor
