// Run options: parsing of "--set key=value", printing, and the dilogarithm used by spenceExact.
#include "radcor/Options.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace radcor {

// "1/true/on/yes" -> true, "0/false/off/no" -> false; anything else is an error.
static bool toBool(const std::string& v) {
  if (v == "1" || v == "true" || v == "on" || v == "yes") return true;
  if (v == "0" || v == "false" || v == "off" || v == "no") return false;
  throw std::runtime_error("bad boolean '" + v + "'");
}
// Comma-separated numbers, e.g. "0.4,0.5,0.5".
static std::vector<double> toList(const std::string& v) {
  std::vector<double> out;
  std::stringstream ss(v);
  std::string t;
  while (std::getline(ss, t, ',')) out.push_back(std::stod(t));
  return out;
}

void setOption(Options& o, const std::string& kv) {
  auto eq = kv.find('=');
  if (eq == std::string::npos) throw std::runtime_error("option must be key=value: " + kv);
  const std::string k = kv.substr(0, eq), v = kv.substr(eq + 1);
  if (k == "internal") o.internal = toBool(v);
  else if (k == "external") o.external = toBool(v);
  else if (k == "collisional_loss") o.collisional_loss = toBool(v);
  else if (k == "finite_acc") o.finite_acc = toBool(v);
  else if (k == "radiate") o.radiate = toBool(v);
  else if (k == "delta" || k == "cute") o.delta = std::stod(v);
  else if (k == "simpson_min" || k == "num") o.simpson_min = std::stoi(v);
  else if (k == "simpson_step" || k == "del") o.simpson_step = std::stod(v);
  else if (k == "nitt") o.nitt_override = std::stoi(v);
  else if (k == "amass") o.amass = std::stod(v);
  else if (k == "deltas") o.deltas = toList(v);
  else if (k == "deltap") o.deltap = toList(v);
  else if (k == "del_overwrite") o.del_overwrite = toBool(v);
  else if (k == "spence_legacy") o.spence_legacy = toBool(v);
  else if (k == "b") o.b_stein = std::stod(v);
  else if (k == "ejsz_cap") o.ejsz_cap = toBool(v);
  else if (k == "spline_bug") o.spline_bug = toBool(v);
  else if (k == "spline_y2zero") o.spline_y2zero = std::stod(v);
  else if (k == "tbta_lag") o.tbta_lag = toBool(v);
  else if (k == "radiate_overwrite") o.radiate_overwrite = toBool(v);
  else if (k == "acc_cells") {
    o.acc_cells.clear();
    std::stringstream ss(v);
    std::string t;
    while (std::getline(ss, t, ',')) o.acc_cells.push_back(t);
  } else if (k == "legacy") {  // legacy=0 switches off every known bug at once
    const bool L = toBool(v);
    o.del_overwrite = L; o.spence_legacy = L; o.spline_bug = L; o.tbta_lag = L; o.radiate_overwrite = L;
  } else throw std::runtime_error("unknown option '" + k + "'");
}

std::string describeOptions(const Options& o) {
  std::ostringstream s;
  auto b = [](bool x) { return x ? "1" : "0"; };
  s << "internal=" << b(o.internal) << " external=" << b(o.external) << " collisional_loss=" << b(o.collisional_loss)
    << " finite_acc=" << b(o.finite_acc) << " radiate=" << b(o.radiate) << " delta=" << o.delta << " simpson_min=" << o.simpson_min
    << " simpson_step=" << o.simpson_step << " amass=" << o.amass << " b=" << o.b_stein << "\n  legacy: del_overwrite="
    << b(o.del_overwrite) << " spence_legacy=" << b(o.spence_legacy) << " ejsz_cap=" << b(o.ejsz_cap)
    << " spline_bug=" << b(o.spline_bug) << " tbta_lag=" << b(o.tbta_lag)
    << " radiate_overwrite=" << b(o.radiate_overwrite);
  return s.str();
}

// Real dilogarithm Li2(x) for x <= 1 (Mo-Tsai II.8 / Stein A48 Spence function).
double dilog(double x) {
  const double pi2_6 = M_PI * M_PI / 6.0;
  if (x == 1.0) return pi2_6;
  if (x == 0.0) return 0.0;
  if (x > 1.0) throw std::domain_error("dilog: x > 1 not supported");
  if (x < 0.0) {  // Li2(x) = -Li2(x/(x-1)) - 1/2 ln^2(1-x), x/(x-1) in (0,1)
    const double l = std::log(1.0 - x);
    return -dilog(x / (x - 1.0)) - 0.5 * l * l;
  }
  if (x > 0.5) return pi2_6 - std::log(x) * std::log(1.0 - x) - dilog(1.0 - x);
  double s = 0, t = x;
  for (int n = 1; n < 200; ++n) {
    const double term = t / (double(n) * n);
    s += term;
    if (term < 1e-18 * s) break;
    t *= x;
  }
  return s;
}

}  // namespace radcor
