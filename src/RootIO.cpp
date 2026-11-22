// ROOT files: the per-point results (TTree "radcor" + TNamed "settings") and the input in
// format "root" (TTree "points" + TNamed "run").
#include "RootIO.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

#include "Details.h"

#ifdef RADCOR_HAS_ROOT
#include <TFile.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TTree.h>
#endif

namespace radcor {

// Write the rows as TTree "radcor" (one entry per point) and the settings as TNamed "settings".
// Returns false, writing nothing, if the library was built without ROOT.
bool writeDetailsRoot(const std::string& file, const std::vector<DetailRow>& rows, const std::string& settings) {
#ifdef RADCOR_HAS_ROOT
  TFile f(file.c_str(), "RECREATE");
  if (f.IsZombie()) throw std::runtime_error("cannot open ROOT output " + file);
  TTree t("radcor", "radcor per-point results (last iteration); energies MeV, xs in input units");
  DetailRow d;
  t.Branch("spectrum", &d.spectrum, "spectrum/I");
  t.Branch("config", &d.config, "config/I");
  t.Branch("point", &d.point, "point/I");
  t.Branch("e0", &d.e0, "e0/D");
  t.Branch("es", &d.es, "es/D");
  t.Branch("ep", &d.ep, "ep/D");
  t.Branch("nu", &d.nu, "nu/D");
  t.Branch("xs_in", &d.xs_in, "xs_in/D");
  t.Branch("tail_es", &d.tail_es, "tail_es/D");
  t.Branch("tail_ep", &d.tail_ep, "tail_ep/D");
  t.Branch("soft_factor", &d.soft_factor, "soft_factor/D");
  t.Branch("rad_model", &d.rad_model, "rad_model/D");
  t.Branch("xs_born", &d.xs_born, "xs_born/D");
  t.Branch("xs_born_prev", &d.xs_born_prev, "xs_born_prev/D");
  t.Branch("rc", &d.rc, "rc/D");
  t.Branch("stat", &d.stat, "stat/D");
  t.Branch("syst", &d.syst, "syst/D");
  t.Branch("acc_factor", &d.acc_factor, "acc_factor/D");
  for (const auto& r : rows) { d = r; t.Fill(); }
  t.Write();
  TNamed("settings", settings.c_str()).Write();
  f.Close();
  return true;
#else
  (void)file; (void)rows; (void)settings;
  return false;
#endif
}

#ifdef RADCOR_HAS_ROOT
namespace {
// Numeric branches of the input tree, in the CSV order (strings are not stored in ROOT).
const char* kCols[] = {"E0", "nu", "xs", "stat", "syst", "t_target", "t_before", "t_after",
                       "rc_sys", "norm", "loss_s", "loss_p"};
}  // namespace
#endif

// Read the input tree (any numeric branch type, via TLeaf::GetValue) into a RunInput by going
// through the CSV logic: rows are grouped into spectra by E0, spectrum columns must be constant.
RunInput readRootInput(const std::string& file, const std::string& kind) {
#ifdef RADCOR_HAS_ROOT
  TFile f(file.c_str(), "READ");
  if (f.IsZombie()) throw std::runtime_error("cannot open ROOT input " + file);
  auto* t = f.Get<TTree>("points");
  if (!t) throw std::runtime_error("no TTree 'points' in " + file);
  // Translate to the CSV text form and parse it with the CSV reader, so both formats follow
  // exactly the same rules (grouping, defaults, consistency checks).
  std::ostringstream csv;
  if (auto* run = f.Get<TNamed>("run")) csv << "#! " << run->GetTitle() << "\n";
  std::vector<std::string> names;
  std::vector<TLeaf*> leaves;
  for (const char* c : kCols)
    if (auto* l = t->GetLeaf(c)) { names.push_back(c); leaves.push_back(l); }
  TLeaf* cfg = t->GetLeaf("config");
  std::string acc, title, *accP = &acc, *titleP = &title;  // optional std::string branches
  const bool hasAcc = t->GetBranch("acc_table") && t->SetBranchAddress("acc_table", &accP) >= 0;
  const bool hasTitle = t->GetBranch("title") && t->SetBranchAddress("title", &titleP) >= 0;
  auto clean = [](std::string x) { std::replace(x.begin(), x.end(), ',', ' '); return x; };
  for (size_t k = 0; k < names.size(); ++k) csv << (k ? "," : "") << names[k];
  if (cfg) csv << ",config";
  if (hasAcc) csv << ",acc_table";
  if (hasTitle) csv << ",title";
  csv << "\n";
  csv.precision(17);
  for (Long64_t e = 0; e < t->GetEntries(); ++e) {
    t->GetEntry(e);
    for (size_t k = 0; k < leaves.size(); ++k) {
      const double v = leaves[k]->GetValue();
      csv << (k ? "," : "");
      if (!std::isnan(v)) csv << v;
    }
    if (cfg) csv << "," << static_cast<int>(cfg->GetValue());
    if (hasAcc) csv << "," << clean(acc);
    if (hasTitle) csv << "," << clean(title);
    csv << "\n";
  }
  std::istringstream is(csv.str());
  return parseCsv(is, file, kind);
#else
  (void)file; (void)kind;
  throw std::runtime_error("format root: radcor was built without ROOT");
#endif
}

// Write the input as TTree "points" and TNamed "run".
void writeRootInput(const RunInput& in, const std::string& file) {
#ifdef RADCOR_HAS_ROOT
  TFile f(file.c_str(), "RECREATE");
  if (f.IsZombie()) throw std::runtime_error("cannot open ROOT output " + file);
  TTree t("points", "radcor input: one entry per point; energies MeV, thicknesses r.l.");
  double v[12];
  int config = 1;
  for (int k = 0; k < 12; ++k) t.Branch(kCols[k], &v[k], (std::string(kCols[k]) + "/D").c_str());
  t.Branch("config", &config, "config/I");
  std::string acc, title;
  t.Branch("acc_table", &acc);
  t.Branch("title", &title);
  for (const auto& s : in.spectra)
    for (size_t c = 0; c < s.configs.size(); ++c) {
      const auto& d = s.configs[c];
      for (size_t j = 0; j < d.nu.size(); ++j) {
        const double row[12] = {s.e0, d.nu[j], d.xs[j], d.stat.empty() ? 0 : d.stat[j],
                                d.syst.empty() ? 0 : d.syst[j], s.tTarget, s.tBefore, s.tAfter,
                                s.rcSys, s.norm, s.lossS, s.lossP};
        std::copy(row, row + 12, v);
        config = static_cast<int>(c) + 1;
        acc = s.accTable;
        title = s.title;
        t.Fill();
      }
    }
  t.Write();
  std::ostringstream run;
  run.precision(17);
  run << "kind=" << in.kind << " nitt=" << in.nitt << " Z=" << in.Z << " A=" << in.A << " theta=" << in.theta
      << " ZW=" << in.ZW << " AW=" << in.AW << " ttheta=" << in.tTheta;
  TNamed("run", run.str().c_str()).Write();
  f.Close();
#else
  (void)in; (void)file;
  throw std::runtime_error("format root: radcor was built without ROOT");
#endif
}

}  // namespace radcor
