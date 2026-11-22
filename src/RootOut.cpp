// Optional ROOT output: a TTree "radcor" with one entry per point (last iteration),
// and a TNamed "settings" with the complete run configuration.
#include "Details.h"

#include <stdexcept>

#ifdef RADCOR_HAS_ROOT
#include <TFile.h>
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

}  // namespace radcor
