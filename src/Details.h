// Per-point output record shared by the text "details" writer and the ROOT writer.
#pragma once
#include <string>
#include <vector>

namespace radcor {

struct DetailRow {
  int spectrum = 0;          // 1-based spectrum index
  int config = 1;            // 1 or 2 (polarized), 1 for unpolarized
  int point = 0;             // 1-based point index
  double e0 = 0;             // beam energy as given [MeV]
  double es = 0, ep = 0;     // kinematics used (after collisional loss) [MeV]
  double nu = 0;             // es - ep [MeV]
  double xs_in = 0;          // measured (radiated) xs, x RNORME
  double tail_es = 0;        // E_s' tail integral (Fortran BBB)
  double tail_ep = 0;        // E_p' tail integral (Fortran CCC)
  double soft_factor = 0;    // soft+virtual factor (Fortran CFACT)
  double rad_model = 0;      // radiated model soft_factor*Born + tails, last pass (Fortran CRS)
  double xs_born = 0;        // unfolded Born xs at these kinematics
  double xs_born_prev = 0;   // Born at pass NITT-1 (convergence check)
  double rc = 0;             // xs_in / xs_born
  double stat = 0;           // statistical error of xs_born (unpol: input stat x xs_born/xs_in; pol: input error)
  double syst = 0;           // systematic error of xs_born (unpol: as the Fortran SYST, incl. RC_SYS; pol: 0)
  double acc_factor = 1;     // finite-acceptance factor at this nu (1 if finite_acc=0); XSfc = xs_born*acc_factor
};

// `settings` is written as a comment header (text) / a TNamed "settings" (ROOT).
void writeDetailsText(const std::string& file, const std::vector<DetailRow>& rows, const std::string& settings);
bool writeDetailsRoot(const std::string& file, const std::vector<DetailRow>& rows,
                      const std::string& settings);  // false if built without ROOT

}  // namespace radcor
