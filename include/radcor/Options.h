// radcor — C++ translation of K. Slifer's nradcor.f (2005), itself derived from
// L. Whitney's radcor.f. Physics: Stein et al. PRD 12, 1884 (1975) Eq. A82,
// iterative unfolding of Mo & Tsai RMP 41, 205 (1969) Eq. IV.2.
//
// All switches that were hard-coded in the Fortran source are collected here.
// Defaults reproduce the Fortran exactly ("legacy" behaviour), including its known
// bugs; each bug can be switched off individually. The bug labels B-n refer to the
// chapter "Known issues" of the manual.
#pragma once
#include <string>
#include <vector>

namespace radcor {

// Every setting of a run. Filled with defaults, then modified by setOption() ("--set key=value").
struct Options {
  // --- physics switches (U:151-158) -------------------------------------------
  bool internal = true;          // internal RC (t_r, vertex/vacuum F-bar)
  bool external = true;          // external RC (t_b, t_a, xi)
  bool collisional_loss = true;  // E94-010 specific most-probable collisional loss shift
  bool finite_acc = true;        // E94-010 specific finite acceptance factor (QFS tables)
  bool radiate = false;          // "radiate a Born xs" mode (needs nitt=1, one spectrum)

  // --- numerical parameters (U:163-167) ----------------------------------------
  double delta = 10.0;           // Delta [MeV], soft/hard photon split          (CUTE)
  int simpson_min = 80;          // Simpson minimum number of intervals          (NUM)
  double simpson_step = 1.0;     // Simpson maximum step [MeV]                   (DEL)
  int nitt_override = -1;        // >0 overrides NITT from the input file

  // --- target mass per nucleon (U: 931.494, P: 938.272) ------------------------
  double amass = 931.494;

  // --- E94-010 collisional losses (U:125-126), MeV, by spectrum index ----------
  std::vector<double> deltas{0.4, 0.5, 0.5, 0.5, 0.5, 0.5};
  std::vector<double> deltap{2.4, 3.3, 2.2, 3.3, 2.5, 2.2};

  // --- legacy-compatibility switches (true = reproduce the Fortran) -------------
  bool del_overwrite = true;     // B-2: DEL=CUTE after the first point (U:494)
  bool spence_legacy = true;     // B-1: Li2(sin^2 th/2) via 50-term series (U:383-389);
                                 //      false = Stein A44 Phi(cos^2 th/2) (exact Li2)
  double b_stein = 4.0 / 3.0;    // b used in the Stein kernels (U:495 "4./3. ! A45")
  bool ejsz_cap = false;         // Whitney's FEP upper limit cap E_s/eta_s - 3 MeV (W:476-478),
                                 // commented out in the shipped U code
  bool spline_bug = true;        // B-3: NR spline with Y2(I-2) instead of Y2(I-1) (U:1491)
  double spline_y2zero = 0.0;    // value read at the out-of-bounds Y2(0) (gfortran build: 0.0)
  bool tbta_lag = true;          // B-9: Tb/Ta (COMMON) are set after the integrals of a point
                                 //      (U:497-503), so FES/FEP see the previous point's values
  bool radiate_overwrite = true; // B-4: in RADIATE mode the top spectrum model is overwritten

  // --- files ---------------------------------------------------------------------
  std::vector<std::string> acc_cells{"8_SY", "1_AR", "2_JN", "3_AR", "4_2T", "5_JN"};
};

// Parse one "key=value" and set the corresponding field of opts. Fortran-style aliases are
// accepted: cute=delta, num=simpson_min, del=simpson_step. Throws std::runtime_error on an
// unknown key or a malformed value.
void setOption(Options& opts, const std::string& keyval);

// All settings as "key=value" text (two lines), for the log and the output headers.
std::string describeOptions(const Options& opts);

}  // namespace radcor
