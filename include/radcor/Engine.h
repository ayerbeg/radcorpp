// Core unfolding engine: translation of the computational part of nradcor.f
// (Slifer 2005; the unpolarized and polarized programs share it).
//
// Names follow the physics, not the Fortran. The equivalence with the Fortran names
// is given in the manual, Appendix "Fortran to C++ names"; the main ones:
//
//   Fortran                    C++                         meaning
//   main loop U:394-606     -> Engine::run                 iterative unfolding (Mo-Tsai IV)
//   CFACT     U:494-537     -> Engine::softFactor          soft+virtual factor of Stein A82
//   FES       U:1114-1192   -> Engine::integrandEs         E_s' tail integrand (photon before scattering)
//   FEP       U:1202-1285   -> Engine::integrandEp         E_p' tail integrand (photon after scattering)
//   BBB, CCC                -> tailEs, tailEp              the two tail integrals
//   FTCS      U:969-1022    -> Engine::bornModel           Born model at any (E_s', E_p')
//   TERP      U:823-890     -> Engine::interpSpectrum      model inside one spectrum, vs nu/E_s
//   SIM       U:1074-1101   -> simpson                     composite Simpson rule
//   phi       U:1196        -> bremsPhi                    bremsstrahlung shape phi(v)
//   SPENCE    U:383-389     -> spenceLegacy / spenceExact
//
// Template parameter T = float reproduces the REAL*4 arithmetic (same operation order as
// the Fortran: do not "simplify" expressions, it breaks the byte-level parity);
// T = double is the production mode. Point arrays are 1-based (index 0 unused).
#pragma once
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "radcor/FortranMath.h"
#include "radcor/Options.h"

namespace radcor {

// One spectrum (fixed beam energy) of one configuration.
template <class T>
struct Spectrum {
  T e0 = 0;           // beam energy as given [MeV]                                   (ES0)
  T es = 0;           // beam energy at the vertex, after collisional loss [MeV]       (ES)
  T etaS = 0;         // recoil factor 1 + 2 E_s sin^2(theta/2)/M_T                   (ETA)
  T mott = 0;         // Mott cross section at E_s, no recoil factor: mottNum/E_s^2   (WMOTT)
  T tHalfTarget = 0;  // half the target thickness TX/2 [r.l.]                         (TXC)
  T tBefore = 0;      // material before the vertex t_b [r.l.]                         (TIW)
  T tAfter = 0;       // material after the vertex t_a [r.l.]                          (TFW)
  int np = 0;         // number of points
  std::vector<T> ep;         // E' at the vertex, after collisional loss [MeV]         (EP)
  std::vector<T> measRatio;  // measured xs / Mott                                     (AR)
  std::vector<T> bornRatio;  // current Born model / Mott (the iterate)                (AT)
  std::vector<T> nuFrac;     // nu/E_s = (E_s - E')/E_s, the interpolation variable    (PA)
  void resize(int n) {
    np = n;
    ep.assign(n + 1, T(0)); measRatio.assign(n + 1, T(0)); bornRatio.assign(n + 1, T(0)); nuFrac.assign(n + 1, T(0));
  }
};

// Per-point quantities of the last iteration (and the Born of the previous one).
template <class T>
struct PointResult {
  T es = 0, ep = 0, nu = 0;  // kinematics used [MeV]
  T tailEs = 0;              // E_s' tail integral, Stein A82 second term      (BBB)
  T tailEp = 0;              // E_p' tail integral, Stein A82 third term       (CCC)
  T softFactor = 0;          // soft+virtual factor multiplying sigma_Born     (CFACT)
  T radModel = 0;            // radiated model softFactor*Born+tails           (CRS)
  T xs_in = 0;               // measured xs                                    (WMOTT*AR)
  T xs_born = 0;             // unfolded Born                                  (WMOTT*AT)
  T xs_born_prev = 0;        // Born at iteration NITT-1
};

// The unfolding engine for one configuration: holds all spectra measured at one angle,
// the current Born model, and the per-point results. Typical use:
//   Engine<double> e(Z, A, thetaDeg, opt);
//   e.addSpectrum(...);   // once per beam energy, in increasing E_s
//   e.run(nitt);          // iterate
//   e.results();          // per-point terms of the last iteration
template <class T>
class Engine {
 public:
  using F = std::function<T(T)>;  // integrand type passed to simpson()

  // Constants exactly as in the Fortran (U:139-145).
  const T PI = T(3.1415926536);
  const T ALPHA = T(1) / T(137.035989);
  const T ME = T(0.510999);      // electron mass [MeV]
  const T HBARC = T(197.32705);  // [MeV fm]

  // Target charge Z and mass number A, scattering angle [deg], run options. Precomputes the
  // angle-only quantities (sin^2(theta/2), Mott numerator, Spence term) as the Fortran does (U:376-389).
  Engine(T Z, T A, T thetaDeg, const Options& o) : opt_(o), Z_(Z), A_(A), thetaDeg_(thetaDeg) {
    using fm::sq;
    mTarget_ = A * T(o.amass);
    thetaRad_ = thetaDeg_ * PI / T(180);
    sin2Half_ = sq(std::sin(thetaRad_ / T(2)));
    T x = HBARC * ALPHA * std::cos(thetaRad_ / T(2)) / T(2) / sin2Half_;
    mottNum_ = sq(x) * T(1.e7);
    spence_ = o.spence_legacy ? spenceLegacy(thetaDeg_) : spenceExact(thetaDeg_);
    simpsonStep_ = T(o.simpson_step);
  }

  // --- accessors ---------------------------------------------------------------
  std::vector<Spectrum<T>>& spectra() { return sp_; }              // spectra in input order (0-based)
  const std::vector<Spectrum<T>>& spectra() const { return sp_; }
  T targetMass() const { return mTarget_; }                        // M_T = A*amass [MeV]
  T sin2Half() const { return sin2Half_; }                         // sin^2(theta/2)
  T mottNumerator() const { return mottNum_; }                     // sigma_Mott * E_s^2
  T spence() const { return spence_; }                             // Spence term used in Fbar
  T thetaDeg() const { return thetaDeg_; }                         // scattering angle [deg]
  // Per-point results of the last iteration: results()[i][j], spectrum i (0-based), point j (1-based).
  const std::vector<std::vector<PointResult<T>>>& results() const { return res_; }

  // Add a spectrum: kinematics from (e0, collisional losses) and points (nu, xs); returns its
  // 0-based index. Mirrors U:296-371 without the hole filling (done by the caller).
  int addSpectrum(T e0, T lossS, T lossP, const std::vector<T>& nu, const std::vector<T>& xs,
                  T tHalfTarget, T tBefore, T tAfter) {
    Spectrum<T> s;
    s.e0 = e0;
    s.es = e0 - lossS;
    s.mott = mottNum_ / fm::sq(s.es);
    s.etaS = T(1) + T(2) * s.es * sin2Half_ / mTarget_;
    s.tHalfTarget = tHalfTarget; s.tBefore = tBefore; s.tAfter = tAfter;
    int n = static_cast<int>(nu.size());
    s.resize(n);
    for (int j = 1; j <= n; ++j) {
      T epNominal = e0 - nu[j - 1];
      s.ep[j] = epNominal + lossP;
      s.measRatio[j] = xs[j - 1] / s.mott;
      s.bornRatio[j] = s.measRatio[j];
      s.nuFrac[j] = (s.es - s.ep[j]) / s.es;
    }
    sp_.push_back(std::move(s));
    return static_cast<int>(sp_.size()) - 1;
  }

  // Iterative unfolding (U:394-606), Gauss-Seidel order: each point is updated as soon as
  // it is computed. In RADIATE mode the Born input is not updated (U:549).
  void run(int nitt) {
    const int ne = static_cast<int>(sp_.size());
    if (ne == 0) throw std::runtime_error("Engine::run: no spectra");
    res_.assign(ne, {});
    for (int i = 0; i < ne; ++i) res_[i].resize(sp_[i].np + 1);
    const T delta = T(opt_.delta);  // soft-photon cut Delta (CUTE)
    for (int pass = 1; pass <= nitt; ++pass) {
      for (int i = 1; i <= ne; ++i) {
        Spectrum<T>& s = sp_[i - 1];
        es_ = s.es;
        const T etaS = s.etaS;
        for (int j = 1; j <= s.np; ++j) {
          ep_ = s.ep[j];
          const T etaP = T(1) / (T(1) - T(2) * ep_ * sin2Half_ / mTarget_);
          // tb_, ta_ play the role of the Fortran COMMON variables read by FES/FEP. The
          // Fortran sets them after the integrals of the current point (U:497-503), so the
          // integrals of point j use the previous point's thicknesses (B-9).
          if (!opt_.tbta_lag) setThickness(s);
          T tailEs = 0, tailEp = 0;
          // E_s' integral: from the elastic line E_s,min(E_p) to E_s - R*Delta
          T upper = es_ - etaS * etaP * delta;
          T lower = ep_ * etaP;
          if (upper <= T(0)) throw std::runtime_error("E_s' integral: upper limit <= 0");
          if (lower <= T(0)) throw std::runtime_error("E_s' integral: lower limit <= 0");
          if (!(lower - upper >= T(0)))
            tailEs = simpson(lower, upper, simpsonStep_, opt_.simpson_min, [this](T v) { return integrandEs(v); });
          // E_p' integral: from E_p + Delta to the elastic line E_p,max(E_s)
          upper = es_ / etaS;
          if (opt_.ejsz_cap) {
            T capEp = es_ / etaS - T(3);  // Whitney's 3 MeV cap (EJSZ)
            if (upper > capEp) upper = capEp;
          }
          lower = ep_ + delta;
          if (!(lower - upper >= T(0)))
            tailEp = simpson(lower, upper, simpsonStep_, opt_.simpson_min, [this](T v) { return integrandEp(v); });

          if (opt_.del_overwrite) simpsonStep_ = delta;  // B-2 (U:494: DEL = CUTE)
          setThickness(s);
          const T soft = softFactor(delta);

          const T nuFrac = (es_ - ep_) / es_;
          const T bornHere = interpSpectrum(i, nuFrac);
          const T radModel = soft * s.mott * bornHere + tailEs + tailEp;
          const bool freeze = opt_.radiate && (i < ne || !opt_.radiate_overwrite);
          if (!freeze) s.bornRatio[j] = (s.measRatio[j] - (tailEs + tailEp) / s.mott) / soft;

          PointResult<T>& r = res_[i - 1][j];
          const T born = s.mott * s.bornRatio[j];
          if (pass == nitt - 1) r.xs_born_prev = born;
          if (pass == nitt) {
            r.es = es_; r.ep = ep_; r.nu = es_ - ep_;
            r.tailEs = tailEs; r.tailEp = tailEp; r.softFactor = soft; r.radModel = radModel;
            r.xs_in = s.mott * s.measRatio[j]; r.xs_born = born;
          }
        }
        if (progress) progress(es_, pass);
      }
    }
  }

  std::function<void(T, int)> progress;  // called after each spectrum ("E_s pass", as the Fortran)

  // Radiation lengths before (t_b) and after (t_a) the vertex for spectrum s: half the target
  // plus the material upstream / downstream (zero if external RC is off).
  void setThickness(const Spectrum<T>& s) {
    if (opt_.external) { tb_ = s.tHalfTarget + s.tBefore; ta_ = s.tHalfTarget + s.tAfter; }
    else { tb_ = T(0); ta_ = T(0); }
  }

  // ---------------------------------------------------------------------------
  // Soft+virtual factor (CFACT, U:494-537): first term of Stein A82 with F-bar of A44,
  //   Fbar (R Delta/E_s)^{b t_b'} (Delta/E_p)^{b t_a'} [1 - xi/Delta/(1 - b(t_b'+t_a'))],
  // t_b' = t_b + t_r, t_a' = t_a + t_r.
  T softFactor(T delta) const {
    using fm::sq;
    const T b = T(opt_.b_stein);
    const T q2 = T(4) * es_ * ep_ * sq(std::sin(thetaRad_ / T(2)));
    T tr, fbar;
    if (opt_.internal) {
      tr = (T(2) * ALPHA / PI) * (std::log(q2 / sq(ME)) - T(1));  // equivalent radiator
      tr = tr / T(2) / b;
      T fbarM1 = T(13) * (std::log(q2 / sq(ME)) - T(1)) / T(12) - T(17) / T(36);  // Fbar - 1
      fbarM1 = fbarM1 + T(0.5) * (sq(PI) / T(6) - spence_);
      fbarM1 = fbarM1 - T(1) / T(4) * sq(std::log(es_ / ep_));
      fbarM1 = fbarM1 * (T(2) * ALPHA / PI);
      fbarM1 = fbarM1 + T(0.5772) * b * (tb_ + ta_);
      fbar = T(1) + fbarM1;
    } else {
      tr = T(0);
      fbar = T(1) + T(0.5772) * b * (tb_ + ta_);
    }
    const T tbEff = tr + tb_;
    const T taEff = tr + ta_;
    const T R = recoilR();
    const T xi = ionisationXi();
    T c = fbar;
    c = c * (std::pow(R * delta / es_, b * tbEff)) * (std::pow(delta / ep_, b * taEff));
    c = c * (T(1) - xi / delta / (T(1) - b * (tbEff + taEff)));
    return c;
  }

  // E_s' tail integrand (FES, U:1114-1192): the incident electron radiates before scattering,
  // and scatters at (E_s', E_p). esP = E_s'.
  T integrandEs(T esP) const {
    using fm::sq;
    const T b = T(opt_.b_stein);
    const T xi = ionisationXi();
    const T R = recoilR();
    const T q2 = T(4) * esP * ep_ * sq(std::sin(thetaRad_ / T(2)));
    T tr, fbarM1;
    if (opt_.internal) {
      tr = T(1) / b * (ALPHA / PI) * (std::log(q2 / sq(ME)) - T(1));
      fbarM1 = T(13) * (std::log(q2 / sq(ME)) - T(1)) / T(12) - T(17) / T(36);
      fbarM1 = fbarM1 - T(1) / T(4) * sq(std::log(es_ / ep_));
      fbarM1 = fbarM1 + T(0.5) * (sq(PI) / T(6) - spence_);
      fbarM1 = fbarM1 * (T(2) * ALPHA / PI);
      fbarM1 = fbarM1 + T(0.5772) * b * (tb_ + ta_);
    } else {
      tr = T(0);
      fbarM1 = T(0.5772) * b * (tb_ + ta_);
    }
    const T tbEff = tr + tb_;
    const T taEff = tr + ta_;
    T f = (b * tbEff / (es_ - esP)) * bremsPhi((es_ - esP) / es_);  // bremsstrahlung
    f = f + xi / (T(2) * sq(es_ - esP));                              // ionisation
    f = f * bornModel(esP, ep_) * (T(1) + fbarM1);
    f = f * std::pow((es_ - esP) / (ep_ * R), b * taEff);
    f = f * std::pow((es_ - esP) / es_, b * tbEff);
    return f;
  }

  // E_p' tail integrand (FEP, U:1202-1285): the scattered electron, of energy E_p', radiates
  // after scattering at (E_s, E_p'). epP = E_p'.
  T integrandEp(T epP) const {
    using fm::sq;
    const T b = T(opt_.b_stein);
    const T xi = ionisationXi();
    const T R = recoilR();
    const T q2 = T(4) * es_ * epP * sq(std::sin(thetaRad_ / T(2)));
    T tr, fbarM1;
    if (opt_.internal) {
      tr = T(1) / b * (ALPHA / PI) * (std::log(q2 / sq(ME)) - T(1));
      fbarM1 = T(13) * (std::log(q2 / sq(ME)) - T(1)) / T(12) - T(17) / T(36);
      fbarM1 = fbarM1 - T(1) / T(4) * sq(std::log(es_ / ep_));
      fbarM1 = fbarM1 + T(0.5) * (sq(PI) / T(6) - spence_);
      fbarM1 = fbarM1 * (T(2) * ALPHA / PI);
      fbarM1 = fbarM1 + T(0.5772) * b * (tb_ + ta_);
    } else {
      tr = T(0);
      fbarM1 = T(0.5772) * b * (tb_ + ta_);
    }
    const T tbEff = tr + tb_;
    const T taEff = tr + ta_;
    T f = (b * taEff / (epP - ep_)) * bremsPhi((epP - ep_) / epP);
    f = f + xi / (T(2) * sq(epP - ep_));
    f = f * bornModel(es_, epP) * (T(1) + fbarM1);
    f = f * std::pow((epP - ep_) / epP, b * taEff);
    f = f * std::pow((epP - ep_) * R / es_, b * tbEff);
    return f;
  }

  // Shape of the complete-screening bremsstrahlung spectrum, phi(v) = 1 - v + 3/4 v^2.
  static T bremsPhi(T v) { return T(1) - v + T(3) / T(4) * fm::sq(v); }

  // Stein's R = (M_T + E_s(1-cos)) / (M_T - E_p(1-cos)) at the current point.
  T recoilR() const {
    return (mTarget_ + es_ * (T(1) - std::cos(thetaRad_))) / (mTarget_ - ep_ * (T(1) - std::cos(thetaRad_)));
  }

  // Ionisation (Landau) parameter xi [MeV] for the total thickness t_b + t_a, target Z.
  T ionisationXi() const {
    const T eta = std::log(T(1440) * std::pow(Z_, T(-2) / T(3))) / std::log(T(183) * std::pow(Z_, T(-1) / T(3)));
    T xi = (PI * ME / T(2) / ALPHA) * (ta_ + tb_);
    xi = xi / ((Z_ + eta) * std::log(T(183) * std::pow(Z_, T(-1) / T(3))));
    return xi;
  }

  // Born model at any (E_s', E_p') (FTCS, U:969-1022): zero beyond the elastic line;
  // linear interpolation in E_s between the two neighbouring spectra at fixed nu/E_s;
  // below the lowest spectrum, its value at the same nu/E_s (scaling).
  T bornModel(T esP, T epP) const {
    if (!(epP < esP / (T(1) + T(2) * esP * sin2Half_ / mTarget_))) return T(0);
    const T nuFrac = (esP - epP) / esP;
    const int ne = static_cast<int>(sp_.size());
    if (!(esP > sp_[0].es)) return interpSpectrum(1, nuFrac) * mottNum_ / esP / esP;
    if (!(sp_[ne - 1].es >= esP)) throw std::runtime_error("ES(NE) IS LESS THAN EO");
    int i = 1;
    while (sp_[i].es - esP < T(0)) ++i;
    const T lo = interpSpectrum(i, nuFrac);
    const T hi = interpSpectrum(i + 1, nuFrac);
    const T esLo = sp_[i - 1].es, esHi = sp_[i].es;
    const T h = (lo * (esHi - esP) + hi * (esP - esLo)) / (esHi - esLo);
    return h * mottNum_ / esP / esP;
  }

  // Born/Mott of spectrum i (1-based) at nu/E_s = x (TERP, U:823-890): zero below the
  // first point, constant above the last, linear on the first two intervals (with the
  // 1e-5 offset of the Fortran), 3-point parabola elsewhere. Binary search kept verbatim.
  T interpSpectrum(int i, T x) const {
    const Spectrum<T>& s = sp_[i - 1];
    const auto& X = s.nuFrac;
    const auto& Y = s.bornRatio;
    const int n = s.np;
    if (!(x >= X[1])) return T(0);
    if (!(x < X[n])) return Y[n];
    int J = n / 2, MAX = n, MIN = 0;
    for (int guard = 0;; ++guard) {
      if (guard > 10000 || J < 1 || J >= n) throw std::runtime_error("interpSpectrum: binary search failed");
      if (x < X[J]) { MAX = J; J = J - (J - MIN + 1) / 2; continue; }
      if (x <= X[J + 1]) break;
      MIN = 1;
      J = J + (MAX - J + 1) / 2;
    }
    if (J <= 2) {
      return (Y[J] * (X[J + 1] - x) + Y[J + 1] * (x - X[J])) / (X[J + 1] - X[J] + T(0.00001));
    }
    const T X1 = X[J - 1], Y1 = Y[J - 1], X2 = X[J], Y2 = Y[J], X3 = X[J + 1], Y3 = Y[J + 1];
    const T dx = x - X2, dp = X3 - X2, dm = X1 - X2, a0 = Y2;
    T a2 = (Y3 - Y1) * (dp + dm) / (dp - dm) - (Y3 + Y1 - T(2) * Y2);
    a2 = a2 / dp / dm / T(2);
    const T a1 = (Y3 - Y1) / (dp - dm) - a2 * (dp + dm);
    return a0 + a1 * dx + a2 * fm::sq(dx);
  }

  // Composite Simpson rule on [a,b] (SIM, U:1074-1101): an even number of intervals, at least
  // minIntervals and with width at most maxStep.
  static T simpson(T a, T b, T maxStep, int minIntervals, const F& f) {
    int n = fm::fint((b - a) / maxStep / T(2));
    n = 2 * n;
    if (n < minIntervals) n = (minIntervals / 2) * 2;
    const T h = (b - a) / T(n);
    T sum = f(a);
    sum = sum + f(b);
    T x = a + h;
    sum = sum + f(x) * T(4);
    for (int k = 2; k < n; k += 2) {
      T xe = T(k) * h;
      xe = xe + a;
      const T xo = xe + h;
      sum = sum + T(2) * f(xe) + T(4) * f(xo);
    }
    return h * sum / T(3);
  }

  // Spence term as coded (U:383-389): equals Li2(sin^2(th/2)) + truncation error (bug B-1).
  T spenceLegacy(T thetaDeg) const {
    const T thr = thetaDeg * PI / T(180);
    const T arg = fm::sq(std::cos(thr / T(2)));
    T S = fm::sq(PI) / T(6) - std::log(arg) * std::log(T(1) - arg);
    for (int n = 1; n <= 50; ++n) S = S - fm::powi(arg, n) / fm::sq(T(n));
    return S;
  }
  // Spence term as printed in Stein A44: Phi(cos^2(th/2)) = Li2(cos^2(th/2)), computed in double.
  T spenceExact(T thetaDeg) const;

 private:
  Options opt_;
  T Z_, A_, thetaDeg_;
  T mTarget_ = 0;      // target mass M_T = A*amass [MeV]               (FMT)
  T thetaRad_ = 0;     //                                               (TH)
  T sin2Half_ = 0;     // sin^2(theta/2)                                (SINSQ)
  T mottNum_ = 0;      // Mott xs times E^2 (no recoil)                 (FMOTT)
  T spence_ = 0;       // Spence term of Fbar                           (SPENCE)
  T simpsonStep_ = 1;  // Simpson maximum step [MeV]                    (DEL)
  // state of the current point (Fortran COMMON)
  T es_ = 0, ep_ = 0;  // E_s, E_p at the vertex                        (EI, EF)
  T tb_ = 0, ta_ = 0;  // radiation lengths before / after              (Tb, Ta)
  std::vector<Spectrum<T>> sp_;
  std::vector<std::vector<PointResult<T>>> res_;
};

double dilog(double x);  // real dilogarithm Li2(x), x <= 1

template <class T>
T Engine<T>::spenceExact(T thetaDeg) const {
  const double thr = double(thetaDeg) * M_PI / 180.0;
  return T(dilog(std::pow(std::cos(thr / 2.0), 2)));
}

}  // namespace radcor
