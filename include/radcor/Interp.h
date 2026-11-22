// Interpolation helpers used only in the output stage of nradcor.f:
//  - splineSecondDerivs / splineEval: Numerical Recipes cubic spline SPLINE / SPLINT
//    (U:1474-1543), with the Y2(I-2) bug switchable (B-3) and the end slopes clamped to
//    YP1=YPN=0 as in the Fortran call.
//  - interpPoly: S. Choi's piecewise-polynomial INTERP (U:1579-1782; polarized variant
//    P:791-994), restricted to the branches the code uses (IFLAG=15: cubic with
//    extrapolation in U; IFLAG=1: linear without extrapolation in P), written generally.
// Arrays are 1-based (index 0 = slot that the buggy spline reads out of bounds).
#pragma once
#include <cmath>
#include <stdexcept>
#include <vector>

namespace radcor {

// Second derivatives of the cubic spline through (X[1..N], Y[1..N]) with end slopes YP1, YPN
// (>= 1e30 means natural spline). bug=true reproduces the Fortran index error (B-3); y2zero is
// the value the buggy code reads at the out-of-bounds slot Y2(0). Returns Y2[0..N].
template <class T>
std::vector<T> splineSecondDerivs(const std::vector<T>& X, const std::vector<T>& Y, int N, T YP1, T YPN,
                        bool bug, T y2zero) {
  std::vector<T> Y2(N + 1, T(0)), U(N + 1, T(0));
  Y2[0] = y2zero;
  if (YP1 > T(.99e30)) { Y2[1] = 0; U[1] = 0; }
  else {
    Y2[1] = T(-0.5);
    U[1] = (T(3) / (X[2] - X[1])) * ((Y[2] - Y[1]) / (X[2] - X[1]) - YP1);
  }
  for (int I = 2; I <= N - 1; ++I) {
    const T SIG = (X[I] - X[I - 1]) / (X[I + 1] - X[I - 1]);
    const T P = SIG * Y2[bug ? I - 2 : I - 1] + T(2);
    Y2[I] = (SIG - T(1)) / P;
    U[I] = (T(6) * ((Y[I + 1] - Y[I]) / (X[I + 1] - X[I]) - (Y[I] - Y[I - 1]) / (X[I] - X[I - 1])) /
                (X[I + 1] - X[I - 1]) - SIG * U[I - 1]) / P;
  }
  T QN, UN;
  if (YPN > T(.99e30)) { QN = 0; UN = 0; }
  else {
    QN = T(0.5);
    UN = (T(3) / (X[N] - X[N - 1])) * (YPN - (Y[N] - Y[N - 1]) / (X[N] - X[N - 1]));
  }
  Y2[N] = (UN - QN * U[N - 1]) / (QN * Y2[N - 1] + T(1));
  for (int K = N - 1; K >= 1; --K) Y2[K] = Y2[K] * Y2[K + 1] + U[K];
  return Y2;
}

// Value at X of the spline defined by (XA, YA) and the second derivatives Y2A from
// splineSecondDerivs (bisection for the interval, then the cubic formula).
template <class T>
T splineEval(const std::vector<T>& XA, const std::vector<T>& YA, const std::vector<T>& Y2A, int N, T X) {
  int KLO = 1, KHI = N;
  while (KHI - KLO > 1) {
    int K = (KHI + KLO) / 2;
    if (XA[K] > X) KHI = K; else KLO = K;
  }
  const T H = XA[KHI] - XA[KLO];
  if (H == T(0)) throw std::runtime_error("BAD XA INPUT IN SPLINT");
  const T A = (XA[KHI] - X) / H;
  const T B = (X - XA[KLO]) / H;
  return A * YA[KLO] + B * YA[KHI] + ((A * A * A - A) * Y2A[KLO] + (B * B * B - B) * Y2A[KHI]) * (H * H) / T(6);
}

// Choi's INTERP (piecewise polynomial). variant=0: unpolarized code (DELTA=1e-3, relative left/right test,
// IERR=0 on extrapolation); variant=1: polarized code (DELTA=1e-5, strict X<XP(1)/X>XP(N)).
// Returns false (IERR<0) when no value is produced; then y is left unchanged, as in Fortran.
template <class T>
bool interpPoly(int N, const std::vector<T>& XP, const std::vector<T>& YP, T X, T& y, int IFLAG,
                int variant, int* ierrOut = nullptr) {
  static const int NPT_MIN[21] = {0, 2, 3, 3, 4, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  const T DELTA = variant == 0 ? T(1.0e-3) : T(1.0e-5);
  int ierr = 0;
  int JFLAG = IFLAG % 10;
  if (N < NPT_MIN[JFLAG]) {
    if (N < 2) { if (ierrOut) *ierrOut = -3; return false; }
    else if (N == 2) JFLAG = 1;
    else if (N == 3 && JFLAG > 2) JFLAG = 2;
    else if (N == 4 && JFLAG > 5) JFLAG = 5;
    else if (N >= 5 && JFLAG > 6) JFLAG = N + 1;
  }
  bool L_LEFT = false, L_RIGHT = false;
  if (variant == 0) {
    if ((X - XP[1]) / X < -DELTA) L_LEFT = true;
    if ((X - XP[N]) / X > DELTA) L_RIGHT = true;
  } else {
    if (X < XP[1]) L_LEFT = true;
    if (X > XP[N]) L_RIGHT = true;
  }
  if (IFLAG < 10) {
    if (L_LEFT) { y = YP[1]; if (ierrOut) *ierrOut = -1; return true; }
    if (L_RIGHT) { y = YP[N]; if (ierrOut) *ierrOut = -2; return true; }
  }
  for (int I = 1; I <= N; ++I) {
    if (std::abs((X - XP[I]) / X) < DELTA) { y = YP[I]; if (ierrOut) *ierrOut = 0; return true; }
  }
  int J;
  if (L_LEFT) J = 0;
  else if (L_RIGHT) J = N;
  else {
    J = -1;
    for (int I = 1; I <= N - 1; ++I) if (X > XP[I] && X < XP[I + 1]) J = I;
  }
  if (J == -1) { if (ierrOut) *ierrOut = -4; return false; }
  if (J == 0 || J == 1) { if (JFLAG == 3 || JFLAG == 4) JFLAG = 2; }
  else if (J == N - 1 || J == N) { if (JFLAG == 2 || JFLAG == 4) JFLAG = 3; }
  int K = 0, L = 0;
  switch (JFLAG) {
    case 1: K = J; L = K + 1; break;
    case 2: K = J; L = K + 2; break;
    case 3: K = J - 1; L = J + 1; break;
    case 4: case 5: K = J - 1; L = J + 2; break;
    case 6: K = J - 1; L = J + 3; break;
    case 7: K = J - 2; L = J + 3; break;
    case 8: K = J - 2; L = J + 4; break;
    case 9: K = J - 3; L = J + 4; break;
    default: throw std::runtime_error("interpPoly: unsupported JFLAG");
  }
  if (K < 1) { L = L + (1 - K); K = 1; }
  if (L > N) { K = K - (L - N); L = N; }
  // Lagrange polynomial through np consecutive points starting at XP[k0], evaluated at X.
  auto lagrange = [&](int k0, int np) {  // Lagrange through XP[k0..k0+np-1]
    T s = 0;
    for (int a = 0; a < np; ++a) {
      T c = 1;
      for (int b = 0; b < np; ++b) if (b != a) c = c * (X - XP[k0 + b]) / (XP[k0 + a] - XP[k0 + b]);
      s = s + c * YP[k0 + a];
    }
    return s;
  };
  if (JFLAG == 1) {
    const T A = (XP[K + 1] - X) / (XP[K + 1] - XP[K]);
    const T B = (X - XP[K]) / (XP[K + 1] - XP[K]);
    y = A * YP[K] + B * YP[K + 1];
  } else if (JFLAG == 2 || JFLAG == 3) {  // INTERPOL_PARA, same algebra as Fortran
    const T D1 = X - XP[K], D2 = X - XP[K + 1], D3 = X - XP[K + 2];
    const T D12 = XP[K] - XP[K + 1], D13 = XP[K] - XP[K + 2], D23 = XP[K + 1] - XP[K + 2];
    const T A = D1 * D2 / (D13 * D23), B = D2 * D3 / (D12 * D13), C = D3 * D1 / (-D23 * D12);
    y = A * YP[K + 2] + B * YP[K] + C * YP[K + 1];
  } else if (JFLAG == 5) {  // INTERPOL_CUBIC, same algebra as Fortran
    const T D1 = X - XP[K], D2 = X - XP[K + 1], D3 = X - XP[K + 2], D4 = X - XP[K + 3];
    const T D12 = XP[K] - XP[K + 1], D13 = XP[K] - XP[K + 2], D14 = XP[K] - XP[K + 3];
    const T D23 = XP[K + 1] - XP[K + 2], D24 = XP[K + 1] - XP[K + 3], D34 = XP[K + 2] - XP[K + 3];
    const T C1 = (D2 * D3 * D4) / (D12 * D13 * D14);
    const T C2 = (D1 * D3 * D4) / (-D12 * D23 * D24);
    const T C3 = (D1 * D2 * D4) / (D13 * D23 * D34);
    const T C4 = (D1 * D2 * D3) / (-D14 * D24 * D34);
    y = T(0);
    y = y + C1 * YP[K]; y = y + C2 * YP[K + 1]; y = y + C3 * YP[K + 2]; y = y + C4 * YP[K + 3];
  } else if (JFLAG == 4) {
    throw std::runtime_error("interpPoly: IFLAG=4 (averaged parabolas) not ported (unused by nradcor)");
  } else {
    y = lagrange(K, JFLAG - 1);  // INTERPOL_NDEG with NDEG=JFLAG-2
  }
  if (variant == 1) ierr = L_LEFT ? 1 : (L_RIGHT ? 2 : 0);
  if (ierrOut) *ierrOut = ierr;
  return true;
}

}  // namespace radcor
