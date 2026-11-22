// Helpers that reproduce Fortran-77 / gfortran arithmetic semantics so that the
// float instantiation of the engine can match the REAL*4 Fortran closely.
#pragma once
#include <cmath>
#include <cstdlib>

namespace radcor::fm {

// X**N with integer N, as libgcc __powisf2 / libgfortran pow_r4_i4 (square-and-multiply).
// Integer power, same multiplication order as gfortran (needed for float parity).
template <class T>
inline T powi(T x, int m) {
  unsigned n = m < 0 ? -m : m;
  T y = (n % 2) ? x : T(1);
  while (n >>= 1) {
    x = x * x;
    if (n % 2) y = y * x;
  }
  return m < 0 ? T(1) / y : y;
}

template <class T> inline T sq(T x) { return x * x; }  // X**2 (gfortran inlines as x*x)

// REAL -> INTEGER assignment: truncation toward zero.
template <class T> inline int fint(T x) { return static_cast<int>(x); }

}  // namespace radcor::fm
