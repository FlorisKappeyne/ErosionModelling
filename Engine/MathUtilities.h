#pragma once
#include "Typedefs.h"
#include <math.h>

#pragma region MathFunctions
/**
* redefinition of some common math functions,
* to make life easier with my Float type
*/
inline Float Sqrt(Float v)
{
#ifdef FLOAT_IS_DOUBLE
  return sqrt(v);
#else
  return sqrtf(v);
#endif
}

inline Float Exp(Float v)
{
#ifdef FLOAT_IS_DOUBLE
  return exp(v);
#else
  return expf(v);
#endif
}
inline Float Log(Float v)
{
#ifdef FLOAT_IS_DOUBLE
  return log(v);
#else
  return logf(v);
#endif
}

inline Float Pow(Float b, Float e)
{
#ifdef FLOAT_IS_DOUBLE
  return pow(b, e);
#else
  return powf(b, e);
#endif
}

inline Float Atan(Float v)
{
#ifdef FLOAT_IS_DOUBLE
  return atan(v);
#else
  return atanf(v);
#endif
}

inline Float Atan2(Float x, Float y)
{
#ifdef FLOAT_IS_DOUBLE
  return atan2(x, y);
#else
  return atan2f(x, y);
#endif
}

inline Float Sin(Float v)
{
#ifdef FLOAT_IS_DOUBLE
  return sin(v);
#else
  return sinf(v);
#endif
}

inline Float ASin(Float v)
{
#ifdef FLOAT_IS_DOUBLE
  return asin(v); // ketjap
#else
  return asinf(v);
#endif
}

inline Float Cos(Float v)
{
#ifdef FLOAT_IS_DOUBLE
  return cos(v);
#else
  return cosf(v);
#endif
}

inline Float ACos(Float v)
{
#ifdef FLOAT_IS_DOUBLE
  return acos(v);
#else
  return acosf(v);
#endif
}
#pragma endregion cos, sin, pow, exp, log etc...

#pragma region Utilities
/**
* some useful math related utilities
*/
inline Float Lerp(Float v1, Float v2, Float t)
{
  return (kOneF - t) * v1 + t * v2;
  // v1 + t * (v2 - v1) would be faster but would be less
  // accurate when the difference between v2 and v1 are high
}
template <typename T>
inline void Swap(T& l, T& r)
{
  T temp = l;
  l = r;
  r = temp;
}
template <typename T>
inline T Clamp(T val, T low, T high)
{
  val = val > high ? high : val;
  val = val < low ? low : val;
  return val;
}
template <typename T>
inline T Min(T a, T b)
{
  return a < b ? a : b;
}
template <typename T>
inline T Max(T a, T b)
{
  return a > b ? a : b;
}
template <typename T>
inline T Abs(T v)
{
  return v < T(0) ? -v : v;
}
inline int32 Mod(int32 a, int32 b)
{
  int32 n = int32(a / b);
  a -= n * b;
  if (a < 0)
    a += b;
  return a;
}
inline Float Deg2Rad(Float deg)
{
  return (kPi / Float(180.0)) * deg;
}
inline Float Rad2Deg(Float rad)
{
  return (Float(180.0) / kPi) * rad;
}
inline Float Log2(Float x)
{
  return Log(x) * kInvLog2;
}
inline int32 Log2Int(Float x)
{
  return (int32)Log2(x);
}
inline bool IsPowerOf2(int32 v)
{
  return (v & (v - 1)) == 0;
}
inline int32 RoundUpPow2(int32 v)
{
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  return v + 1;
}
template <typename Predicate> 
inline size_t FindInterval(size_t size, const Predicate& pred)
{
  size_t first = 0, len = size;
  while (len > 0)
  {
    size_t half = len >> 1, middle = first + half;
    if (pred(middle))
    {
      first = middle + 1;
      len -= half + 1;
    }
    else
    {
      len = half;
    }
  }
  return Clamp(first - 1, (size_t)0, size - 2);
}
#pragma endregion min max, mod, log2, deg2rad etc...

#pragma region Solvers
/**
* functions that can solve simple mathematical situations
*/
inline bool SolveQuadratic(Float A, Float B, Float C, Float* t0, Float* t1)
{
  Float discrim = B * B - 4.0f * A * C;
  if (discrim < 0.0f)
  {
    return false;
  }
  Float rootDiscrim = Sqrt(discrim);
  Float q;
  if (B < 0.0f)
  {
    q = -0.5f * (B - rootDiscrim);
  }
  else
  {
    q = -0.5f * (B + rootDiscrim);
  }
  *t0 = q / A;
  *t1 = C / q;
  if (*t0 > *t1)
  {
    Swap(*t0, *t1);
  }
  return true;
}
inline bool SolveLinearSystem(const Float A[2][2], const Float B[2],
  Float& x0, Float& x1)
{
  Float det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
  if (fabs(det) < 1e-10f)
  {
    return false;
  }
  x0 = (A[1][1] * B[0] - A[0][1] * B[1]) / det;
  x1 = (A[0][0] * B[1] - A[1][0] * B[0]) / det;
  if (isnan(x0) || isnan(x1))
  {
    return false;
  }
  return true;
}
#pragma endregion quadratic formulas and linear systems
