#pragma once
#include <stdint.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <assert.h>

/**
* This file contains useful aliases, constants and defines that are used to
* standardize some common data types
*/

// signed and unsigned numbers with their size in bits for clarity in code
using uint8 = uint8_t;
using int8 = int8_t;

using uint16 = uint16_t;
using int16 = int16_t;

using uint32 = uint32_t;
using int32 = int32_t;

using uint64 = uint64_t;
using int64 = int64_t;

using f32 = float;
using f64 = double;

// if defined, we use 64 bit as standard. else, we use 32 bit
#define FLOAT_IS_DOUBLE

#ifdef FLOAT_IS_DOUBLE
///////////////////////////////////////////////////////////////////////////////
// defines for 64 bit
///////////////////////////////////////////////////////////////////////////////
using Float = double;
using Int = int64;
using HalfInt = int32;
using UInt = uint64;
using QF = __m256d;

// some useful SIMD functions
inline QF mm_set(Float d, Float c, Float b, Float a)
{
  return _mm256_set_pd(d, c, b, a);
}
inline QF mm_mul(QF lhs, QF rhs)
{
  return _mm256_mul_pd(lhs, rhs);
}
inline QF mm_sub(QF lhs, QF rhs)
{
  return _mm256_sub_pd(lhs, rhs);
}
inline QF mm_max(QF lhs, QF rhs)
{
  return _mm256_max_pd(lhs, rhs);
}
inline QF mm_min(QF lhs, QF rhs)
{
  return _mm256_min_pd(lhs, rhs);
}

// error allowance
static constexpr Float kEpsilon = 1e-11;
#else
///////////////////////////////////////////////////////////////////////////////
// defines for 32 bit
///////////////////////////////////////////////////////////////////////////////
using Float = float;
using Int = int32;
using HalfInt = int16;
using UInt = uint32;
using QF = __m128;

// some useful SIMD functions
inline QF mm_set(Float d, Float c, Float b, Float a)
{
  return _mm_set_ps(d, c, b, a);
}
inline QF mm_mul(QF lhs, QF rhs)
{
  return _mm_mul_ps(lhs, rhs);
}
inline QF mm_sub(QF lhs, QF rhs)
{
  return _mm_sub_ps(lhs, rhs);
}
inline QF mm_max(QF lhs, QF rhs)
{
  return _mm_max_ps(lhs, rhs);
}
inline QF mm_min(QF lhs, QF rhs)
{
  return _mm_min_ps(lhs, rhs);
}

// error allowance
static constexpr Float kEpsilon = 1e-4f;
#endif

// some common constants in enough precision
static constexpr Float kOneF = Float(1);
static constexpr Float kZeroF = Float(0);
static constexpr Float kHalfF = Float(0.5);
static constexpr Float kTwoF = Float(2);

static constexpr Float OneMinusEpsilon = kOneF - kEpsilon;

static constexpr Float kPi = Float(3.1415926535897932384626433832795);
static constexpr Float kTwoPi = Float(6.283185307179586476925286766559);
static constexpr Float kTwoOverPi = Float(1.5707963267948966192313216916398);
static constexpr Float kFourOverPi = Float(0.78539816339744830961);
static constexpr Float kInvPi = Float(0.31830988618379067153776752674503);
static constexpr Float kInvTwoPi = Float(0.15915494309189533576888376337251);
static constexpr Float kInvFourPi = Float(0.07957747154594766788444188168626);

static constexpr Float kInvLog2 = Float(3.3219280948873623478703194294894);

// quick and dirty keycodes
enum class KeyCode
{
  zero = 0x30,
  one = 0x31,
  two = 0x32,
  three = 0x33,
  four = 0x34,
  five = 0x35,
  six = 0x36,
  seven = 0x37,
  eight = 0x38,
  nine = 0x39,
  A = 0x41,
  B = 0x42,
  C = 0x43,
  D = 0x44,
  E = 0x45,
  F = 0x46,
  G = 0x47,
  H = 0x48,
  I = 0x49,
  J = 0x4A,
  K = 0x4B,
  L = 0x4C,
  M = 0x4D,
  N = 0x4E,
  O = 0x4F,
  P = 0x50,
  Q = 0x51,
  R = 0x52,
  S = 0x53,
  T = 0x54,
  U = 0x55,
  V = 0x56,
  W = 0x57,
  X = 0x58,
  Y = 0x59,
  Z = 0x5A
};

// a define to help make my intentions clear
#define MATH_INLINE _forceinline
