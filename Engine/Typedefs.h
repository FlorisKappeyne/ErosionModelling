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
using QI = __m256i;

// some useful SIMD functions
inline QF mm_set(Float d, Float c, Float b, Float a)
{
  return _mm256_set_pd(d, c, b, a);
}
inline QI mm_set(Int d, Int c, Int b, Int a)
{
	return _mm256_set_epi64x(d, c, b, a);
}
inline QF mm_set(Float a)
{
	return _mm256_set_pd(a, a, a, a);
}
inline QI mm_set(Int a)
{
	return _mm256_set1_epi64x(a);
}
inline QF mm_mul(QF lhs, QF rhs)
{
  return _mm256_mul_pd(lhs, rhs);
}
inline QF mm_div(QF lhs, QF rhs)
{
	return _mm256_div_pd(lhs, rhs);
}
inline QF mm_sub(QF lhs, QF rhs)
{
	return _mm256_sub_pd(lhs, rhs);
}
inline QI mm_sub(QI lhs, QI rhs)
{
	return _mm256_sub_epi64(lhs, rhs);
}
inline QF mm_add(QF lhs, QF rhs)
{
	return _mm256_add_pd(lhs, rhs);
}
inline QF mm_max(QF lhs, QF rhs)
{
  return _mm256_max_pd(lhs, rhs);
}
inline QF mm_min(QF lhs, QF rhs)
{
  return _mm256_min_pd(lhs, rhs);
}
inline QF mm_andnot(QF lhs, QF rhs)
{
	return _mm256_andnot_pd(lhs, rhs);
}
inline QF mm_hadd(QF lhs, QF rhs)
{
	return _mm256_hadd_pd(lhs, rhs);
}
inline QF mm_and(QF mask, QF rhs)
{
	return _mm256_and_pd(mask, rhs);
}
inline QF mm_xor(QF lhs, QF rhs)
{
	return _mm256_xor_pd(lhs, rhs);
}
inline QF mm_blend(QF a, QF b, QF mask)
{
	return _mm256_blendv_pd(a, b, mask);
}
inline QI mm_cmpeq(QI lhs, QI rhs)
{
	return _mm256_cmpeq_epi64(lhs, rhs);
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
using QI = __m128i;

// some useful SIMD functions
inline QF mm_set(Float d, Float c, Float b, Float a)
{
  return _mm_set_ps(d, c, b, a);
}
inline QI mm_set(Int d, Int c, Int b, Int a)
{
	return _mm_set_epi32(d, c, b, a);
}
inline QF mm_set(Float a)
{
	return _mm_set_ps1(a);
}
inline QI mm_set(Int a)
{
	return _mm_set1_epi32(a);
}
inline QF mm_mul(QF lhs, QF rhs)
{
  return _mm_mul_ps(lhs, rhs);
}
inline QF mm_div(QF lhs, QF rhs)
{
	return _mm_div_ps(lhs, rhs);
}
inline QF mm_sub(QF lhs, QF rhs)
{
	return _mm_sub_ps(lhs, rhs);
}
inline QI mm_sub(QI lhs, QI rhs)
{
	return _mm_sub_epi32(lhs, rhs);
}
inline QF mm_add(QF lhs, QF rhs)
{
	return _mm_add_ps(lhs, rhs);
}
inline QF mm_max(QF lhs, QF rhs)
{
  return _mm_max_ps(lhs, rhs);
}
inline QF mm_min(QF lhs, QF rhs)
{
  return _mm_min_ps(lhs, rhs);
}
inline QF mm_andnot(QF lhs, QF rhs)
{
	return _mm_andnot_ps(lhs, rhs);
}
inline QF mm_hadd(QF lhs, QF rhs)
{
	return _mm_hadd_ps(lhs, rhs);
}
inline QF mm_and(QI mask, QF rhs)
{
	return _mm_and_ps(*(QF*)&mask, rhs);
}
inline QF mm_xor(QF lhs, QF rhs)
{
	return _mm_xor_ps(lhs, rhs);
}
inline QF mm_blend(QF a, QF b, QF mask)
{
	return _mm_blendv_ps(a, b, mask);
}
inline QI mm_cmpeq(QI lhs, QI rhs)
{
	return _mm_cmpeq_epi32(lhs, rhs);
}

// error allowance
static constexpr Float kEpsilon = 1e-4f;
#endif

// some handy dandy overloads
inline QF operator*(QF lhs, QF rhs)
{
	return mm_mul(lhs, rhs);
}
inline QF operator/(QF lhs, QF rhs)
{
	return mm_div(lhs, rhs);
}
inline QF operator+(QF lhs, QF rhs)
{
	return mm_add(lhs, rhs);
}
inline QF operator-(QF lhs, QF rhs)
{
	return mm_sub(lhs, rhs);
}
inline QI operator==(QI lhs, QI rhs)
{
	return mm_cmpeq(lhs, rhs);
}
inline QF mm_abs(QF x) {
	static const QF sign_mask = mm_set(-Float(0)); // -0.f = 1 << 31
	return mm_andnot(sign_mask, x);
}
inline Float mm_sum_of_elements(QF vec)
{
	vec = mm_hadd(vec, vec);
	vec = mm_hadd(vec, vec);
	return ((Float*)&vec)[0];
}

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
