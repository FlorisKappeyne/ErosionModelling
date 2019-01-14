#pragma once
#include "MathUtilities.h"

/*
* Vec3
*   A 3d vector
*
*   A mathematical structure with both direction and magnitude.
*   Contains all data and functions that a Vector3 should have such as
*   Dot, Cross, Magnitude and others.
*/

// class template
template <typename T>
class Vec3T
{
public:

  // Constructors

  // sets all values to zero
  Vec3T()
    :
    x(T(0)),
    y(T(0)),
    z(T(0))
  {
  }
  // sets all values to equal the given value
  Vec3T(const T v)
    :
    x(v),
    y(v),
    z(v)
  {
  }
  // sets all values to equal the given values
  Vec3T(const T X, const T Y, const T Z)
    :
    x(X),
    y(Y),
    z(Z)
  {
  }


  // Helper functions

  // Get the upward vector
  static inline const Vec3T GetUp()
  {
    return Vec3T(T(0), T(1), T(0));
  }
  // Get the downward vector
  static inline const Vec3T GetDown()
  {
    return Vec3T(T(0), T(-1), T(0));
  }
  // Get the righthandvector
  static inline const Vec3T GetRight()
  {
    return Vec3T(T(1), T(0), T(0));
  }
  // Get the lefthand vector
  static inline const Vec3T GetLeft()
  {
    return Vec3T(T(-1), T(0), T(0));
  }
  // Get the forward vector
  static inline const Vec3T GetForward()
  {
    return Vec3T(T(0), T(0), T(1));
  }
  // Get the backward vector
  static inline const Vec3T GetBackward()
  {
    return Vec3T(T(0), T(0), T(-1));
  }


  // Math functions

  // returns the square magnitude of the given vector
  inline const T SqrMagnitude() const
  {
    return
      x * x +
      y * y +
      z * z;
  }
  // return the magnitude of the given vector
  inline const T Magnitude() const
  {
    return Sqrt(SqrMagnitude());
  }
  // returns the normalized version of this vector
  inline const Vec3T Normalized() const
  {
    Float invMag = T(1) / Magnitude();
    return Vec3T(
      x * invMag,
      y * invMag,
      z * invMag);
  }


  // overloaded operators

  // adds this vector with another and returns it
  inline const Vec3T operator+(const Vec3T& rhs) const
  {
    return Vec3T(x + rhs.x, y + rhs.y, z + rhs.z);
  }
  // sets this vector to this vector plus another, and returns it
  inline Vec3T& operator+=(const Vec3T& rhs)
  {
    x += rhs.x, y += rhs.y, z += rhs.z;
    return *this;
  }
  // returns this vector with a Float added to it
  inline const Vec3T operator+(const T rhs) const
  {
    return Vec3T(x + rhs, y + rhs, z + rhs);
  }
  // sets this vector to this vector with a Float added to it, and returns it
  inline Vec3T& operator+=(const T rhs)
  {
    x += rhs, y += rhs, z += rhs;
    return *this;
  }
  // negates the current vector
  inline const Vec3T operator-() const
  {
    return Vec3T(-x, -y, -z);
  }
  // substracts this vector from another and returns it
  inline const Vec3T operator-(const Vec3T& rhs) const
  {
    return Vec3T(x - rhs.x, y - rhs.y, z - rhs.z);
  }
  // sets this vector to this vector minus another vector, and returns it
  inline Vec3T& operator-=(const Vec3T& rhs)
  {
    x -= rhs.x, y -= rhs.y, z -= rhs.z;
    return *this;
  }
  // substracts this vector with a Float and returns it
  inline const Vec3T operator-(const T rhs) const
  {
    return Vec3(x - rhs, y - rhs, z - rhs);
  }
  // sets this vector to this vector minus the given Float, and returns it
  inline Vec3T& operator-=(const T rhs)
  {
    x -= rhs, y -= rhs, z -= rhs;
    return *this;
  }
  // multiplies this vector with a Float and returns it
  inline const Vec3T operator*(const T rhs) const
  {
    return Vec3T(x * rhs, y * rhs, z * rhs);
  }
  // set this vector to this vector multiplied by the given Float, and returns it
  inline Vec3T& operator*=(const T rhs)
  {
    x *= rhs, y *= rhs, z *= rhs;
    return *this;
  }
  // divides the vector by the given Float and returns it
  inline const Vec3T operator/(const T rhs) const
  {
    T inv = T(1) / rhs;
    return Vec3T(x * inv, y * inv, z * inv);
  }
  // sets this vector to this vector divided by the given Float, and returns it
  inline Vec3T& operator/=(const T rhs)
  {
    T inv = T(1) / rhs;
    x *= inv, y *= inv, z *= inv;
    return *this;
  }
  // returns true if these vectors are equal, and false if not
  inline const bool operator==(const Vec3T& rhs) const
  {
    return
      x == rhs.x &&
      y == rhs.y &&
      z == rhs.z;
  }
  // returns true if these vectors are not equal, and false if they are
  inline const bool operator!=(const Vec3T& rhs) const
  {
    return
      x != rhs.x ||
      y != rhs.y ||
      z != rhs.z;
  }
  // returns a const reference to the value of the corresponding axis
  inline const T& operator[](const uint32 i) const
  {
    return (&x)[i];
  }
  // returns a reference to the value of the corresponding axis
  inline T& operator[](const uint32 i)
  {
    return (&x)[i];
  }

public:
  // Member variables
  T x, y, z;
};


// Global operators

template<typename T>
inline const Vec3T<T> operator*(const T& lhs, const Vec3T<T>& rhs)
{
  return rhs * lhs;
}


// Global functions

// returns the dot product of the given vectors
template<typename T>
inline const T Dot(const Vec3T<T>& lhs, const Vec3T<T>& rhs)
{
  return
    lhs.x * rhs.x +
    lhs.y * rhs.y +
    lhs.z * rhs.z;
}
// returns the cross product of the given vectors
template <typename T>
inline const Vec3T<T> Cross(const Vec3T<T>& lhs, const Vec3T<T>& rhs)
{
  return Vec3T<T>(
    lhs.y * rhs.z - lhs.z * rhs.y,
    lhs.z * rhs.x - lhs.x * rhs.z,
    lhs.x * rhs.y - lhs.y * rhs.x);
}
// returns the squared distance between two vectors
template <typename T>
inline const T SqrDistance(const Vec3T<T>& lhs, const Vec3T<T>& rhs)
{
  Vec3 diff = lhs - rhs;
  return diff.SqrMagnitude();
}
// returns the distance between two vectors
template <typename T>
inline const T Distance(const Vec3T<T>& lhs, const Vec3T<T>& rhs)
{
  Vec3 diff = lhs - rhs;
  return diff.Magnitude();
}
// returns a coordinate system (frame) around the given vector
template <typename T>
inline void CoordinateSystem(const Vec3T<T>& v1, Vec3T<T>* v2, Vec3T<T>* v3)
{
  //* test
  if (Abs(v1.x) > Abs(v1.y))
  {
    *v2 = Vec3T<T>(-v1.z, T(0), v1.x) / Sqrt(v1.x * v1.x + v1.z * v1.z);
  }
  else
  {
    *v2 = Vec3T<T>(T(0), v1.z, -v1.y) / Sqrt(v1.y * v1.y + v1.z * v1.z);
  }
  *v3 = Cross(v1, *v2);
}
// returns the given vector flipped to be in the hemisphere of the given normal
template <typename T>
inline const Vec3T<T> FaceForward(const Vec3T<T>& w, const Vec3T<T>& n)
{
  return (Dot(w, n) > T(0)) ? -w : w;
}


// usefull aliases

using Vec3 = Vec3T<Float>;
using Vec3I = Vec3T<int32>;
