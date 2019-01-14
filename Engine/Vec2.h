#pragma once
#include "MathUtilities.h"

/*
* Vec2
*   a 2d Vector
*
*   A mathematical structure with both direction and magnitude.
*   Contains all data and functions that a Vector2 should have such as
*   dot, scalar multiplication and division, addition and substraction
*/

// class template
template <typename T>
class Vec2T
{
public:

  // Constructors
  Vec2T()
    :
    x(T(0)),
    y(T(0))
  {
  }
  Vec2T(const T v)
    :
    x(v),
    y(v)
  {
  }
  Vec2T(const T X, const T Y)
    :
    x(X),
    y(Y)
  {
  }


  // Helper functions
  inline Vec2T GetUp()
  {
    return Vec2T(T(0), T(1));
  }
  inline Vec2T GetDown()
  {
    return Vec2T(T(0), T(-1));
  }
  inline Vec2T GetRight()
  {
    return Vec2T(T(1), T(0));
  }
  inline Vec2T GetLeft()
  {
    return Vec2T(T(-1), T(0));
  }


  // Math functions
 
  // returns the square magnitude of the given vector
  inline T SqrMagnitude()
  {
    return 
      x * x +
      y * y;
  }
  // return the magnitude of the given vector
  inline T Magnitude()
  {
    return Sqrt(SqrMagnitude());
  }
  // returns the normalized version of this vector
  inline Vec2T Normalized()
  {
    T invMag = T(1) / Magnitude();
    return Vec2T(
      x * invMag,
      y * invMag);
  }


  // operators

  // adds this vector with another and returns it
  inline const Vec2T operator+(const Vec2T& rhs) const
  {
    return Vec2T(x + rhs.x, y + rhs.y);
  }
  // sets this vector to this vector plus another, and returns it
  inline const Vec2T& operator+=(const Vec2T& rhs)
  {
    x += rhs.x, y += rhs.y;
    return *this;
  }
  // returns this vector with a Float added to it
  inline const Vec2T operator+(const T& rhs) const
  {
    return Vec2T(x + rhs, y + rhs);
  }
  // sets this vector to this vector with a Float added to it, and returns it
  inline const Vec2T& operator+=(const T& rhs)
  {
    x += rhs, y += rhs;
    return *this;
  }
  // substracts this vector from another and returns it
  inline const Vec2T operator-(const Vec2T& rhs) const
  {
    return Vec2T(x - rhs.x, y - rhs.y);
  }
  // sets this vector to this vector minus another vector, and returns it
  inline const Vec2T& operator-=(const Vec2T& rhs)
  {
    ; x -= rhs.x, y -= rhs.y;
    return *this;
  }
  // substracts this vector with a Float and returns it
  inline const Vec2T operator-(const T& rhs) const
  {
    return Vec2T(x - rhs, y - rhs);
  }
  // sets this vector to this vector minus the given Float, and returns it
  inline const Vec2T& operator-=(const T& rhs)
  {
    x -= rhs, y -= rhs;
    return *this;
  }
  // multiplies this vector with a Float and returns it
  inline const Vec2T operator*(const T& rhs) const
  {
    return Vec2T(x * rhs, y * rhs);
  }
  // set this vector to this vector multiplied by the given Float, and returns it
  inline const Vec2T& operator*=(const T& rhs)
  {
    x *= rhs, y *= rhs;
    return *this;
  }
  // divides the vector by the given Float and returns it
  inline const Vec2T operator/(const T& rhs) const
  {
    Float inv = T(1) / rhs;
    return Vec2T(x * inv, y * inv);
  }
  // sets this vector to this vector divided by the given Float, and returns it
  inline const Vec2T& operator/=(const T& rhs)
  {
    T inv = T(1) / rhs;
    x *= inv, y *= inv;
    return *this;
  }
  // sets this vector to equal the given value, and returns it
  inline Vec2T operator=(const Vec2T& rhs)
  {
    x = rhs.x;
    y = rhs.y;
    return *this;
  }
  // returns true if these vectors are equal, and false if not
  inline const bool operator==(const Vec2T& rhs) const
  {
    return  
      x == rhs.x &&
      y == rhs.y;
  }
  // returns true if these vectors are not equal, and false if they are
  inline const bool operator!=(const Vec2T& rhs) const
  {
    return  
      x != rhs.x ||
      y != rhs.y;
  }

public:
  // Member variables
  T x, y;
};


// Global operator
template<typename T>
inline const Vec2T<T> operator*(const T& lhs, const Vec2T<T>& rhs)
{
  return rhs * lhs;
}


// Global functions

// returns the dot product of the given vectors
template<typename T>
inline const T Dot(const Vec2T<T>& lhs, const Vec2T<T>& rhs)
{
  return
    lhs.x * rhs.x +
    lhs.y * rhs.y;
}
// returns the cross product of the given vectors
template <typename T>
inline const T Cross(const Vec2T<T>& lhs, const Vec2T<T>& rhs)
{
  return lhs.x * rhs.y - lhs.y * rhs.x;
}
// returns the squared distance between two vectors
template <typename T>
inline const T SqrDistance(const Vec2T<T>& lhs, const Vec2T<T>& rhs)
{
  Vec2T<T> diff = lhs - rhs;
  return diff.SqrMagnitude();
}
// returns the distance between two vectors
template <typename T>
inline const T Distance(const Vec2T<T>& lhs, const Vec2T<T>& rhs)
{
  Vec3 diff = lhs - rhs;
  return diff.Magnitude();
}
// returns the given vector flipped to be in the hemisphere of the given normal
template <typename T>
inline const Vec2T<T> FaceForward(const Vec2T<T>& w, const Vec2T<T>& n)
{
  return (Dot(w, n) > T(0)) ? -w : w;
}


// usefull aliases
using Vec2 = Vec2T<Float>;
using Vec2I = Vec2T<int32>;
