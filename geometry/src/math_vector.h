//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef VECTOR_H_
#define VECTOR_H_

#include <memory.h>
#include <math.h>
#include <vector>
#include <assert.h>
#include "raw_math_vector.h"

namespace geometry_utils
{
  template<class T, size_t Dim>
  class MathVector
  {
  public:

    MathVector();
    explicit MathVector(const T data[Dim]);
    MathVector(const MathVector& toCopy);
    MathVector& operator= (const MathVector& toCopy);
    ~MathVector();

    double getLength() const;

    void normalize();

    void setCoord(size_t index, T val);

    T getCoord(size_t index) const;

  private:
    T data_[Dim];
  };

  // Dim == 2
  template<class T>
  class MathVector<T, 2>
  {
    T m_data[2];
  public:

    MathVector()
    {
      memset(&m_data, 0, sizeof(T) * 2);
    }

    MathVector(const T& x, const T& y)
    {
      m_data[0] = x;
      m_data[1] = y;
    }

    explicit MathVector(const T data[2])
    {
      memcpy(m_data, data, sizeof(T) * 2);
    }

    MathVector(const MathVector& toCopy)
    {
      if (this == &toCopy) return;

      memcpy(m_data, toCopy.m_data, sizeof(T) * 2);
    }

    MathVector<T,2>& operator= (const MathVector<T,2>& toCopy)
    {
      if (this == &toCopy)
          return *this;

      memcpy(m_data, toCopy.m_data, sizeof(T) * 2);

      return *this;
    }

    void setCoord(size_t index, T val)
    {
      m_data[index] = val;
    }

    T getCoord(size_t index) const
    {
      return m_data[index];
    }

    double getLength() const
    {
      double result = 0;
      for (size_t i= 0; i < 2; ++i)
        result += m_data[i]*m_data[i];

      return sqrt(result);
    }

    void normalize()
    {
      (*this) /= getLength();
    }

    T getX() const
    {
      return m_data[0];
    }

    T getY() const
    {
      return m_data[1];
    }

    void setX(const T& value)
    {
      m_data[0] = value;
    }

    void setY(const T& value)
    {
      m_data[1] = value;
    }

    ~MathVector()
    {
    }
  };

  typedef MathVector<double, 2> Vector2D;

  template<class T>
  class MathVector<T, 3>
  {
    T m_data[3];
  public:

    MathVector()
    {
      raw_math_vector::zero(m_data);
    }

    explicit MathVector(const T& x)
    {
      m_data[0] = x;
      m_data[1] = x;
      m_data[2] = x;
    }

    MathVector(const T& x, const T& y, const T& z)
    {
      m_data[0] = x;
      m_data[1] = y;
      m_data[2] = z;
    }

    MathVector(const MathVector& toCopy)
    {
      raw_math_vector::copy(m_data, toCopy.m_data);
    }

    MathVector<T,3>& operator= (const MathVector<T,3>& toCopy)
    {
      if (this == &toCopy)
        return *this;

      raw_math_vector::copy(m_data, toCopy.m_data);

      return *this;
    }

    void setCoord(size_t index, T val)
    {
      m_data[index] = val;
    }

    T getCoord(size_t index) const
    {
      return m_data[index];
    }

    double getLength() const
    {
      return raw_math_vector::length(m_data);
    }

    void normalize()
    {
      raw_math_vector::normalize(m_data);
    }

    T getX() const
    {
      return m_data[0];
    }

    T getY() const
    {
      return m_data[1];
    }

    T getZ() const
    {
      return m_data[2];
    }

    void setX(const T& value)
    {
      m_data[0] = value;
    }

    void setY(const T& value)
    {
      m_data[1] = value;
    }

    void setZ(const T& value)
    {
      m_data[2] = value;
    }
  };

  typedef MathVector<double, 3> Vector3D;

  // Operations on vectors

  template<class T, size_t Dim>
  MathVector<T,Dim>::MathVector()
  {
    memset(&data_, 0, sizeof(T) * Dim);
  }

  template<class T, size_t Dim>
  MathVector<T,Dim>::MathVector(const T data[Dim])
  {
    memcpy(data_, data, sizeof(T)*Dim);
  }

  template<class T, size_t Dim>
  MathVector<T,Dim>::MathVector(const MathVector<T,Dim>& toCopy)
  {
    memcpy(data_, toCopy.data_, sizeof(T) * Dim);
  }

  template<class T, size_t Dim>
  MathVector<T,Dim>& MathVector<T,Dim>::operator= (const MathVector<T,Dim>& toCopy)
  {
    if (this == &toCopy)
        return *this;

    memcpy(data_, toCopy.data_, sizeof(T)*Dim);

    return *this;
  }

  template<class T, size_t Dim>
  void MathVector<T,Dim>::setCoord(size_t index, T val)
  {
    assert(index < Dim);
    data_[index] = val;
  }

  template<class T, size_t Dim>
  T MathVector<T,Dim>::getCoord(size_t index) const
  {
    assert(index < Dim);
    return data_[index];
  }

  template<class T, size_t Dim>
  double MathVector<T,Dim>::getLength() const
  {
    double result = 0.0;
    for (size_t i = 0; i < Dim; ++i)
        result += data_[i] * data_[i];

    return sqrt(result);
  }

  template<class T, size_t Dim>
  void MathVector<T,Dim>::normalize()
  {
    (*this) /= getLength();
  }

  inline Vector3D getCrossProduct(const Vector3D& left, const Vector3D& right)
  {
    return Vector3D(left.getCoord(1)* right.getCoord(2) - left.getCoord(2) * right.getCoord(1),
        left.getCoord(2) * right.getCoord(0) - left.getCoord(0) * right.getCoord(2),
        left.getCoord(0) * right.getCoord(1) - left.getCoord(1) * right.getCoord(0));
  }

  //Unary operators
  template<class T, size_t Dim>
  MathVector<T, Dim>& operator *= (MathVector<T, Dim>& left, const T& scalar)
  {
    for (size_t i = 0; i < Dim; ++i)
    {
        left.setCoord(i, left.getCoord(i )* scalar);
    }
    return left;
  }

  template<class T, size_t Dim>
  MathVector<T, Dim>& operator /= (MathVector<T, Dim>& left, const T& scalar)
  {
    for (size_t i = 0; i < Dim; ++i)
    {
        left.setCoord(i, left.getCoord(i)/scalar);
    }
    return left;
  }

  template<class T, size_t Dim>
  MathVector<T, Dim>& operator += (MathVector<T, Dim>& left, const MathVector<T, Dim>& right)
  {
    for (size_t i = 0; i < Dim; ++i)
    {
        left.setCoord(i, left.getCoord(i) + right.getCoord(i));
    }
    return left;
  }

  template<class T, size_t Dim>
  MathVector<T, Dim>& operator -= (MathVector<T, Dim>& left, const MathVector<T, Dim>& right)
  {
    if (&left == &right)
        return left;

    for (size_t i = 0; i < Dim; ++i)
    {
        left.setCoord(i, left.getCoord(i) - right.getCoord(i));
    }
    return left;
  }

  //Binary operators
  template<class T, size_t Dim>
  inline T operator * (const MathVector<T, Dim>& left, const MathVector<T, Dim>& right)
  {
    T rv = 0;
    for (size_t i = 0; i < Dim ; ++i)
        rv += left.getCoord(i) * right.getCoord(i);
    return rv;
  }

  template<class T, size_t Dim>
  const MathVector<T, Dim> operator * (const MathVector<T, Dim>& left, const T& scalar)
  {
    MathVector<T, Dim> result = left;
    return result *= scalar;
  }

  template<class T, size_t Dim>
  const MathVector<T, Dim> operator * (const T& scalar, const MathVector<T, Dim>& left)
  {
    MathVector<T, Dim> result = left;
    return result *= scalar;
  }

  template<class T, size_t Dim>
  const MathVector<T, Dim> operator / (const MathVector<T, Dim>& left, const T& scalar)
  {
    MathVector<T, Dim> result = left;
    return result /= scalar;
  }

  template<class T, size_t Dim>
  const MathVector<T, Dim> operator + (const MathVector<T, Dim>& left, const MathVector<T, Dim>& right)
  {
    MathVector<T, Dim> result = left;
    return result += right;
  }

  template<class T, size_t Dim>
  const MathVector<T, Dim> operator - (const MathVector<T, Dim>& left, const MathVector<T, Dim>& right)
  {
    MathVector<T, Dim> result = left;
    return result -= right;
  }

  template<class T, size_t Dim>
  const MathVector<T, Dim> operator - (const MathVector<T, Dim>& right)
  {
    MathVector<T, Dim> result = right;
    return -1.0 * result;
  }

  template<class T, size_t Dim>
  bool operator == (const MathVector<T, Dim>& left, const MathVector<T, Dim>& right)
  {
    for (size_t i = 0; i < Dim; ++i)
    if (!raw_math_vector::isEqualLinear(left.getCoord(i), right.getCoord(i)))
      return false;
    return true;
  }

  template<class T, size_t Dim>
  bool operator != (const MathVector<T, Dim>& left, const MathVector<T, Dim>& right)
  {
    return !(left == right);
  }
}

#endif 
