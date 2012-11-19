//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include "box.h"
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <stdexcept>
#include "raw_math_vector.h"

namespace geometry_utils
{
  Box::Box()
  {
    raw_math_vector::zero(low);
    raw_math_vector::zero(top);
  }

  Box::Box(double domainSz)
  {
    double half = domainSz / 2.0;
    for (int i = 0; i < 3; ++i) {
      low[i] = -half;
      top[i] = half;
    }
  }

  Box::Box(const double* domainSz)
  {
    for (int i = 0; i < 3; ++i) {
      double half = domainSz[i] / 2.0;
      low[i] = -half;
      top[i] = half;
    }
  }

  Box::Box(const double* low, const double* top)
  {
    if (low[0] > top[0] || low[1] >= top[1] || low[2] >= top[2])
      throw std::logic_error("low must be bottomLeft point, while top - upper right");
    raw_math_vector::copy(this->low, low);
    raw_math_vector::copy(this->top, top);
  }

  Box::Box::Box(const Box& anotherBox)
  {
    raw_math_vector::copy(this->low, anotherBox.low);
    raw_math_vector::copy(this->top, anotherBox.top);
  }

  Box& Box::operator= (const Box& anotherBox)
  {
    if (this == &anotherBox)
      return *this;

    raw_math_vector::copy(this->low, anotherBox.low);
    raw_math_vector::copy(this->top, anotherBox.top);

    return *this;
  }

  double Box::getIthSize(size_t i) const
  {
    assert(i < 3);
    return fabs(top[i] - low[i]);
  }

  double Box::getSizeX() const
  {
    return getIthSize(0);
  }

  double Box::getSizeY() const
  {
    return getIthSize(1);
  }

  double Box::getSizeZ() const
  {
    return getIthSize(2);
  }

  bool Box::inside(double* point) const
  {
    if (point[0] >= low[0] && point[0] <= top[0]
     && point[1] >= low[1] && point[1] <= top[1]
     && point[2] >= low[2] && point[2] <= top[2])
      return true;
    return false;
  }

  void Box::getCenter(double* center) const
  {
    for (int i =0 ; i < 3; ++i)
      center[i] = 0.5 * (top[i] + low[i]);
  }

  double Box::getVolume() const
  {
    return getSizeX() * getSizeY() * getSizeZ();
  }

  void Box::shift(double* newOrigin)
  {
    raw_math_vector::add(low, newOrigin);
    raw_math_vector::add(top, newOrigin);
  }

  double Box::getMaxSize() const
  {
    return std::max(getSizeX(), std::max(getSizeY(), getSizeZ()));
  }

  bool Box::isTrivial() const
  {
    //TODO use tolerance
    return getVolume() == 0.0;
  }

  const double* Box::getLow() const
  {
    return low;
  }

  const double* Box::getTop() const
  {
    return top;
  }
}
