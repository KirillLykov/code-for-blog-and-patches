//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef BOX_H_
#define BOX_H_

#include <cstddef>

namespace geometry_utils
{
  class Box
  {
    double low[3];
    double top[3];
  public:
    Box();
    /**
     * Cube with center at origin
     */
    explicit Box(double domainSz);

    /**
     * Cuboid with center at origin
     */
    Box(const double* domainSz);

    Box(const double* low, const double* top);

    Box(const Box& anotherCube);
    Box& operator= (const Box& anotherCube);

    double getIthSize(size_t index) const;
    double getSizeX() const;
    double getSizeY() const;
    double getSizeZ() const;

    bool inside(double* point) const;
    void getCenter(double* center) const;
    double getVolume() const;
    void shift(double* newOrigin);
    double getMaxSize() const;

    /**
     * return true if box volume is 0
     */
    bool isTrivial() const;

    //unsafe functions because a user can try to modify internal data
    const double* getLow() const;
    const double* getTop() const;
  };
}

#endif /* BOX_H_ */
