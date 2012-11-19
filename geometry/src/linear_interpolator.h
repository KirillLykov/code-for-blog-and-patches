//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef LINEAR_INTERPOLATOR_H_
#define LINEAR_INTERPOLATOR_H_

#include <cstddef>
#include <algorithm>
#include <grid.h>
#include "box.h"
#include "raw_math_vector.h"

class BasicAccessStrategy
{
protected:
  typedef containers::Grid3D<double> _Grid;
  geometry_utils::Box m_bbox; /* bounding box for the grid */
  const _Grid& m_grid;

  BasicAccessStrategy(geometry_utils::Box box, const _Grid& grid)
  : m_bbox(box), m_grid(grid)
  {
  }

  explicit BasicAccessStrategy(const _Grid& grid)
  : m_grid(grid)
  {
  }

  double getValue(size_t i, size_t j, size_t k) const
  {
    //i %= m_grid.size(0);
    //j %= m_grid.size(1);
    //k %= m_grid.size(2);
    // handle this condition to correctly compute value on the border
    if (i >= m_grid.size(0) || j >= m_grid.size(1) || k >= m_grid.size(2)) {
      return 0.0;
    }
    return m_grid(i, j, k);
  }

  /**
   * return coordinates of the left down corner
   * if the domain was divided into sub-domains, it may have sense
   */
  double getOrigin(size_t i) const
  {
    assert(i < 3);
    return 0;
  }
};

template <class AccessStrategy = BasicAccessStrategy>
class LinearInterpolator : public AccessStrategy
{
  typedef AccessStrategy _AS;
  typedef typename _AS::_Grid _Grid;

  double h[3];

  /**
   * Interpolate using Newton interpolation polynomial (4 points are used):
   * p(x, y, z) = f0 + (f(x1, y0, z0) - f0)/(x1 - x0)*(x - x0) + (f(x0, y1, z0) - f0)/(y1 - y0)*(y - y0) +
   *  (f(x0, y0, z1) - f0)/(z1 - z0)*(z - z0)
   */
  double newton_lin_interp(const double* inputPoint, const size_t* index) const
  {
    double origin = _AS::getValue(index[0], index[1], index[2]);
    double res = origin;
    res += (_AS::getValue(index[0] + 1, index[1], index[2]) - origin) * ((double)_AS::m_grid.size(0) - 1.0)
        * (inputPoint[0] - (double)index[0] * h[0]);
    res += (_AS::getValue(index[0], index[1] + 1, index[2]) - origin) * ((double)_AS::m_grid.size(1) - 1.0)
        * (inputPoint[1] - (double)index[1] * h[1]);
    res += (_AS::getValue(index[0], index[1], index[2] + 1) - origin) * ((double)_AS::m_grid.size(2) - 1.0)
        * (inputPoint[2] - (double)index[2] * h[2]);

    return res;
  }

  /**
   * Interpolate using trilinear interpolation method
   * names of variables are from http://en.wikipedia.org/wiki/Trilinear_interpolation
   */
  double trilin_interp(const double* inputPoint, const size_t* index) const
  {
    double x0[] = {index[0] * h[0], index[1] * h[1], index[2] * h[2]};
    double xd[3];
    geometry_utils::raw_math_vector::substract(xd, inputPoint, x0);
    double c[2][2];
    for (size_t i = 0; i < 2; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        c[i][j] = _AS::getValue(index[0], index[1] + i, index[2] + j) * (h[0] - xd[0]) +
            _AS::getValue(index[0] + 1, index[1] + i, index[2] + j) * xd[0];
      }
    }

    double c0 = c[0][0] * (h[1] - xd[1]) + c[1][0] * xd[1];
    double c1 = c[0][1] * (h[1] - xd[1]) + c[1][1] * xd[1];

    double res = 1.0 / h[0] / h[1] / h[2] * (c0 * (h[2] - xd[2]) + c1 * xd[2]);
    return res;
  }

public:
  LinearInterpolator(geometry_utils::Box box, const _Grid& grid)
  : _AS(box, grid)
  {
    for (size_t i = 0; i < 3; ++i)
      h[i] = 1.0 / (_AS::m_grid.size(i) - 1.0);
  }

  double run(double x, double y, double z) const
  {
    double point[3] = {x, y, z};
    return run(point);
  }

  void computeIndex(const double* relativePosition, size_t* index) const
  {
    assert(relativePosition[0] >= 0.0 && relativePosition[1] >= 0.0 && relativePosition[2] >= 0.0);
    // index_x = floor( p_x / h_x )
    for (size_t i = 0; i < 3; ++i) {
      double coef = (_AS::m_grid.size(i) - 1.0) / static_cast<double>(_AS::m_bbox.getIthSize(i));
      index[i] = static_cast<size_t>(relativePosition[i] * coef);
    }
  }

  double run(const double* point) const
  {
    // work with Cartesian with origin in left bottom point of the domain
    // thus shift the input point. Then fin index of the cell where the point is.
    // After that interpolate function value for the point.
    double relativePosition[3];
    geometry_utils::raw_math_vector::copy(relativePosition, point);
    geometry_utils::raw_math_vector::substract(relativePosition, _AS::m_bbox.getLow());

    size_t index[3];
    computeIndex(relativePosition, index);

    return trilin_interp(relativePosition, index);
  }
};

#endif /* LINEAR_INTERPOLATOR_H_ */
