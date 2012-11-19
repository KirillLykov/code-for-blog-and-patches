//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>
#include <grid.h>
#include <iostream>
#include <random> //c++11
#include "linear_interpolator.h"

using namespace geometry_utils;
using namespace containers;

const double tolerance = 0.00000001;

namespace
{
  double nonSumFunciton(const double* p)
  {
    return p[0] + p[1] * p[1] + p[2] * p[2] * p[2];
  }

  double nonSumFunciton2(const double* p)
  {
    return 1.0 + p[0] + sin(p[1])  + log( fabs(p[2]) + 1.0);
  }
}

TEST(InterpolatorTest, cubicGrid_222)
{
  Box box(1.0);

  double h = 1.0;

  size_t n = 2, m = 2, w = 2;
  Grid3D<double> grid(n, m, w);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        grid(i, j, k) = (k == 0) ? 1.0 : -1.0;
      }
    }
  }

  LinearInterpolator<_AccessStrategy> li(box, grid);

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        EXPECT_TRUE( fabs(li.run(i * h - 0.5, j * h - 0.5, k * h - 0.5) - grid(i, j, k)) < tolerance );
      }
    }
  }

  EXPECT_TRUE( fabs(li.run(0.0, 0.0, 0.0)) < tolerance );
  EXPECT_TRUE( fabs(li.run(0.25, 0.25, 0.25) + 0.5) < tolerance );
}

TEST(InterpolatorTest, cubicGrid_333)
{
  Box box(1.0);

  double h = 0.5;

  size_t n = 3, m = 3, w = 3;
  Grid3D<double> grid(n, m, w);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double p[] = {i * h - 0.5, j * h - 0.5, k * h - 0.5};
        grid(i, j, k) = nonSumFunciton(p);
      }
    }
  }

  LinearInterpolator<_AccessStrategy> li(box, grid);

  EXPECT_TRUE( fabs(li.run(0.0, 0.0, 0.0)) < tolerance );

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        EXPECT_TRUE( fabs(li.run(i * h - 0.5, j * h - 0.5, k * h - 0.5) - grid(i, j, k)) < tolerance );
      }
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-0.5, 0.5);

  for (size_t i = 0; i < 10; ++i) {
    double p[] = {distribution(generator), distribution(generator), distribution(generator)};
    EXPECT_TRUE( fabs(li.run(p) - nonSumFunciton(p)) < h );
  }
}

TEST(InterpolatorTest, cubicGrid_567)
{
  Box box(1.0);

  size_t n = 5, m = 6, w = 7;

  double h[] = {1.0 / (n - 1.0), 1.0 / (m - 1.0), 1.0 / (w - 1.0)};

  Grid3D<double> grid(n, m, w);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double p[] = {i * h[0] - 0.5, j * h[1] - 0.5, k * h[2] - 0.5};
        grid(i, j, k) = nonSumFunciton2(p);
      }
    }
  }

  LinearInterpolator<_AccessStrategy> li(box, grid);

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double p[] = {i * h[0] - 0.5, j * h[1] - 0.5, k * h[2] - 0.5};
        EXPECT_TRUE( fabs(li.run(p) - grid(i, j, k)) < tolerance );
      }
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-0.5, 0.5);

  double error = raw_math_vector::length(h);

  for (size_t i = 0; i < 10; ++i) {
    double p[] = {distribution(generator), distribution(generator), distribution(generator)};
    EXPECT_TRUE( fabs(li.run(p) - nonSumFunciton2(p)) < error );
  }
}

