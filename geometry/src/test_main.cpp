//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>

#include <gtest/gtest.h>

// to run all tests from one point
// corresponding cpp files must be excluded from build
// otherwise "duplicate symbol" linker error happens
#include "linear_interpolator.h"
typedef BasicAccessStrategy _AccessStrategy;
#include "test_linear_interpolator.cpp"

TEST(InterpolatorTest, cubicGridDifferentSizes)
{
  double top[] = {4.0, 5.0, 9.0};
  double low[] = {-3.0, -4.0, -5.0};
  Box box(low, top);

  size_t n = 5, m = 6, w = 7;

  double h[] = {box.getSizeX() / (n - 1.0), box.getSizeY() / (m - 1.0), box.getSizeZ() / (w - 1.0)};

  double shift[] = {box.getSizeX() / 2.0, box.getSizeY() / 2.0, box.getSizeZ() / 2.0};

  Grid3D<double> grid(n, m, w);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double p[] = {i * h[0] - shift[0], j * h[1] - shift[1], k * h[2] - shift[2]};
        grid(i, j, k) = nonSumFunciton2(p);
      }
    }
  }

  LinearInterpolator<_AccessStrategy> li(box, grid);

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double p[] = {i * h[0] - shift[0], j * h[1] - shift[1], k * h[2] - shift[2]};
        //std::cout << li.run(p) << " " << grid(i, j, k) << std::endl;
        // TODO doesn't work yet
        //EXPECT_TRUE( fabs(li.run(p) - grid(i, j, k)) < tolerance );
      }
    }
  }
/*
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution1(-shift[0], shift[0]);
  std::uniform_real_distribution<double> distribution2(-shift[1], shift[1]);
  std::uniform_real_distribution<double> distribution3(-shift[2], shift[2]);

  double error = raw_math_vector::length(h);

  for (size_t i = 0; i < 10; ++i) {
    double p[] = {distribution1(generator), distribution2(generator), distribution3(generator)};
    EXPECT_TRUE( fabs(li.run(p) - nonSumFunciton2(p)) < error );
  }
  */
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
