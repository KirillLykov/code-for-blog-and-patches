//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>
#include "array_memory_helper.h"

using namespace AllocationUtils;

TEST(TwoDArrayTest, allocate)
{
  ArrayMemoryHelper<double> memUtils;
  size_t n = 2, m = 3;
  double** twoDArray = memUtils.allocate(n, m);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      twoDArray[i][j] = i + j;
    }
  }

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      EXPECT_EQ(twoDArray[i][j], i + j);
    }
  }

  memUtils.deallocate(twoDArray);
}

struct TestClass
{
  double val;

  TestClass() : val (0.0)
  {
  }

  ~TestClass()
  {
    val = -1.0;
  }
};

TEST(TwoDArrayTest, create)
{
  ArrayMemoryHelper<TestClass> memUtils;
  size_t n = 2, m = 3;
  TestClass** twoDArray = memUtils.allocate(n, m);
  memUtils.create(twoDArray,n, m);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      EXPECT_EQ(twoDArray[i][j].val, 0.0);
    }
  }
  memUtils.destroy(twoDArray, n, m);
  memUtils.deallocate(twoDArray);
}

TEST(TwoDArrayTest, destroy)
{
  ArrayMemoryHelper<TestClass> memUtils;
  size_t n = 2, m = 3;
  TestClass** twoDArray = memUtils.allocate(n, m);
  memUtils.create(twoDArray, n, m);
  memUtils.destroy(twoDArray, n, m);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      EXPECT_EQ(twoDArray[i][j].val, -1.0);
    }
  }
  memUtils.deallocate(twoDArray);
}

TEST(ThreeDArrayTest, allocate)
{
  ArrayMemoryHelper<double> memUtils;
  size_t n = 2, m = 3, w = 4;
  double*** threeDArray = memUtils.allocate(n, m, w);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        threeDArray[i][j][k] = (i + 1) * (j + 1) * (k + 1);
      }
    }
  }

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        EXPECT_EQ(threeDArray[i][j][k], (i + 1) * (j + 1) * (k + 1));
      }
    }
  }

  memUtils.deallocate(threeDArray);
}

TEST(ThreeDArrayTest, create)
{
  ArrayMemoryHelper<TestClass> memUtils;
  size_t n = 2, m = 3, w = 4;
  TestClass*** twoDArray = memUtils.allocate(n, m, w);
  memUtils.create(twoDArray, n, m, w);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        EXPECT_EQ(twoDArray[i][j][k].val, 0.0);
      }
    }
  }
  memUtils.destroy(twoDArray, n, m, w);
  memUtils.deallocate(twoDArray);
}

TEST(ThreeDArrayTest, destroy)
{
  ArrayMemoryHelper<TestClass> memUtils;
  size_t n = 2, m = 3, w = 4;
  TestClass*** twoDArray = memUtils.allocate(n, m, w);
  memUtils.create(twoDArray, n, m, w);
  memUtils.destroy(twoDArray, n, m, w);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        EXPECT_EQ(twoDArray[i][j][k].val, -1.0);
      }
    }
  }
  memUtils.deallocate(twoDArray);
}


