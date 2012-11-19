//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>
#include "math_vector.h"

//TODO write more tests
TEST(VectorTest, vector3D)
{
  using namespace geometry_utils;

  Vector3D v1(0.0, 0.0, 0.0);
  Vector3D v2(1.0, 2.0, 3.0);
  Vector3D v3 = v1 + v2;
  EXPECT_EQ(v2, v3);

  v3 *= 0.0;
  EXPECT_EQ(v1, v3);

  v3 += Vector3D(1.0, 0.0, 1.0);
  EXPECT_EQ(v3 * v2, 4.0);
}

TEST(VectorTest, vector2D)
{
  using namespace geometry_utils;

  Vector2D v1(0.0, 0.0);
  Vector2D v2(1.0, 2.0);
  Vector2D v3 = v1 + v2;
  EXPECT_EQ(v2, v3);

  v3 *= 0.0;
  EXPECT_EQ(v1, v3);

  v3 += Vector2D(1.0, 0.0);
  EXPECT_EQ(v3 * v2, 1.0);
}
