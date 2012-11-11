#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>

#include <gtest/gtest.h>

#include "grid.h"
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

TEST(TwoDGridTest, construction)
{
	grid2D<double> grid(1, 2);
}

TEST(TwoDGridTest, copyConstruction)
{
//TODO check a valgrind error:
//Conditional jump or move depends on uninitialised value in equal operator
	grid2D<double> grid1(1, 2);
	grid2D<double> grid2(grid1);
	EXPECT_EQ(grid1, grid2);
}

TEST(TwoDGridTest, assignment)
{
	grid2D<double> grid1(1, 2);
	grid2D<double> grid2(1, 2);
	grid2 = grid1;
	EXPECT_EQ(grid1, grid2);
}

TEST(TwoDGridTest, assignmentBrake)
{
  grid2D<double> grid1(1, 2);
  grid2D<double> grid2(2, 3);
  try {
    grid2 = grid1;
  }
  catch (const AllocationUtils::array_size_error& er)
  {
    return;
  }
  EXPECT_TRUE(false);
}

TEST(TwoDGridTest, access)
{
  size_t n = 2, m = 4;
  grid2D<double> grid(n, m);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      grid(i, j) = i * j;
    }
  }

  const grid2D<double> grid2 = grid;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      EXPECT_EQ(grid2(i, j), i * j);
    }
  }
}

TEST(TwoDGridTest, swap)
{
  size_t n1 = 3, m1 = 2;
  grid2D<double> grid1(n1, m1);
  for (size_t i = 0; i < n1; ++i) {
    for (size_t j = 0; j < m1; ++j) {
      grid1(i, j) = i * j;
    }
  }

  size_t n2 = 5, m2 = 3;
  grid2D<double> grid2(n2, m2);
  for (size_t i = 0; i < n2; ++i) {
    for (size_t j = 0; j < m2; ++j) {
      grid2(i, j) = i + j;
    }
  }

  std::swap(grid1, grid2);

  for (size_t i = 0; i < n1; ++i) {
    for (size_t j = 0; j < m1; ++j) {
      EXPECT_EQ(grid2(i, j), i * j);
    }
  }

  for (size_t i = 0; i < n2; ++i) {
    for (size_t j = 0; j < m2; ++j) {
      EXPECT_EQ(grid1(i, j), i + j);
    }
  }

}

TEST(ThreeDGridTest, construction)
{
	grid3D<double> grid(1, 2, 3);
}

TEST(ThreeDGridTest, copyConstruction)
{
	grid3D<double> grid1(1, 2, 3);
	grid3D<double> grid2(grid1);
	EXPECT_EQ(grid1, grid2);
}

TEST(ThreeDGridTest, assignment)
{
	grid3D<double> grid1(2, 2, 3);
	grid3D<double> grid2(2, 2, 3);
	grid2 = grid1;
	EXPECT_EQ(grid1, grid2);
}

TEST(ThreeDGridTest, assignmentBrake)
{
	grid3D<double> grid1(1, 2, 5);
	grid3D<double> grid2(2, 3, 4);
	try {
		grid2 = grid1;
	} 
	catch (const AllocationUtils::array_size_error& er)
	{
		return;
	}
	EXPECT_TRUE(false);
}

TEST(ThreeDGridTest, access)
{
	size_t n = 3, m = 4, w = 7;
	grid3D<double> grid(n, m, w);
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			for (size_t k = 0; k < w; ++k) {
				grid(i, j, k) = i * j * k;
			}
		}
	}

	const grid3D<double> grid2 = grid;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			for (size_t k = 0; k < w; ++k) {
			  EXPECT_EQ(grid2(i, j, k), i * j * k);
			}
		}
	}
}


TEST(ThreeDGridTest, swap)
{
  size_t n1 = 3, m1 = 2, w1 = 4;
  grid3D<double> grid1(n1, m1, w1);
  for (size_t i = 0; i < n1; ++i) {
    for (size_t j = 0; j < m1; ++j) {
      for (size_t k = 0; k < w1; ++k) {
        grid1(i, j, k) = i * j * k;
      }
    }
  }

  size_t n2 = 5, m2 = 3, w2 = 5;
  grid3D<double> grid2(n2, m2, w2);
  for (size_t i = 0; i < n2; ++i) {
    for (size_t j = 0; j < m2; ++j) {
      for (size_t k = 0; k < w2; ++k) {
        grid2(i, j, k) = i + j + k;
      }
    }
  }

  std::swap(grid1, grid2);

  for (size_t i = 0; i < n1; ++i) {
    for (size_t j = 0; j < m1; ++j) {
      for (size_t k = 0; k < w1; ++k) {
        EXPECT_EQ(grid2(i, j, k), i * j * k);
      }
    }
  }

  for (size_t i = 0; i < n2; ++i) {
    for (size_t j = 0; j < m2; ++j) {
      for (size_t k = 0; k < w2; ++k) {
        EXPECT_EQ(grid1(i, j, k), i + j + k);
      }
    }
  }

}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
