//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef ARRAY_MEMORY_HELPER_H_
#define ARRAY_MEMORY_HELPER_H_

#include <memory>
#include "allocation_utils_common.h"

namespace allocation_utils
{
 /**
   * @class
   *  help creating 2D and 3D arrays such that data can be accessed by [],
   *  for instance,
   *  double*** threeDArray = myArrayMemoryHelper.allocate(n, m, w);
   *  threeDArray[i][j][k] = 5;
   */
  template< class T, template<typename X> class allocator = std::allocator >
  class ArrayMemoryHelper
  {
    allocator<T> m_alloc;
    allocator<T*> m_pAlloc;
    allocator<T**> m_ppAlloc;
  public:
    // 2D
    T** allocate(size_t n1, size_t n2)
    {
      size_t totalBytes = n1 * n2;
      T* data = m_alloc.allocate(totalBytes);
      T** block = m_pAlloc.allocate(n1);

      for (size_t i = 0, n = 0; i < n1; ++i, n += n2)
      {
        block[i] = &data[n];
      }
      return block;
    }

    void deallocate(T** block, size_t n1 = 0, size_t n2 = 0)
    {
      // usually, size parameter is ignored
      m_alloc.deallocate(block[0], n1 * n2);
      m_pAlloc.deallocate(block, n1);
    }

    void create(T** block, size_t n1, size_t n2)
    {
      allocation_utils::uninitialized_fill_n_a(block[0], n1 * n2, T(), m_alloc);
    }

    void destroy(T** block, size_t n1, size_t n2)
    {
      allocation_utils::destroy(block[0], block[0] + n1 * n2, m_alloc);
    }

    void copy(T** dest, T** source, size_t n1, size_t n2)
    {
      std::copy(source[0], source[0] + n1 * n2, dest[0]);
    }

    void uninitCopy(T** dest, T** source, size_t n1, size_t n2)
    {
      std::uninitialized_copy(source[0], source[0] + n1 * n2, dest[0]);
    }

    bool equal(T** left, T** right, size_t n1, size_t n2)
    {
      return std::equal(left[0], left[0] + n1 * n2, right[0]);
    }

    // 3d
    T*** allocate(size_t n1, size_t n2, size_t n3)
    {
      size_t totalBytes = n1 * n2 * n3;
      T* data = m_alloc.allocate(totalBytes);
      size_t twoDBytes = n1 * n2;
      T** twoDblock = m_pAlloc.allocate(twoDBytes);
      T*** three3Dblock = m_ppAlloc.allocate(n1);

      size_t n = 0;
      for (size_t i = 0; i < n1; ++i) {
        size_t m = i * n2;
        three3Dblock[i] = &twoDblock[m];
        for (size_t j = 0; j < n2; ++j, n += n3) {
          twoDblock[m + j] = &data[n];
        }
      }
      return three3Dblock;
    }

    void deallocate(T*** block, size_t n1 = 0, size_t n2 = 0, size_t n3 = 0)
    {
      // usually, size parameter is ignored
      m_alloc.deallocate(block[0][0], n1 * n2 * n3);
      m_pAlloc.deallocate(block[0], n1 * n2);
      m_ppAlloc.deallocate(block, n1);
    }

    void create(T*** block, size_t n1, size_t n2, size_t n3)
    {
      allocation_utils::uninitialized_fill_n_a(block[0][0], n1 * n2 * n3, T(), m_alloc);
    }

    void destroy(T*** block, size_t n1, size_t n2, size_t n3)
    {
      allocation_utils::destroy(block[0][0], block[0][0] + n1 * n2 * n3, m_alloc);
    }

    void copy(T*** dest, T*** source, size_t n1, size_t n2, size_t n3)
    {
      std::copy(source[0][0], source[0][0] + n1 * n2 * n3, dest[0][0]);
    }

    void uninitCopy(T*** dest, T*** source, size_t n1, size_t n2, size_t n3)
    {
      std::uninitialized_copy(source[0][0], source[0][0] + n1 * n2 * n3, dest[0][0]);
    }

    bool equal(T*** left, T*** right, size_t n1, size_t n2, size_t n3)
    {
      return std::equal(left[0][0], left[0][0] + n1 * n2 * n3, right[0][0]);
    }
  };
}

#endif /* ARRAY_MEMORY_HELPER_H_ */
