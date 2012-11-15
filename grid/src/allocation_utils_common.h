//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef ALLOCATION_UTILS_COMMON_H_
#define ALLOCATION_UTILS_COMMON_H_

#include <algorithm>
#include <stdexcept>

// Algorithm must contain specialized algorithms (20.6.12):
// addressof, uninitialized_copy

namespace AllocationUtils
{
  template<typename _ForwardIterator, typename _Allocator>
  void destroy(_ForwardIterator first, _ForwardIterator last, _Allocator& alloc)
  {
    for (; first != last; ++first)
      alloc.destroy( std::addressof(*first) );
  }

  // Don't want to use uninitialized_fill_n because
  // it requires default value and a copy-constructor
  template<typename _ForwardIterator, typename _Size, typename _T, typename _Allocator>
  void uninitialized_fill_n_a(_ForwardIterator first, _Size n, const _T& t, _Allocator& alloc)
  {
    _ForwardIterator cur = first;
    try
    {
      for (; n > 0; --n, ++cur)
        alloc.construct( std::addressof(*cur) );
    }
    catch (...)
    {
      destroy(first, cur, alloc);
      throw;
    }
  }

  class array_size_error : public std::logic_error
  {
  public:
    array_size_error()
      : std::logic_error("sizes of source and destination arrays are not equal")
    {
    }
  };
}

#endif /* ALLOCATION_UTILS_COMMON_H_ */
