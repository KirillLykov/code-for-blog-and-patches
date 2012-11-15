//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef GRID_H_
#define GRID_H_

#include "allocation_utils_common.h"
#include <stdexcept>

// Algorithm must contain specialized algorithms (20.6.12):
// addressof, uninitialized_copy

namespace AllocationUtils
{
  template< class T, class allocator >
  class grid_impl
  {
  protected:

    typedef typename allocator::pointer _pointer;
    typedef typename allocator::size_type _size_type;
    typedef grid_impl<T, allocator> _TGridImpl;

    mutable allocator m_allocator;
    _size_type m_linearSz;
    typename allocator::pointer m_data;

    grid_impl(_size_type linesize)
       : m_linearSz(linesize), m_data(0)
     {
       m_data = m_allocator.allocate(m_linearSz);
       AllocationUtils::uninitialized_fill_n_a(m_data, m_linearSz, T(), m_allocator);
     }

    grid_impl(const _TGridImpl& another)
      : m_linearSz(another.m_linearSz)
    {
      m_data = m_allocator.allocate(m_linearSz);
      std::uninitialized_copy(another.m_data, another.m_data + m_linearSz, this->m_data);
    }

    virtual ~grid_impl()
    {
      destroy(m_data, m_data + m_linearSz, m_allocator);
      m_allocator.deallocate(m_data, m_linearSz);
    }

    _TGridImpl& operator= (const _TGridImpl& another)
    {
      if (this == &another)
       return *this;

      std::copy(another.m_data, another.m_data + m_linearSz, this->m_data);

      return *this;
    }

    bool operator== (const _TGridImpl& another) const
    {
      return std::equal(m_data, m_data + m_linearSz, another.m_data);
    }

    bool operator!= (const _TGridImpl& another) const
    {
      return !this->operator== (another);
    }

    void swap(_TGridImpl& another)
    {
      std::swap(m_allocator, another.m_allocator);
      std::swap(m_data, another.m_data);
      std::swap(m_linearSz, another.m_linearSz);
    }
  };

  template< class T, class allocator = std::allocator<T> >
  class grid2D : public grid_impl<T, allocator>
  {
    typedef grid_impl<T, allocator> _TGridImpl;
    typedef typename _TGridImpl::_size_type _size_type;
    _size_type m_n1, m_n2;
  public:

    typedef grid2D<T, allocator> _TGrid;

    grid2D(_size_type n1, _size_type n2)
      : _TGridImpl(n1 * n2), m_n1(n1), m_n2(n2)
    {
    }

    grid2D(const _TGrid& another)
      : _TGridImpl(another), m_n1(another.m_n1), m_n2(another.m_n2)
    {
    }

    _TGrid& operator= (const _TGrid& another)
    {
      if (this == &another)
        return *this;

      if (m_n1 != another.m_n1 || m_n2 != another.m_n2)
        throw AllocationUtils::array_size_error();

      _TGridImpl::operator= (another);

      return *this;
    }

    bool operator== (const _TGrid& another) const
    {
      if (m_n1 != another.m_n1 || m_n2 != another.m_n2)
        throw AllocationUtils::array_size_error();

      return _TGridImpl::operator== (another);
    }

    bool operator!= (const _TGrid& another) const
    {
      return !this->operator== (another);
    }

    const T& operator() (_size_type i, _size_type j) const
    {
      return _TGridImpl::m_data[i + m_n1 * j];
    }

    T& operator() (_size_type i, _size_type j)
    {
      return _TGridImpl::m_data[i + m_n1 * j];
    }

    void swap(_TGrid& another)
    {
      _TGridImpl::swap(another);
      std::swap(m_n1, another.m_n1);
      std::swap(m_n2, another.m_n2);
    }
  };

  template< class T, class allocator = std::allocator<T> >
  class grid3D : public grid_impl<T, allocator>
  {
    typedef grid_impl<T, allocator> _TGridImpl;
    typedef typename _TGridImpl::_size_type _size_type;
    _size_type m_n1, m_n2, m_n3;
  public:

    typedef grid3D<T, allocator> _TGrid;

    grid3D(_size_type n1, _size_type n2, _size_type n3)
      : _TGridImpl(n1 * n2 * n3), m_n1(n1), m_n2(n2), m_n3(n3)
    {
    }

    grid3D(const _TGrid& another)
      :  _TGridImpl(another), m_n1(another.m_n1), m_n2(another.m_n2), m_n3(another.m_n3)
    {
    }

    _TGrid& operator= (const _TGrid& another)
    {
      if (this == &another)
        return *this;

      if (m_n1 != another.m_n1 || m_n2 != another.m_n2 || m_n3 != another.m_n3)
        throw AllocationUtils::array_size_error();

      _TGridImpl::operator= (another);

      return *this;
    }

    bool operator== (const _TGrid& another) const
    {
      if (m_n1 != another.m_n1 || m_n2 != another.m_n2 || m_n3 != another.m_n3)
        throw AllocationUtils::array_size_error();

      return _TGridImpl::operator== (another);
    }

    bool operator!= (const _TGrid& another) const
    {
      return !this->operator== (another);
    }

    const T& operator() (_size_type i, _size_type j, _size_type k) const
    {
      return _TGridImpl::m_data[i + m_n1 * (j + m_n2 * k)];
    }

    T& operator() (_size_type i, _size_type j, _size_type k)
    {
      return _TGridImpl::m_data[i + m_n1 * (j + m_n2 * k)];
    }

    void swap(_TGrid& another)
    {
      _TGridImpl::swap(another);
      std::swap(m_n1, another.m_n1);
      std::swap(m_n2, another.m_n2);
      std::swap(m_n3, another.m_n3);
    }
  };
}

namespace std {
  template<typename _Tp, typename _Alloc>
  inline void
  swap(AllocationUtils::grid2D<_Tp, _Alloc>& __x, AllocationUtils::grid2D<_Tp, _Alloc>& __y)
  { __x.swap(__y); }

  template<typename _Tp, typename _Alloc>
  inline void
  swap(AllocationUtils::grid3D<_Tp, _Alloc>& __x, AllocationUtils::grid3D<_Tp, _Alloc>& __y)
  { __x.swap(__y); }
}

#endif /* GRID_H_ */
