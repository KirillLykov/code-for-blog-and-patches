#include <algorithm>
#include <stdexcept>

namespace AllocationUtils
{
  template<typename T>
  T* addressof(T& r) throw()
  {
    return reinterpret_cast<T*>(&const_cast<char&>(reinterpret_cast<const volatile char&>(r)));
  }

  template<typename T>
  inline void destroy(T* p) { p->~T(); }

  template<typename _ForwardIterator>
  void destroy(_ForwardIterator first, _ForwardIterator last)
  {
    for (; first != last; ++first)
      destroy(addressof(*first));
  }

  template<typename _T1>
  inline void construct(_T1* p)
  {
    ::new(static_cast<void*>(p)) _T1;
  }

  // Don't want to use uninitialized_fill_n because
  // it requires default value and a copy-constructor
  template<typename _ForwardIterator, typename _Size>
  void uninit_fill_n(_ForwardIterator first, _Size n)
  {
    _ForwardIterator cur = first;
    try
    {
      for (; n > 0; --n, ++cur)
        construct(addressof(*cur));
    }
    catch (...)
    {
      destroy(first, cur);
      throw;
    }
  }

  /**
   * @class
   *  help creating 2D and 3D arrays such that data can be accessed by [],
   *  for instance,
   *  double*** threeDArray = myArrayMemoryHelper.allocate(n, m, w);
   *  threeDArray[i][j][k] = 5
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
      uninit_fill_n(block[0], n1 * n2);
    }

    void destroy(T** block, size_t n1, size_t n2)
    {
      AllocationUtils::destroy(block[0], block[0] + n1 * n2);
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
      uninit_fill_n(block[0][0], n1 * n2 * n3);
    }

    void destroy(T*** block, size_t n1, size_t n2, size_t n3)
    {
      AllocationUtils::destroy(block[0][0], block[0][0] + n1 * n2 * n3);
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

  class array_size_error : public std::logic_error
  {
  public:
    array_size_error()
      : std::logic_error("sizes of source and destination arrays are not eqaul")
    {
    }
  };

  template< class T, template<typename X> class allocator >
  class grid_impl
  {
  protected:
    mutable allocator<T> m_allocator;
    size_t m_linearSz;
    T* m_data;

    typedef grid_impl<T, allocator> _TGridImpl;

    grid_impl(size_t linesize)
       : m_linearSz(linesize), m_data(0)
     {
       m_data = m_allocator.allocate(m_linearSz);
       uninit_fill_n(m_data, m_linearSz);
     }

    grid_impl(const _TGridImpl& another)
      : m_linearSz(another.m_linearSz)
    {
      m_data = m_allocator.allocate(m_linearSz);
      std::uninitialized_copy(another.m_data, another.m_data + m_linearSz, this->m_data);
    }

    virtual ~grid_impl()
    {
      destroy(m_data, m_data + m_linearSz);
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
  };

  template< class T, template<typename X> class allocator = std::allocator >
  class grid2D : public grid_impl<T, allocator>
  {
    size_t m_n1, m_n2;
  public:

    typedef grid2D<T, allocator> _TGrid;
    typedef grid_impl<T, allocator> _TGridImpl;

    grid2D(size_t n1, size_t n2)
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

    const T& operator() (size_t i, size_t j) const
    {
      return _TGridImpl::m_data[i + m_n1 * j];
    }

    T& operator() (size_t i, size_t j)
    {
      return _TGridImpl::m_data[i + m_n1 * j];
    }
  };

  template< class T, template<typename X> class allocator = std::allocator >
  class grid3D : public grid_impl<T, allocator>
  {
    size_t m_n1, m_n2, m_n3;
  public:

    typedef grid3D<T, allocator> _TGrid;
    typedef grid_impl<T, allocator> _TGridImpl;

    grid3D(size_t n1, size_t n2, size_t n3)
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

    const T& operator() (size_t i, size_t j, size_t k) const
    {
      return _TGridImpl::m_data[i + m_n1 * (j + m_n2 * k)];
    }

    T& operator() (size_t i, size_t j, size_t k)
    {
      return _TGridImpl::m_data[i + m_n1 * (j + m_n2 * k)];
    }
  };

}

