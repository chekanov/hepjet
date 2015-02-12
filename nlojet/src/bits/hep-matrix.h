//  Copyright (C) 2002 Zoltan Nagy
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#ifndef __HEP_MATRIX_H__
#define __HEP_MATRIX_H__ 1

#include <memory>
#include <vector>

namespace nlo {
  

  template<class _Matrix>
  class _Matrix_row 
  {
  public:
    //   types
    typedef typename _Matrix::reference reference;
    typedef typename _Matrix::index_type index_type;
    
    //   constructor
    _Matrix_row(_Matrix& __m, index_type __r) 
      : _M_m(__m), _M_r(__r) {}
    
    //   element acces
    reference operator[] (index_type __c) {
      return _M_m(_M_r, __c);
    }

  private:
    _Matrix&   _M_m;
    index_type _M_r;
  };

  template<class _Matrix>
  class _Const_Matrix_row 
  {
  public:
    //   types
    typedef typename _Matrix::const_reference const_reference;
    typedef typename _Matrix::index_type index_type;
    
    //   constructor
    _Const_Matrix_row(const _Matrix& __m, index_type __r) 
      : _M_m(__m), _M_r(__r) {}
    
    //   element acces
    const_reference operator[] (index_type __c) const {
      return _M_m(_M_r, __c);
    }
    
  private:
    const _Matrix&   _M_m;
    index_type _M_r;
  };
  

  template <typename _Tp, class _Shape, class _Alloc>
  class _Matrix_base 
  {
    //   private types
    typedef std::vector<_Tp, _Alloc> _Container;
    typedef _Matrix_base<_Tp, _Shape, _Alloc> _Self;
    typedef _Matrix_row<_Self> matrix_row;
    typedef _Const_Matrix_row<_Self> const_matrix_row;
    
  public:
    //   shape type
    typedef _Shape shape_type;

    //   inhereted types form the _Container
    typedef typename _Container::value_type value_type;
    typedef typename _Container::size_type size_type;
    typedef typename _Container::size_type index_type;
    typedef typename _Container::allocator_type allocator_type;
    typedef typename _Container::iterator iterator;
    typedef typename _Container::pointer pointer;
    typedef typename _Container::reference reference;
    typedef typename _Container::reverse_iterator reverse_iterator;
    typedef typename _Container::const_iterator const_iterator;
    typedef typename _Container::const_pointer const_pointer;
    typedef typename _Container::const_reference const_reference;
    typedef typename _Container::const_reverse_iterator const_reverse_iterator;
    
  protected:
    //   constructors
    explicit _Matrix_base(const allocator_type& __a)
      : _M_data(__a), _M_rows(0), _M_cols(0) {}   
    
    _Matrix_base(size_type n, const _Tp& val, const _Alloc& a) 
      : _M_data(shape_type::size(n, n), val, a), _M_rows(n), _M_cols(n) {}

    _Matrix_base(size_type n1, size_type n2, const _Tp& val, const _Alloc& a) 
      : _M_data(shape_type::size(n1, n2), val, a), _M_rows(n1), _M_cols(n2) {}
    
    //   resize the matrix
    void _M_resize(size_type n, const _Tp& val) {
      _M_rows = _M_cols = n;
      _M_data.resize(shape_type::size(n, n), val);
    }
    
    void _M_resize(size_type n1, size_type n2, const _Tp& val) {
      _M_rows = n1; _M_cols = n2;
      _M_data.resize(shape_type::size(n1, n2), val);
    }
    
  public:    
    //   get the allocator object
    const allocator_type& get_allocator() const { 
      return _M_data.get_allocator(); 
    }
    
    //   iterator operations
    iterator begin() { return _M_data.begin();}
    iterator end()   { return _M_data.end();}
    
    const_iterator begin() const { return _M_data.begin();}
    const_iterator end()   const { return _M_data.end();}
    
    //   reverse iterator operations
    typename _Self::reverse_iterator rbegin() { return _M_data.rbegin();}
    typename _Self::reverse_iterator rend()   { return _M_data.rend();}
    
    const_reverse_iterator rbegin() const { return _M_data.rbegin();}
    const_reverse_iterator rend()   const { return _M_data.rend();}
    
    //  number of roows and coloums
    size_type nrows() const { return _M_rows;}
    size_type ncols() const { return _M_cols;}
    
    //   other structure informations
    size_type size() const { return _M_data.size();}
    size_type max_size() const { return _M_data.max_size();}
    size_type capacity() const { return _M_data.capacity();}
    
    //   true if the matrix is empty
    bool empty() const { return _M_data.empty();}
    
    //    element access
    reference operator() (size_type __r, size_type __c) { 
      return _M_data[shape_type::index(_M_rows, _M_cols, __r, __c)];
    }
    
    const_reference operator() (size_type __r, size_type __c) const {
      return _M_data[shape_type::index(_M_rows, _M_cols, __r, __c)];
    }
    
    matrix_row operator[] (size_type __r) { 
      return matrix_row(*this, __r);
    } 
    
    const_matrix_row operator[] (size_type __r) const { 
      return const_matrix_row(*this, __r);
    }  
    
  protected:
    //    data members
    _Container _M_data;  
    size_type  _M_rows;
    size_type  _M_cols;
  };
  
  
  struct rectangle {
    static size_t size(size_t m, size_t n) {
      return m*n;
    }
    
    static size_t index(size_t, size_t n, size_t i, size_t j) {
      return n*i + j;
    }
  };
 
  struct symmetric {
    static size_t size(size_t m, size_t) {
      return (m*(m + 1))/2;
    }
    
    static size_t index(size_t, size_t, size_t r, size_t c) {
      return  r >= c ? (r*(r+1))/2 + c : (c*(c+1))/2 + r;
    }
  };
  
  struct diagonal {
    static size_t size(size_t m, size_t) {
      if(m < 2) return m;
      else return m + 1;
    }
    
    static size_t index(size_t n, size_t, size_t r, size_t c) {
      if(r == c && n > 1) return r+1;
      else return 0;
    }
  };

  template<typename _Tp>
	struct upper_triangle {
	  static size_t size(size_t n, size_t) {
		if(n < 2) return n;
		else return (n*(n-1))/2;
	  }
	  
	  static size_t index(size_t n, size_t, size_t i, size_t j) {
		if(i < j) return 1 + n*i - (i*(i+1))/2 + j;
		else return 0;
	  }
	  
	  //   zero element
	  static const _Tp zero;
	};
  
  template<typename _Tp> const _Tp upper_triangle<_Tp>::zero = _Tp();
  
  //template<> const float upper_triangle<float>::zero = 0.0f;
  //template<> const double upper_triangle<double>::zero = 0.0;
  //template<> const long double upper_triangle<long double>::zero = 0.0L;
  
  
  
//   struct lower_triangle {
//     static size_t size(size_t n, size_t) {
//       if(n < 2) return n;
//       else return (n*(n-1))/2;
//     }
    
//     static size_t index(size_t, size_t, size_t i, size_t j) {
//       if(i > j) return 1 + (i(i-1))/2 + j;
//       else return 0;
//     }
//   };
  
  template<typename _Tp, class _Shape = rectangle, class _Alloc = std::allocator<_Tp> >
  class matrix : public _Matrix_base<_Tp, _Shape, _Alloc>
  {
    //  private types
    typedef _Matrix_base<_Tp, _Shape, _Alloc> _Base;

  public:
    //   inhereted types form the _Container
    typedef typename _Base::value_type value_type;
    typedef typename _Base::size_type size_type;
    typedef typename _Base::size_type index_type;
    typedef typename _Base::allocator_type allocator_type;
    typedef typename _Base::iterator iterator;
    typedef typename _Base::pointer pointer;
    typedef typename _Base::reference reference;
    typedef typename _Base::reverse_iterator reverse_iterator;
    typedef typename _Base::const_iterator const_iterator;
    typedef typename _Base::const_pointer const_pointer;
    typedef typename _Base::const_reference const_reference;
    typedef typename _Base::const_reverse_iterator const_reverse_iterator;

    //   constructors
    explicit matrix(const allocator_type& __a = allocator_type())
      : _Matrix_base<_Tp, _Shape, _Alloc>(__a) {}

    matrix(size_type __n, const _Tp& __value = _Tp(),
	   const allocator_type& __a = allocator_type()) 
      : _Matrix_base<_Tp, _Shape, _Alloc>(__n, __value, __a) {}

    matrix(size_type __n1, size_type __n2, const _Tp& __value = _Tp(),
	   const allocator_type& __a = allocator_type()) 
      : _Matrix_base<_Tp, _Shape, _Alloc>(__n1, __n2, __value, __a) {}

    //   resize the matrix
    void resize(size_type __n, const _Tp& __value = _Tp()) {
      this -> _M_resize(__n, __value);
    }
    
    void resize(size_type __n1, size_type __n2, const _Tp& __value = _Tp()) {
      this -> _M_resize(__n1, __n2, __value);
    }
  };

  template<typename _Tp, class _Alloc>
  class matrix<_Tp, symmetric, _Alloc> 
    : public _Matrix_base<_Tp, symmetric, _Alloc>
  {
    //  private types
    typedef _Matrix_base<_Tp, symmetric, _Alloc> _Base;

  public:
    //   inhereted types form the _Container
    typedef typename _Base::value_type value_type;
    typedef typename _Base::size_type size_type;
    typedef typename _Base::size_type index_type;
    typedef typename _Base::allocator_type allocator_type;
    typedef typename _Base::iterator iterator;
    typedef typename _Base::pointer pointer;
    typedef typename _Base::reference reference;
    typedef typename _Base::reverse_iterator reverse_iterator;
    typedef typename _Base::const_iterator const_iterator;
    typedef typename _Base::const_pointer const_pointer;
    typedef typename _Base::const_reference const_reference;
    typedef typename _Base::const_reverse_iterator const_reverse_iterator;

    //   constructors
    explicit matrix(const allocator_type& __a = allocator_type())
      : _Matrix_base<_Tp, symmetric, _Alloc>(__a) {}

    matrix(size_type __n, const _Tp& __value = _Tp(),
	   const allocator_type& __a = allocator_type()) 
      : _Matrix_base<_Tp, symmetric, _Alloc>(__n, __value, __a) {}

  public:
    //   resize the matrix
    void resize(size_type __n, const _Tp& __value = _Tp()) {
      this -> _M_resize(__n, __value);
    }
  };

  template<typename _Tp, class _Alloc>
  class matrix<_Tp, diagonal, _Alloc> 
    : public _Matrix_base<_Tp, diagonal, _Alloc>
  {
    //  private types
    typedef _Matrix_base<_Tp, diagonal, _Alloc> _Base;

  public:
    //   inhereted types form the _Container
    typedef typename _Base::value_type value_type;
    typedef typename _Base::size_type size_type;
    typedef typename _Base::size_type index_type;
    typedef typename _Base::allocator_type allocator_type;
    typedef typename _Base::iterator iterator;
    typedef typename _Base::pointer pointer;
    typedef typename _Base::reference reference;
    typedef typename _Base::reverse_iterator reverse_iterator;
    typedef typename _Base::const_iterator const_iterator;
    typedef typename _Base::const_pointer const_pointer;
    typedef typename _Base::const_reference const_reference;
    typedef typename _Base::const_reverse_iterator const_reverse_iterator;

    //   constructors
    explicit matrix(const allocator_type& __a = allocator_type())
      : _Matrix_base<_Tp, diagonal, _Alloc>(__a) {}

    matrix(size_type __n, const _Tp& __value = _Tp(),
	   const allocator_type& __a = allocator_type()) 
      : _Matrix_base<_Tp, diagonal, _Alloc>(__n, __value, __a) {}

  public:
    //   resize the matrix
    void resize(size_type __n, const _Tp& __value = _Tp()) {
      this -> _M_resize(__n, __value);
    }
  };
}    //  namespace nlo   



#endif
