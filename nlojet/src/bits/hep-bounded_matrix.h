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
#ifndef __HEP_BOUNDED_MATRIX_H__
#define __HEP_BOUNDED_MATRIX_H__ 1

#include <bits/hep-matrix.h>


namespace nlo {



  template<typename _Tp, class _Shape, class _Alloc = std::allocator<_Tp> >
  class _Bounded_matrix_base : public matrix<_Tp, _Shape, _Alloc>
  {
    //    private types
    typedef matrix<_Tp, _Shape, _Alloc> _Base;
    typedef _Bounded_matrix_base<_Tp, _Shape, _Alloc> _Self;
    typedef _Matrix_row<_Self> matrix_row;
    typedef _Const_Matrix_row<_Self> const_matrix_row;
    
  public:
    //   types
    typedef int index_type;
    typedef _Base matrix_type;
    
    //   inhereted types form the _Base
    typedef typename _Base::value_type value_type;
    typedef typename _Base::size_type size_type;
    typedef typename _Base::allocator_type allocator_type;
    typedef typename _Base::iterator iterator;
    typedef typename _Base::pointer pointer;
    typedef typename _Base::reference reference;
    typedef typename _Base::reverse_iterator reverse_iterator;
    typedef typename _Base::const_iterator const_iterator;
    typedef typename _Base::const_pointer const_pointer;
    typedef typename _Base::const_reference const_reference;
    typedef typename _Base::const_reverse_iterator const_reverse_iterator;

  protected:
    //   constructors
    explicit _Bounded_matrix_base(const allocator_type& __a)
      : _Base(__a), _M_rlow(0), _M_clow(0) {}
    
    _Bounded_matrix_base(index_type rl, index_type rh, index_type cl, 
			 index_type ch, const _Tp& val, const allocator_type& a) 
      : _Base((size_type) (rh-rl+1), (size_type) (ch-cl+1), val, a),
    	_M_rlow(rl), _M_clow(cl) {}

    _Bounded_matrix_base(index_type low, index_type high, const _Tp& val,
			 const allocator_type& a) 
      : _Base((size_type) (high-low+1), val, a), _M_rlow(low), _M_clow(low) {}
    
    //   resize the matrix
    void _M_resize(index_type low, index_type high, const _Tp& val) {
      _M_rlow = _M_clow = low;
      _Base::resize((size_type) (high-low+1), val);
    }

    void _M_resize(index_type rl, index_type rh, index_type cl, index_type ch,
		   const _Tp& val) 
    {
      _M_rlow = rl; _M_clow = cl;
      _Base::resize((size_type) (rh-rl+1), (size_type) (ch-cl+1), val);
    }
    
  public:
    //   element access
    reference operator() (index_type r, index_type c) { 
      return _Base::operator()((size_type) (r-_M_rlow), (size_type) (c-_M_clow));
    }
    
    const_reference operator() (index_type r, index_type c) const { 
      return _Base::operator()((size_type) (r-_M_rlow), (size_type) (c-_M_clow));
    }
    
    matrix_row operator[] (index_type r) { 
      return matrix_row(*this, r);
    } 
    
    const_matrix_row operator[] (index_type r) const { 
      return const_matrix_row(*this, r);
    }  
    
    //   structural informations
    index_type lowrow() const { return _M_rlow;}
    index_type lowcol() const { return _M_clow;}
    index_type uprow() const { return _Base::nrows() + _M_rlow - 1;}
    index_type upcol() const { return _Base::ncols() + _M_clow - 1;}
    
    
  protected:
    //  lower indices
    index_type _M_rlow;
    index_type _M_clow;
  };
  
  
  template<typename _Tp, class _Shape = rectangle, class _Alloc = std::allocator<_Tp> >
  class bounded_matrix : public _Bounded_matrix_base<_Tp, _Shape, _Alloc>
  {
    //    private types
    typedef _Bounded_matrix_base<_Tp, _Shape, _Alloc> _Base;
    
  public:
    //   inhereted types form the _Base
    typedef typename _Base::index_type index_type;
    typedef typename _Base::matrix_type matrix_type;

    typedef typename _Base::value_type value_type;
    typedef typename _Base::size_type size_type;
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
    explicit bounded_matrix(const allocator_type& a = allocator_type())
      : _Base(a) {}
    
    bounded_matrix(index_type low, index_type high, const _Tp& val = _Tp(),
		   const allocator_type& a = allocator_type()) 
      : _Base(low, high, val, a) {}

    bounded_matrix(index_type rl, index_type rh, index_type cl,
		   index_type ch, const _Tp& val = _Tp(),
		   const allocator_type& a = allocator_type()) 
      : _Base(rl, rh, cl, ch, val, a) {}
    
    //   resize the matrix
    void resize(index_type low, index_type high, const _Tp& val = _Tp()) {
      _Base::_M_resize(low, high, val);
    }
    
    void resize(index_type rl, index_type rh, index_type cl, index_type ch, 
		const _Tp& val = _Tp()) 
    {
      _Base::_M_resize(rl, rh, cl, ch, val);
    }
  };


  template<typename _Tp, class _Alloc>
  class bounded_matrix<_Tp, symmetric, _Alloc> 
    : public _Bounded_matrix_base<_Tp, symmetric, _Alloc>
  {
    //    private types
    typedef _Bounded_matrix_base<_Tp, symmetric, _Alloc> _Base;
    
  public:
    //   inhereted types form the _Base
    typedef typename _Base::index_type index_type;
    typedef typename _Base::matrix_type matrix_type;

    typedef typename _Base::value_type value_type;
    typedef typename _Base::size_type size_type;
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
    explicit bounded_matrix(const allocator_type& a = allocator_type())
      : _Base(a) {}
    
    bounded_matrix(index_type low, index_type high, const _Tp& val = _Tp(),
		   const allocator_type& a = allocator_type()) 
      : _Base(low, high, val, a) {}
    
    //   resize the matrix
    void resize(index_type low, index_type high, const _Tp& val = _Tp()) {
      _Base::_M_resize(low, high, val);
    }
  };

  template<typename _Tp, class _Alloc>
  class bounded_matrix<_Tp, diagonal, _Alloc> 
    : public _Bounded_matrix_base<_Tp, diagonal, _Alloc>
  {
    //    private types
    typedef _Bounded_matrix_base<_Tp, diagonal, _Alloc> _Base;
    
  public:
    //   inhereted types form the _Base
    typedef typename _Base::index_type index_type;
    typedef typename _Base::matrix_type matrix_type;

    typedef typename _Base::value_type value_type;
    typedef typename _Base::size_type size_type;
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
    explicit bounded_matrix(const allocator_type& a = allocator_type())
      : _Base(a) {}
    
    bounded_matrix(index_type low, index_type high, const _Tp& val = _Tp(),
		   const allocator_type& a = allocator_type()) 
      : _Base(low, high, val, a) {}
    
    //   resize the matrix
    void resize(index_type low, index_type high, const _Tp& val = _Tp()) {
      _Base::_M_resize(low, high, val);
    }
  };
  

}  //  namespace nlo

#endif
    
