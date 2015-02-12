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
#ifndef __NLO_HEP_BOUNDED_VECTOR_H__
#define __NLO_HEP_BOUNDED_VECTOR_H__ 1

//   Standard includes
#include <vector>


namespace nlo {
  
  
  template<typename _Tp, class _Alloc = std::allocator<_Tp> >
  class bounded_vector : public std::vector<_Tp, _Alloc>
  {
    typedef std::vector<_Tp, _Alloc> _Base;
    
  public:
    //   types
    typedef int index_type;
    typedef typename _Base::size_type size_type;
    typedef typename _Base::value_type value_type;
    typedef typename _Base::iterator iterator;
    typedef typename _Base::pointer pointer;
    typedef typename _Base::reference reference;
    typedef typename _Base::reverse_iterator reverse_iterator;
    typedef typename _Base::const_iterator const_iterator;
    typedef typename _Base::const_reference const_reference;
    typedef typename _Base::const_pointer const_pointer;
    typedef typename _Base::const_reverse_iterator const_reverse_iterator;
    typedef typename _Base::allocator_type allocator_type;
    
    //   constructors
    bounded_vector(const allocator_type& __a = allocator_type()) 
      : _Base(__a), _M_low(0) {}
    
    bounded_vector(index_type __low, index_type __high, 
		   const _Tp& __value = _Tp(),
		   const allocator_type& __a = allocator_type()) 
      : _Base((size_type) (__high - __low + 1), __value, __a), _M_low(__low) {}
    
    bounded_vector(const bounded_vector& __x) 
      : _Base(__x), _M_low(__x.lower()) {}
    
    template <class _InIter>
    bounded_vector(_InIter __first, _InIter __last, index_type __low,
		   const allocator_type& __a = allocator_type())
      : _Base(__first, __last, __a), _M_low(__low) {}
    
    //   assignments
    bounded_vector& operator=(const bounded_vector& __x) {
      _Base::operator=(__x);
      _M_low = __x.lower();
      return *this;
    }
    
    template <class _InIter>
    void assign(_InIter __first, _InIter __last, index_type __low) {
      _Base::assign(__first, __last);
      _M_low = __low;
    }
    
    //   element access
    reference operator[](index_type __n) { 
      return _Base::operator[]((size_type) (__n - _M_low)); 
    }
    
    const_reference operator[](index_type __n) const {
      return _Base::operator[]((size_type) (__n - _M_low));
    }

    //   resize operation
    void resize(index_type __low, index_type __up, const _Tp& __val = _Tp()) {
      _M_low = __low;
      _Base::resize(__up - __low + 1, __val);
    }
    
    //   structural informations
    index_type lower() const { return _M_low;}
    index_type upper() const { return _Base::size() + _M_low - 1;}
    
  private:
    index_type _M_low;
  };
} //  namespace nlo

#endif
