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
#ifndef __NLO_NLO_BASIC_EVENT_H__
#define __NLO_NLO_BASIC_EVENT_H__ 1

//   Standard includes
#include <iostream>

//   nlo includes
#include <bits/hep-bounded_vector.h>


namespace nlo {
  
  /// \brief This class stores the four-vector of the incoming 
  /// and outgoing QCD partons and non-QCD external particles. 
  ///
  /// Actually, this object is inhereted from the bounded_vector class. 
  /// This is a bounded vector with some restriction which stores lorentz 
  /// vector type objects. 
  /// This class doesn't have any public constructor. This type object can be 
  /// created only by an inhereted object. We don't want to use this class 
  /// directly because it doesn't carry any information about the type of the event.
  /// \param _Lv should be a lorentz vector type. This type should fulfill
  /// the the requirements of Lorentz vector type \ref lorentzvector .
  /// \see hadronic_event, lorentzvector
  template<class _Lv>
  class partonic_event : public bounded_vector<_Lv>
  {
    /// This is just an alias for the base type.
    typedef bounded_vector<_Lv> _Base;

  public:
    /// @name Index, size and value types
    //@{
    typedef typename _Base::index_type index_type;
    typedef typename _Base::size_type size_type;
    typedef typename _Base::value_type value_type;
    //@}

    /// @name Iterator types
    //@{
    typedef typename _Base::iterator iterator;
    typedef typename _Base::reverse_iterator reverse_iterator;
    typedef typename _Base::const_iterator const_iterator;
    typedef typename _Base::const_reverse_iterator const_reverse_iterator;
    //@}

    /// @name Pointer an reference types
    //@{
    typedef typename _Base::pointer pointer;
    typedef typename _Base::const_pointer const_pointer;
    typedef typename _Base::reference reference;
    typedef typename _Base::const_reference const_reference;
    //@}

    ///   Lorentz vector type is same like the value_type.
    typedef _Lv  lorentzvector_type;
    
  protected:
    ///   default constructor
    partonic_event() {}
    
    /// \brief Construct an object with given lower and upper indices
    partonic_event(index_type __low, index_type __high)
      : _Base(__low, __high) {}
  };
 
  ///  \brief Structure for tagging the hardonic_event type objects.
  /// 
  ///  This structure carry static information about the event. In the NLO
  ///  jet calculations the processes can be classify by number of incoming 
  ///  and outgoing identified hadrons and outgoing non QCD particles. This 
  ///  class store these informatios as static member of this class and as 
  ///  template parameters of it.
  ///  \param _Nqcd an unsigned int. This is the number of non QCD outgoing particles.
  ///  \param _Inh an unsigned int. This is the number of the incoming hadrons.
  ///  \param _Idh an unsigned int. This is the number of the outgoing idetified hadrons
  template<unsigned int _Nqcd, unsigned int _Inh, unsigned int _Idh>
  struct hadronic_event_traits {
    /// Number of non QCD outgoing particles
    static const unsigned int non_qcd_outgoings  = _Nqcd;
    /// Number of the incoming hadrons
    static const unsigned int incoming_hadrons   = _Inh;
    /// Number of the outgoing idetified hadrons
    static const unsigned int identified_hadrons = _Idh;
  };  

  /// \defgroup group1 Index classes for hadronic_event
  /// Classes for implementing the subscript operator of 
  /// class hadronic_event.
  //@{
  /// \brief This class labels the incoming and outgoing hardrons.
  /// 
  ///  Using this helper classes the acces to the hadron momenta can be 
  ///  computed in the following way:
  ///  \verbatim 
  ///  typedef lorentzvector<double> _Lv;
  ///  typedef hadronic_event_traits<0,2,3> _EventTraits;  
  ///  hadronic_event_traits<_Lv,_EventTraits> p(4);  
  ///  lorentzvector<double> Ia, Ib, H1, H2;
  ///
  ///  //  acces to the incoming hadrons
  ///  Ia = p[hadron(-1)]; 
  ///  Ib = p[hadron(hadron::B)];
  ///
  ///  //  or some example for the outgoings
  ///  H1 = p[hadron(1)];
  ///  H2 = p[hadron(hadron::H2)];
  ///  \endverbatim
  ///  The index runs on the -1,0,1,... The -1 and 0 always label the incoming 
  ///  hadrons and the positive integers label the outgoing hadrons. This class 
  ///  doesn't do any check on the bound this type is just a kind of \em tagged 
  ///  \em integer.
  ///  There are some other information about this subscript operator in 
  ///  documentation of class \ref hadronic_event. 
  ///  \see hadronic_event
  struct hadron {
    /// Aliases to the incoming (A,B) and outgoing (H1,H2,...) hadrons
    enum hadron_label { B = -1, A = 0, H1, H2, H3, H4, H5};
    
    /// Creates an object using the hadron aliases
    hadron(hadron_label hdr) : idx(hdr) {}

    /// Creates an object using the int type indices
    hadron(int hdr) : idx((int) hdr) {}
    
    /// Index of the hadron
    int idx;
  };
  
  /// \brief This class labels the outgoing non QCD particles.
  /// 
  ///  Using this helper classes the acces to the momenta of non QCD particle 
  ///  can be computed in the following simple way:
  ///  \verbatim 
  ///  typedef lorentzvector<double> _Lv;
  ///  typedef hadronic_event_traits<2,2,0> _EventTraits;  
  ///  hadronic_event_traits<_Lv,_EventTraits> p(4);  
  ///  lorentzvector<double> f1, f2;
  ///
  ///  //  acces to the incoming hadrons
  ///  f1 = p[non_qcd(1)]; 
  ///  f2 = p[non_qcd(2)];
  ///  \endverbatim 
  ///  The index runs on the positive integer numbers (1,2,3,..). This class 
  ///  doesn't do any check on the bound this type is just a kind of \em tagged 
  ///  \em integer.
  ///  There are some other information about this subscript operator in 
  ///  documentation of class \ref hadronic_event. 
  ///  \see hadronic_event
  struct non_qcd {
    /// Creates an object using the unsigned int type indices
    non_qcd(unsigned int __nqcd) : idx(__nqcd) {} 
    
    /// Index of the non QCD particle
    unsigned int idx;
  };

  /// \brief This class labels the incoming and outgoing partons.
  /// 
  ///  Using this helper classes the acces to the identidied parton 
  ///  momenta can be computed in the following way:
  ///  \verbatim 
  ///  typedef lorentzvector<double> _Lv;
  ///  typedef hadronic_event_traits<0,2,3> _EventTraits;  
  ///  hadronic_event_traits<_Lv,_EventTraits> p(4);  
  ///  lorentzvector<double> Ia, Ib, H1, H2;
  ///
  ///  //  acces to the incoming hadrons
  ///  Ia = p[hadron(-1)]; 
  ///  Ib = p[hadron(hadron::B)];
  ///
  ///  //  or some example for the outgoings
  ///  H1 = p[hadron(1)];
  ///  H2 = p[hadron(hadron::H2)];
  ///  \endverbatim
  ///  The index runs on the -1,0,1,... The -1 and 0 always label the incoming 
  ///  partons and the positive integers labels the outgoing identified partons.
  ///  This class doesn't do any check on the bound this type is just a kind of \em tagged 
  ///  \em integer.
  ///  There are some other information about this subscript operator in 
  ///  documentation of class \ref hadronic_event. 
  ///  \see hadronic_event
  struct identified {
    /// an alias to the parton
    enum in_parton_type { B = -1, A = 0, H1, H2, H3, H4, H5};
    
    /// Create an object using the parton aliases
    identified(in_parton_type __id) : idx(__id) {}

    /// Create an object using the int type indices
    identified(int __id) : idx(__id) {}
    
    /// Index of the identified parton
    int idx;
  };

  /// \brief This class labels the outgoing nonidentified partons.
  /// 
  ///  Using this helper classes the acces to the momenta of the 
  ///  non-identidied (always outgoing) parton can be computed.
  ///  Here is a simple example:
  ///  \verbatim 
  ///  typedef lorentzvector<double> _Lv;
  ///  typedef hadronic_event_traits<0,2,3> _EventTraits;  
  ///  hadronic_event_traits<_Lv,_EventTraits> p(4);  
  ///  lorentzvector<double> p1,p2,p3;
  ///
  ///  //  acces to the the incoming hadrons
  ///  p1 = p[non_identified(1)];
  ///  p2 = p[non_identified(2)];
  ///  p3 = p[non_identified(3)];
  ///  \endverbatim
  ///  The index runs on the positive integer numbers (1,2,3,...).
  ///  This class doesn't do any check on the bound this type is just a kind of \em tagged 
  ///  \em integer.
  ///  There are some other information about this subscript operator in 
  ///  documentation of class \ref hadronic_event. 
  /// \see hadronic_event
  struct non_identified {
    /// Create an object using the unsigned int type indices  
    non_identified(unsigned int prt) : idx(prt) {}
    
    /// Index of the nonidentified parton
    unsigned int idx;
  };
  // @}
  
  template<class _Lv, class _Traits>
  class hadronic_event : public partonic_event<_Lv>
  {
    //   private types
    typedef _Traits _Tr;
    typedef partonic_event<_Lv> _Base;

  public:
    //   inherited types
    typedef typename _Base::index_type index_type;
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

    //   public types
    typedef _Traits             traits_type;
    typedef partonic_event<_Lv> partonic_event_type;
    typedef typename _Base::lorentzvector_type lorentzvector_type;
    
    //    constructors
    hadronic_event() {}

    explicit hadronic_event(size_type __n) 
      : _Base(-_Tr::non_qcd_outgoings - 1, _Tr::identified_hadrons + __n),
	_M_hadron(1 - _Tr::incoming_hadrons, _Tr::identified_hadrons) {}
    
    
    //  hadron iterators
    iterator hadron_begin() { return _M_hadron.begin();}
    iterator hadron_end() { return _M_hadron.end();}
    
    reverse_iterator hadron_rbegin() { return _M_hadron.rbegin();}
    reverse_iterator hadron_rend() { return _M_hadron.rend();}
    
    const_iterator hadron_begin() const {  return _M_hadron.begin();}
    const_iterator hadron_end() const { return _M_hadron.end();}
    
    const_reverse_iterator hadron_rbegin() const { return _M_hadron.rbegin();}
    const_reverse_iterator hadron_rend() const { return _M_hadron.rend();}
        
    //       element access
    reference operator[](index_type __n) { 
      return _Base::operator[](__n); 
    }
    
    const_reference operator[](index_type __n) const {
      return _Base::operator[](__n);
    }

    reference operator[](const hadron& __hdr) {
      return _M_hadron[__hdr.idx];
    }
    
    const_reference operator[](const hadron& __hdr) const {
      return _M_hadron[__hdr.idx];
    }

    reference operator[](const non_qcd& __nqcd) {
      return *(_Base::begin() + __nqcd.idx - 1);
    }
    
    const_reference operator[](const non_qcd& __nqcd) const {
      return *(_Base::begin() + __nqcd.idx - 1);
    }
    
    reference operator[](const identified& __idp) {
      return _Base::operator[](__idp.idx);
    }
    
    const_reference operator[](const identified& __idp) const {
      return _Base::operator[](__idp.idx);
    }
    
    reference operator[](const non_identified& __par) {
      return _Base::operator[](_Tr::identified_hadrons + __par.idx);
    }
    
    const_reference operator[](const non_identified& __par) const {
      return _Base::operator[](_Tr::identified_hadrons + __par.idx);
    }

    //    member functions
    void resize(size_type __n) {
      _Base::resize(-((int) _Tr::non_qcd_outgoings) - 1,
		    _Tr::identified_hadrons + __n);
    }

    //   structural informations
    size_type non_qcd_outgoings () const { return _Tr::non_qcd_outgoings;}
    size_type incoming_hadrons  () const { return _Tr::incoming_hadrons;}
    size_type identified_hadrons() const { return _Tr::identified_hadrons;}
    size_type idendified_partons() const { return _M_hadron.size();}
    
    size_type non_idendified_partons() const {
      return _Base::upper() - _Tr::indentified_hadrons;
    }

    index_type hadron_upper() const { return _M_hadron.upper();}
    index_type hadron_lower() const { return _M_hadron.lower();}

    index_type identified_upper() const { return _M_hadron.upper();}
    index_type identified_lower() const { return _M_hadron.lower();}

  protected:
    //   incoming and outgoing hadrons
    bounded_vector<lorentzvector_type> _M_hadron;  
  };


  //   I/O operations
  template<class _Lv, class _Tr> std::ostream& 
  operator<<(std::ostream& __os, const hadronic_event<_Lv, _Tr>& __a) 
  {
    for(int __j = __a.hadron_lower(); __j <= __a.hadron_upper(); __j++)
      __os<<"H["<<__j<<"] : "<<__a[hadron(__j)]
	  <<"  "<<__a[hadron(__j)].mag2()<<std::endl;
    
    for(int __i = __a.lower(); __i <= __a.upper(); __i++)
      __os<<"P["<<__i<<"] : "<<__a[__i]<<"  "<<__a[__i].mag2()<<std::endl;
    
    return __os;
  }

  //   Copies the basic_event object to a pre-allocated array. 
  template<class _Lv, class _Tr>
  void copy_to_array(typename _Lv::value_type *__p, 
		     const hadronic_event<_Lv, _Tr>& __ev)
  {
    unsigned int __j = 0U;
    
    typename hadronic_event<_Lv, _Tr>::const_iterator __h_iter, 
      __h_low = __ev.hadron_begin(), __h_up = __ev.hadron_end();
    
    for(__h_iter = __h_low; __h_iter < __h_up; __h_iter++) {
      __p[__j]    = __h_iter -> X();
      __p[__j+1U] = __h_iter -> Y();
      __p[__j+2U] = __h_iter -> Z();
      __p[__j+3U] = __h_iter -> T();
      __j += 4U;
    }
    
    typename hadronic_event<_Lv, _Tr>::const_iterator __iter, 
      __low = __ev.begin(), __up = __ev.end();
    
    for(__iter = __low; __iter < __up; __iter++) {
      __p[__j]    = __iter -> X();
      __p[__j+1U] = __iter -> Y();
      __p[__j+2U] = __iter -> Z();
      __p[__j+3U] = __iter -> T();
      __j += 4U;
    }
  }

  //
  //    Specialization
  //
  template<class _Lv>
  class hadronic_event<_Lv, hadronic_event_traits<0U,0U,0U> > 
    : public partonic_event<_Lv>
  {
    //   private types
    typedef partonic_event<_Lv> _Base;

  public:
    //   inherited types
    typedef typename _Base::index_type index_type;
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

    //   public types
    typedef hadronic_event_traits<0U,0U,0U> traits_type;
    typedef partonic_event<_Lv> partonic_event_type;
    typedef typename _Base::lorentzvector_type lorentzvector_type;

    //    constructors
    hadronic_event() {}

    explicit hadronic_event(size_type __n) : _Base(- 1, __n) {}    
    
    //       element access
    reference operator[](index_type __n) { 
      return _Base::operator[](__n); 
    }
    
    const_reference operator[](index_type __n) const {
      return _Base::operator[](__n);
    }

    reference operator[](const non_identified& __par) {
      return _Base::operator[](__par.idx);
    }
    
    const_reference operator[](const non_identified& __par) const {
      return _Base::operator[](__par.idx);
    }

    //    member functions
    void resize(size_type __n) { _Base::resize(-1, __n);}

    //   structural informations
    size_type non_idendified_partons() const { return _Base::upper();}
  };
  
  //   I/O operations
  template<class _Lv> std::ostream& 
  operator<<(std::ostream& __os, 
	     const hadronic_event<_Lv, hadronic_event_traits<0U,0U,0U> >& __a) 
  {
    for(int __i = __a.lower(); __i <= __a.upper(); __i++)
      __os<<"P["<<__i<<"] : "<<__a[__i]<<"  "<<__a[__i].mag2()<<std::endl;
    
    return __os;
  }

  //   Copies the basic_event object to a pre-allocated array. 
  template<class _Lv>
  void copy_to_array(typename _Lv::value_type *__p, 
		     const hadronic_event<_Lv, hadronic_event_traits<0U,0U,0U> >& __ev)
  {
    unsigned int __j = 0U;
    typename hadronic_event<_Lv, hadronic_event_traits<0U,0U,0U> >::const_iterator __iter, 
      __low = __ev.begin(), __up = __ev.end();
    
    for(__iter = __low; __iter < __up; __iter++) {
      __p[__j]    = __iter -> X();
      __p[__j+1U] = __iter -> Y();
      __p[__j+2U] = __iter -> Z();
      __p[__j+3U] = __iter -> T();
      __j += 4U;
    }
  }
  
}  // namespace nlo

#endif
