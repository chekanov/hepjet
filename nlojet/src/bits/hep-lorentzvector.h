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
#ifndef __NLO_HEP_LORENTZVECTOR_H__
#define __NLO_HEP_LORENTZVECTOR_H__ 1

//   Standard includes
#include <complex>

//   nlo includes
#include <bits/hep-threevector.h>


namespace nlo {
  
  
  template<typename _Tp>
  class lorentzvector : public threevector<_Tp> 
  {
  public:
    //  types
    typedef threevector<_Tp> threevector_type;
    typedef typename threevector_type::value_type value_type;

    //       constructors
    lorentzvector();
    lorentzvector(const _Tp&, const _Tp&, const _Tp&, const _Tp&);
    lorentzvector(const threevector_type&, const _Tp&);
    lorentzvector(const lorentzvector&);
    template<typename _Xp> explicit lorentzvector(const lorentzvector<_Xp>&);
    
    //   assignments
    lorentzvector& operator=(const lorentzvector&);
    
    template<typename _Xp> 
    lorentzvector<_Tp>& operator=(const lorentzvector<_Xp>&);
    
    //   computed assignments
    lorentzvector& operator+=(const lorentzvector&);
    lorentzvector& operator-=(const lorentzvector&);
    lorentzvector& operator*=(const value_type&);
    lorentzvector& operator/=(const value_type&);
    
    template<typename _Xp> 
    lorentzvector<_Tp>& operator+=(const lorentzvector<_Xp>&);
    
    template<typename _Xp> 
    lorentzvector<_Tp>& operator-=(const lorentzvector<_Xp>&);
    
    //   elements access
    const value_type& T() const;
    value_type& T();
  };
  
  
  //   Spacializations
  
  template<>
  class lorentzvector<float> : public threevector<float>
  {
  public:
    //  types
    typedef threevector<float> threevector_type;
    typedef threevector_type::value_type value_type;
    
    //       constructors
    lorentzvector() : _M_t(0.0f) {}
    
    lorentzvector(const value_type& x, const value_type& y, const value_type& z, const value_type& t)
      : threevector_type(x, y, z), _M_t(t) {}
    
    lorentzvector(const threevector_type& v, const value_type& t)
      : threevector_type(v), _M_t(t) {} 
 
    lorentzvector(const lorentzvector& a)
      : threevector_type(a), _M_t(a._M_t) {} 
    
    template<typename _Xp> 
    explicit lorentzvector(const lorentzvector<_Xp>& a) 
	  : threevector_type(a), _M_t(a.T()) {} 
    
    //   assignments
    lorentzvector& operator=(const lorentzvector& a) {
      if(this != &a) { 
		threevector_type::operator=(a); _M_t = a._M_t;
	  }
      return *this;
    }
    
    template<typename _Xp> 
    lorentzvector<value_type>& operator=(const lorentzvector<_Xp>& a) {
      threevector_type::operator=(a); _M_t = a.T();
      return *this;
    }
   
    //   computed assignments
    lorentzvector& operator+=(const lorentzvector& a) {
      threevector_type::operator+=(a); _M_t += a._M_t;
      return *this;
    }

    lorentzvector& operator-=(const lorentzvector& a) {
      threevector_type::operator-=(a); _M_t -= a._M_t;
      return *this;
    }

    lorentzvector& operator*=(const value_type& a) {
      threevector_type::operator*=(a); _M_t *= a;
      return *this;
    }

    lorentzvector& operator/=(const value_type& a) {
      threevector_type::operator/=(a); _M_t /= a;
      return *this;
    }
        
    template<typename _Xp> 
    lorentzvector<value_type>& operator+=(const lorentzvector<_Xp>& a) {
      threevector_type::operator+=(a); _M_t += a.T();
      return *this;
    }
  
    template<typename _Xp> 
    lorentzvector<value_type>& operator-=(const lorentzvector<_Xp>& a) {
      threevector_type::operator-=(a); _M_t -= a.T();
      return *this;
    }
    
    //   elements access
    const value_type& T() const { return _M_t;}
    value_type&  T() { return _M_t;}

    //  member functions
    value_type plus () const { return _M_t + _M_z;}
    value_type minus() const { return _M_t - _M_z;}

    value_type rapidity() const 
      { return 0.5f*std::log(plus()/minus());}

    value_type prapidity() const 
      { return -std::log(std::tan(0.5f*theta()));}

    value_type mag2() const 
      { return _M_t*_M_t - threevector_type::mag2();}  
   
    value_type mag() const {
      value_type p =  mag2();
      return p >= 0.0f ? std::sqrt(p) : -std::sqrt(-p);
    }
    
    threevector_type boostVector() const { 
      return threevector_type(*this) /= _M_t;
    }
    
    //  Lorentz transformations
    void boost(const value_type&, const value_type&, const value_type&); 
    void boost(const threevector_type& a) { 
      boost(a.X(), a.Y(), a.Z());
    }

  protected:
    value_type  _M_t;
  };
  
  
  template<>
  class lorentzvector<double> : public threevector<double>
  {
  public:
	//  types
	typedef threevector<double> threevector_type;
	typedef threevector_type::value_type value_type;
	
	//       constructors
	lorentzvector() : _M_t(0.0) {}
	
	lorentzvector(const value_type& x, const value_type& y, const value_type& z, const value_type& t)
	  : threevector_type(x, y, z), _M_t(t) {}
	
	lorentzvector(const threevector_type& v, const value_type& t)
	  : threevector_type(v), _M_t(t) {} 
	
	lorentzvector(const lorentzvector& a)
	  : threevector_type(a), _M_t(a._M_t) {} 
	
	lorentzvector(const lorentzvector<float>& a) 
	  : threevector_type(a), _M_t(a.T()) {} 

	template<typename _Xp> 
	explicit lorentzvector(const lorentzvector<_Xp>& a) 
	  : threevector_type(a), _M_t(a.T()) {} 
	
	//   assignments
	lorentzvector& operator=(const lorentzvector& a) {
	  if(this != &a) { 
		threevector_type::operator=(a); _M_t = a._M_t;
	  }
	  return *this;
	}
	
	template<typename _Xp> 
	  lorentzvector<value_type>& operator=(const lorentzvector<_Xp>& a) {
		threevector_type::operator=(a); _M_t = a.T();
		return *this;
	  }
	
	//   computed assignments
	lorentzvector& operator+=(const lorentzvector& a) {
	  threevector_type::operator+=(a); _M_t += a._M_t;
	  return *this;
	}
	
	lorentzvector& operator-=(const lorentzvector& a) {
	  threevector_type::operator-=(a); _M_t -= a._M_t;
	  return *this;
	}
	
	lorentzvector& operator*=(const value_type& a) {
	  threevector_type::operator*=(a); _M_t *= a;
	  return *this;
	}
	
	lorentzvector& operator/=(const value_type& a) {
	  threevector_type::operator/=(a); _M_t /= a;
	  return *this;
	}
	
	template<typename _Xp> 
	lorentzvector<value_type>& operator+=(const lorentzvector<_Xp>& a) {
	  threevector_type::operator+=(a); _M_t += a.T();
	  return *this;
	}
	
	template<typename _Xp> 
	lorentzvector<value_type>& operator-=(const lorentzvector<_Xp>& a) {
	  threevector_type::operator-=(a); _M_t -= a.T();
	  return *this;
	}
	
	//   elements access
	const value_type& T() const { return _M_t;}
	value_type&  T() { return _M_t;}
	
	//  member functions
	value_type plus () const { return _M_t + _M_z;}
	value_type minus() const { return _M_t - _M_z;}
	
	value_type rapidity() const 
	  { return 0.5*std::log(plus()/minus());}
	
	value_type prapidity() const 
	  { return -std::log(std::tan(0.5*theta()));}
	
	value_type mag2() const 
	  { return _M_t*_M_t - threevector_type::mag2();}  
	
	value_type mag() const {
	  value_type p =  mag2();
	  return p >= 0.0 ? std::sqrt(p) : -std::sqrt(-p);
	}
	
	threevector_type boostVector() const { 
	  return threevector_type(*this) /= _M_t;
	}
	
	//  Lorentz transformations
	void boost(const value_type&, const value_type&, const value_type&); 
	void boost(const threevector_type& a) { 
	  boost(a.X(), a.Y(), a.Z());
	}
	
  protected:
	value_type  _M_t;
  };
  
  
  
  template<>
  class lorentzvector<long double> : public threevector<long double>
  {
  public:
	//  types
	typedef threevector<long double> threevector_type;
	typedef threevector_type::value_type value_type;
	  
	//       constructors
	lorentzvector() : _M_t(0.0L) {}
	
	lorentzvector(const value_type& x, const value_type& y, const value_type& z, const value_type& t)
	  : threevector_type(x, y, z), _M_t(t) {}
	
	lorentzvector(const threevector_type& v, const value_type& t)
	  : threevector_type(v), _M_t(t) {} 
	
	lorentzvector(const lorentzvector& a)
	  : threevector_type(a), _M_t(a._M_t) {} 
	
	lorentzvector(const lorentzvector<float>& a) 
	  : threevector_type(a), _M_t(a.T()) {} 

	lorentzvector(const lorentzvector<double>& a) 
	  : threevector_type(a), _M_t(a.T()) {} 
	
	template<typename _Xp> 
	  explicit lorentzvector(const lorentzvector<_Xp>& a) 
	  : threevector_type(a), _M_t(a.T()) {} 
	
	//   assignments
	lorentzvector& operator=(const lorentzvector& a) {
	  if(this != &a) { 
		threevector_type::operator=(a); _M_t = a._M_t;
	  }
	  return *this;
	}
	
	template<typename _Xp> 
	lorentzvector<value_type>& operator=(const lorentzvector<_Xp>& a) {
	  threevector_type::operator=(a); _M_t = a.T();
	  return *this;
	}
	
	//   computed assignments
	lorentzvector& operator+=(const lorentzvector& a) {
	  threevector_type::operator+=(a); _M_t += a._M_t;
	  return *this;
	}
	
	lorentzvector& operator-=(const lorentzvector& a) {
	  threevector_type::operator-=(a); _M_t -= a._M_t;
	  return *this;
	}
	
	lorentzvector& operator*=(const value_type& a) {
	  threevector_type::operator*=(a); _M_t *= a;
	  return *this;
	}
	
	lorentzvector& operator/=(const value_type& a) {
	  threevector_type::operator/=(a); _M_t /= a;
	  return *this;
	}
	
	template<typename _Xp> 
	lorentzvector<value_type>& operator+=(const lorentzvector<_Xp>& a) {
	  threevector_type::operator+=(a); _M_t += a.T();
	  return *this;
	}
	
	template<typename _Xp> 
	lorentzvector<value_type>& operator-=(const lorentzvector<_Xp>& a) {
	  threevector_type::operator-=(a); _M_t -= a.T();
	  return *this;
	}
	
	//   elements access
	const value_type& T() const { return _M_t;}
	value_type&  T() { return _M_t;}
	
	//  member functions
	value_type plus () const { return _M_t + _M_z;}
	value_type minus() const { return _M_t - _M_z;}
	
	value_type rapidity() const 
	  { return 0.5L*std::log(plus()/minus());}
	
	value_type prapidity() const 
	  { return -std::log(std::tan(0.5L*theta()));}
	
	value_type mag2() const 
	  { return _M_t*_M_t - threevector_type::mag2();}  
	
	value_type mag() const {
	  value_type p =  mag2();
	  return p >= 0.0L ? std::sqrt(p) : -std::sqrt(-p);
	}
	
	threevector_type boostVector() const { 
	  return threevector_type(*this) /= _M_t;
	}
	
	//  Lorentz transformations
	void boost(const value_type&, const value_type&, const value_type&); 
	void boost(const threevector_type& a) { 
	  boost(a.X(), a.Y(), a.Z());
	}
	
  protected:
	value_type  _M_t;
  };
  

  //   unary operations
  template<typename _Tp> 
  inline lorentzvector<_Tp> operator+(const lorentzvector<_Tp>& a) {
    return a;
  }

  template<typename _Tp> 
  inline lorentzvector<_Tp> operator-(const lorentzvector<_Tp>& a) {
    return lorentzvector<_Tp>(-a.X(), -a.Y(), -a.Z(), -a.T());
  }
  
  //     other operators
  template<typename _Tp> inline lorentzvector<_Tp> 
  operator+(const lorentzvector<_Tp>& a, const lorentzvector<_Tp>& b) {
    return lorentzvector<_Tp>(a) += b;
  }
  
  template<typename _Tp> inline lorentzvector<_Tp> 
  operator-(const lorentzvector<_Tp>& a, const lorentzvector<_Tp>& b) {
    return lorentzvector<_Tp>(a) -= b;
  }
  
  template<typename _Tp> inline lorentzvector<_Tp> 
  operator*(const lorentzvector<_Tp>& a, const _Tp & b) {
    return lorentzvector<_Tp>(a) *= b;
  }

  template<typename _Tp> inline lorentzvector<_Tp> 
  operator*(const _Tp & b, const lorentzvector<_Tp>& a) {
    return lorentzvector<_Tp>(a) *= b;
  }
  
  template<typename _Tp> inline lorentzvector<_Tp> 
  operator/(const lorentzvector<_Tp>& a, const _Tp & b) {
    return lorentzvector<_Tp>(a) /= b;
  }

  template<typename _Tp> inline _Tp 
  operator*(const lorentzvector<_Tp>& a, const lorentzvector<_Tp>& b) {
    return a.T()*b.T() - a.X()*b.X() - a.Y()*b.Y() - a.Z()*b.Z();
  }

  template<typename _Tp> inline bool 
  operator==(const lorentzvector<_Tp>& a, const lorentzvector<_Tp>& b) {
    return a.X() == b.X() && a.Y() == b.Y() && a.Z()== b.Z() && a.T()== b.T();
  }
  
  
  template<typename _Tp> inline bool 
  operator!=(const lorentzvector<_Tp>& a, const lorentzvector<_Tp>& b) {
    return a.X() != b.X() || a.Y() != b.Y() || a.Z() != b.Z() || a.T() != b.T();
  }
  
  template<typename _Tp> inline _Tp 
  dot(const lorentzvector<_Tp>& a, const lorentzvector<_Tp>& b) {
    return a.T()*b.T() - a.X()*b.X() - a.Y()*b.Y() - a.Z()*b.Z();
  }
 
  //    spinor operations
  std::complex<float> scalar_mp(lorentzvector<float>, lorentzvector<float>);
  std::complex<double> scalar_mp(lorentzvector<double>, lorentzvector<double>);
  std::complex<long double> scalar_mp(lorentzvector<long double>, lorentzvector<long double>);

  inline std::complex<float> 
  scalar_pm(const lorentzvector<float>& a, const lorentzvector<float>& b) {
    return (2.0f*dot(a, b))/scalar_mp(b, a);
  }

  //   spinor products
  inline std::complex<double> 
  scalar_pm(const lorentzvector<double>& a, const lorentzvector<double>& b) {
    return (2.0*dot(a, b))/scalar_mp(b, a);
  }

  inline std::complex<long double> 
  scalar_pm(const lorentzvector<long double>& a, const lorentzvector<long double>& b) {
    return (2.0L*dot(a, b))/scalar_mp(b, a);
  }
  
  //    I/O operations  
  template<typename _Tp>
  std::ostream& operator<<(std::ostream& os, const lorentzvector<_Tp>& q)  {
	return os<<"("<<q.X()<<", "<<q.Y()<<", "<<q.Z()<<"; "<<q.T() <<")";
  }
  
  template<typename _Tp>
  std::istream& operator>>(std::istream& is, lorentzvector<_Tp>& a) {
	char delin;
	is>>delin>>a.X()>>delin>>a.Y()>>delin>>a.Z()>>delin>>a.T()>>delin;
	return is;
  }
}   // namespace nlo

#endif
