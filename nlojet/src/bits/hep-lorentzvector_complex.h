//  Copyright (C) 2005 Zoltan Nagy
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
#ifndef __NLO_HEP_LORENTZVECTOR_COMPLEX_H__
#define __NLO_HEP_LORENTZVECTOR_COMPLEX_H__ 1

//   Standard includes
#include <complex>

//   nlo includes
#include <bits/hep-threevector_complex.h>


namespace nlo {

  //   Spacializations
  class lorentzvector<std::complex<float> >;
  class lorentzvector<std::complex<double> >;
  class lorentzvector<std::complex<long double> >;
  
  template<>
  class lorentzvector<std::complex<float> > 
	: public threevector<std::complex<float> >
  {
  public:
	//  types
	typedef threevector<std::complex<float> > threevector_type;
	typedef threevector_type::value_type value_type;
	
	//       constructors
	lorentzvector() : _M_t(0.0f) {}
	
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

	lorentzvector& operator*=(const float& a) {
	  threevector_type::operator*=(a); _M_t *= a;
	  return *this;
	}
	
	lorentzvector& operator/=(const float& a) {
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
		
  protected:
	value_type  _M_t;
  };
  
  
  template<>
  class lorentzvector<std::complex<double> > 
	: public threevector<std::complex<double> >
  {
  public:
	//  types
	typedef threevector<std::complex<double> > threevector_type;
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

	lorentzvector(const lorentzvector<std::complex<float> >& a) 
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
	
	lorentzvector& operator*=(const double& a) {
	  threevector_type::operator*=(a); _M_t *= a;
	  return *this;
	}
	
	lorentzvector& operator/=(const double& a) {
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
	
  protected:
	value_type  _M_t;
  };
  
  
  
  template<>
  class lorentzvector<std::complex<long double> > 
	: public threevector<std::complex<long double> >
  {
  public:
	//  types
	typedef threevector<std::complex<long double> > threevector_type;
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
	
	lorentzvector(const lorentzvector<std::complex<float> >& a) 
	  : threevector_type(a), _M_t(a.T()) {} 
	
	lorentzvector(const lorentzvector<double>& a) 
	  : threevector_type(a), _M_t(a.T()) {} 

	lorentzvector(const lorentzvector<std::complex<double> >& a) 
	  : threevector_type(a), _M_t(a.T()) {} 

	lorentzvector(const lorentzvector<long double>& a) 
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
	
	lorentzvector& operator*=(const long double& a) {
	  threevector_type::operator*=(a); _M_t *= a;
	  return *this;
	}
	
	lorentzvector& operator/=(const long double& a) {
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
	
  protected:
	value_type  _M_t;
  };
    
  //     complex operators
  template<typename _Tp> 
  inline lorentzvector<std::complex<_Tp> > real(const lorentzvector<std::complex<_Tp> >& a) {
	return lorentzvector<std::complex<_Tp> >(std::real(a.X()), std::real(a.Y()), std::real(a.Z()), std::real(a.T()));
  }
  
  template<typename _Tp> 
  inline lorentzvector<std::complex<_Tp> > imag(const lorentzvector<std::complex<_Tp> >& a) {
	return lorentzvector<std::complex<_Tp> >(std::imag(a.X()), std::imag(a.Y()), std::imag(a.Z()), std::imag(a.T()));
  }
  
  template<typename _Tp> 
  inline lorentzvector<std::complex<_Tp> > conj(const lorentzvector<std::complex<_Tp> >& a) {
	return lorentzvector<std::complex<_Tp> >(std::conj(a.X()), std::conj(a.Y()), std::conj(a.Z()), std::conj(a.T()));
  }
  
  //     aditional operetors
  template<typename _Tp> inline lorentzvector<std::complex<_Tp> >
  operator+(const lorentzvector<std::complex<_Tp> >& a, const lorentzvector<_Tp>& b) {
	return lorentzvector<std::complex<_Tp> >(a) += b;
  }
  
  template<typename _Tp> inline lorentzvector<std::complex<_Tp> >
  operator+(const lorentzvector<_Tp>& b, const lorentzvector<std::complex<_Tp> >& a) {
	return lorentzvector<std::complex<_Tp> >(a) += b;
  }
  
  template<typename _Tp> inline lorentzvector<std::complex<_Tp> >
  operator-(const lorentzvector<std::complex<_Tp> >& a, const lorentzvector<_Tp>& b) {
	return lorentzvector<std::complex<_Tp> >(a) -= b;
  }
  
  template<typename _Tp> inline lorentzvector<std::complex<_Tp> >
  operator-(const lorentzvector<_Tp>& b, const lorentzvector<std::complex<_Tp> >& a) {
	return lorentzvector<std::complex<_Tp> >(a) -= b;
  }
  
  template<typename _Tp> inline lorentzvector<std::complex<_Tp> >
  operator*(const lorentzvector<std::complex<_Tp> >& a, const _Tp & b) {
	return lorentzvector<std::complex<_Tp> >(a) *= b;
  }
  
  template<typename _Tp> inline lorentzvector<std::complex<_Tp> >
  operator*(const _Tp & b, const lorentzvector<std::complex<_Tp> >& a) {
	return lorentzvector<std::complex<_Tp> >(a) *= b;
  }
  
  template<typename _Tp> inline lorentzvector<std::complex<_Tp> >
  operator/(const lorentzvector<std::complex<_Tp> >& a, const _Tp & b) {
	return lorentzvector<std::complex<_Tp> >(a) /= b;
  }
  
  template<typename _Tp> inline std::complex<_Tp> 
  operator*(const lorentzvector<std::complex<_Tp> >& a, const lorentzvector<_Tp>& b) {
	return a.T()*b.T() - a.X()*b.X() - a.Y()*b.Y() - a.Z()*b.Z();
  }
  
  template<typename _Tp> inline std::complex<_Tp> 
  operator*(const lorentzvector<_Tp>& b, const lorentzvector<std::complex<_Tp> >& a) {
	return a.T()*b.T() - a.X()*b.X() - a.Y()*b.Y() - a.Z()*b.Z();
  }
  
  template<typename _Tp> inline std::complex<_Tp> 
  dot(const lorentzvector<std::complex<_Tp> >& a, const lorentzvector<_Tp>& b) {
	return a.T()*b.T() - a.X()*b.X() - a.Y()*b.Y() - a.Z()*b.Z();
  }
  
  template<typename _Tp> inline std::complex<_Tp> 
  dot(const lorentzvector<_Tp>& b, const lorentzvector<std::complex<_Tp> >& a) {
	return a.T()*b.T() - a.X()*b.X() - a.Y()*b.Y() - a.Z()*b.Z();
  }
  
  template<typename _Tp> inline bool 
  operator==(const lorentzvector<std::complex<_Tp> >& a, const lorentzvector<_Tp>& b) {
	return real(a) == b && imag(a) == lorentzvector<_Tp>();
  }
  
  template<typename _Tp> inline bool 
  operator==(const lorentzvector<_Tp>& b, const lorentzvector<std::complex<_Tp> >& a) {
	return real(a) == b && imag(a) == lorentzvector<_Tp>();
  }
  
  template<typename _Tp> inline bool 
  operator!=(const lorentzvector<std::complex<_Tp> >& a, const lorentzvector<_Tp>& b) {
	return real(a) != b || imag(a) != lorentzvector<_Tp>();
  }
  
  template<typename _Tp> inline bool 
  operator!=(const lorentzvector<_Tp>& b, const lorentzvector<std::complex<_Tp> >& a) {
	return real(a) != b || imag(a) != lorentzvector<_Tp>();
  }

  //    spinor functions  <p|\gamma^\mu|q>
  lorentzvector<std::complex<float> > vector_pp(lorentzvector<float>, lorentzvector<float>);
  lorentzvector<std::complex<double> > vector_pp(lorentzvector<double>, lorentzvector<double>);
  lorentzvector<std::complex<long double> > vector_pp(lorentzvector<long double>, lorentzvector<long double>);
    
  inline lorentzvector<std::complex<float> > 
  vector_mm(const lorentzvector<float>& a, const lorentzvector<float>& b) {
	return vector_pp(b, a);
  }
  
  inline lorentzvector<std::complex<double> > 
  vector_mm(const lorentzvector<double>& a, const lorentzvector<double>& b) {
	return vector_pp(b, a);
  }

  inline lorentzvector<std::complex<long double> > 
  vector_mm(const lorentzvector<long double>& a, const lorentzvector<long double>& b) {
	return vector_pp(b, a);
  }
}

#endif

