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
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#ifndef __NLO_HEP_THREEVECTOR_COMPLEX_H__
#define __NLO_HEP_THREEVECTOR_COMPLEX_H__ 1


//   Standard includes
#include <cmath>
#include <iostream>
#include <complex>


namespace nlo {
  
  //    Specializations
  class threevector<std::complex<float> >;
  class threevector<std::complex<double> >;
  class threevector<std::complex<long double> >;
  
  template<>
  class threevector<std::complex<float> > 
  {
  public:
	//        types
	typedef std::complex<float> value_type;
	
	//  constructors
	threevector() 
	  : _M_x(0.0f), _M_y(0.0f), _M_z(0.0f) {}
	
	threevector(const value_type& x, const value_type& y, const value_type& z)
	  : _M_x(x), _M_y(y), _M_z(z) {}
	
	threevector(const threevector& x)
	  : _M_x(x._M_x), _M_y(x._M_y), _M_z(x._M_z) {}
	
	threevector(const threevector<float>& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}

	template<typename _Xp> 
	explicit threevector(const threevector<_Xp>& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}
	
	//    assignments
	threevector& operator=(const threevector& a) 
	{
	  if(this != &a) {
		_M_x = a._M_x; _M_y = a._M_y; _M_z = a._M_z;
	  }
	  return *this;
	}
	
	template<typename _Xp> 
	threevector<value_type>& operator=(const threevector<_Xp>& x) {
	  _M_x = x.X(); _M_y = x.Y(); _M_z = x.Z();
	  return *this;
	}
	
	//  elements access
	const value_type& X() const { return _M_x;}
	const value_type& Y() const { return _M_y;}
	const value_type& Z() const { return _M_z;} 
	
	value_type& X() { return _M_x;}
	value_type& Y() { return _M_y;}
	value_type& Z() { return _M_z;} 
	
	//   computed assignments
	threevector& operator+=(const threevector& a) {
	  _M_x += a._M_x; _M_y += a._M_y; _M_z += a._M_z;
	  return *this;
	}
	
	threevector& operator-=(const threevector& a) {
	  _M_x -= a._M_x; _M_y -= a._M_y; _M_z -= a._M_z;
	  return *this;
	}
	
	threevector& operator*=(const value_type& a) {
	  _M_x *= a; _M_y *= a; _M_z *= a;
	  return *this;
	}
	
	threevector& operator/=(const value_type& a) {
	  _M_x /= a; _M_y /= a; _M_z /= a; 
	  return *this ;
	}

	threevector& operator*=(const float& a) {
	  _M_x *= a; _M_y *= a; _M_z *= a;
	  return *this;
	}
	
	threevector& operator/=(const float& a) {
	  _M_x /= a; _M_y /= a; _M_z /= a; 
	  return *this ;
	}
	
	template<typename _Xp> 
	threevector<value_type>& operator+=(const threevector<_Xp>& x) {
	  _M_x += x.X(); _M_y += x.Y(); _M_z += x.Z();
	  return *this;
	}
	
	template<typename _Xp> 
	threevector<value_type>& operator-=(const threevector<_Xp>& x) {
	  _M_x -= x.X(); _M_y -= x.Y(); _M_z -= x.Z();
	  return *this;
	}
	
  protected:
	//        data member
	value_type _M_x, _M_y, _M_z; 
  };
  
  
  
  template<>
  class threevector<std::complex<double> > 
  {
  public:
	//        types
	typedef std::complex<double> value_type;
	
	//  constructors
	threevector() 
	  : _M_x(0.0), _M_y(0.0), _M_z(0.0) {}
	
	threevector(const value_type& x, const value_type& y, const value_type& z)
	  : _M_x(x), _M_y(y), _M_z(z) {}
	
	threevector(const threevector& x)
	  : _M_x(x._M_x), _M_y(x._M_y), _M_z(x._M_z) {}
	
	threevector(const threevector<float>& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}

	threevector(const threevector<std::complex<float> >& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}

	threevector(const threevector<double>& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}
	
	template<typename _Xp> 
	explicit threevector(const threevector<_Xp>& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}
	
	//    assignments
	threevector& operator=(const threevector& a) 
	{
	  if(this != &a) {
		_M_x = a._M_x; _M_y = a._M_y; _M_z = a._M_z;
	  }
	  return *this;
	}
	
	template<typename _Xp> 
	threevector<value_type>& operator=(const threevector<_Xp>& x) {
	  _M_x = x.X(); _M_y = x.Y(); _M_z = x.Z();
	  return *this;
	}
	
	//  elements access
	const value_type& X() const { return _M_x;}
	const value_type& Y() const { return _M_y;}
	const value_type& Z() const { return _M_z;} 
	
	value_type& X() { return _M_x;}
	value_type& Y() { return _M_y;}
	value_type& Z() { return _M_z;} 
	
	//   computed assignments
	threevector& operator+=(const threevector& a) {
	  _M_x += a._M_x; _M_y += a._M_y; _M_z += a._M_z;
	  return *this;
	}
	
	threevector& operator-=(const threevector& a) {
	  _M_x -= a._M_x; _M_y -= a._M_y; _M_z -= a._M_z;
	  return *this;
	}
	
	threevector& operator*=(const value_type& a) {
	  _M_x *= a; _M_y *= a; _M_z *= a;
	  return *this;
	}
	
	threevector& operator/=(const value_type& a) {
	  _M_x /= a; _M_y /= a; _M_z /= a; 
	  return *this ;
	}
	
	threevector& operator*=(const double& a) {
	  _M_x *= a; _M_y *= a; _M_z *= a;
	  return *this;
	}
	
	threevector& operator/=(const double& a) {
	  _M_x /= a; _M_y /= a; _M_z /= a; 
	  return *this ;
	}
	
	template<typename _Xp> 
	threevector<value_type>& operator+=(const threevector<_Xp>& x) {
	  _M_x += x.X(); _M_y += x.Y(); _M_z += x.Z();
	  return *this;
	}
	
	template<typename _Xp> 
	threevector<value_type>& operator-=(const threevector<_Xp>& x) {
	  _M_x -= x.X(); _M_y -= x.Y(); _M_z -= x.Z();
	  return *this;
	}
	
  protected:
	//        data member
	value_type _M_x, _M_y, _M_z; 
  };
  
  
  template<>
  class threevector<std::complex<long double> > 
  {
  public:
	//        types
	typedef std::complex<long double> value_type;
	
	//  constructors
	threevector() 
	  : _M_x(0.0L), _M_y(0.0L), _M_z(0.0L) {}
	
	threevector(const value_type& x, const value_type& y, const value_type& z)
	  : _M_x(x), _M_y(y), _M_z(z) {}
	
	threevector(const threevector& x)
	  : _M_x(x._M_x), _M_y(x._M_y), _M_z(x._M_z) {}
	
	threevector(const threevector<float>& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}
	
	threevector(const threevector<std::complex<float> >& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}
	
	threevector(const threevector<double>& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}
	
	threevector(const threevector<std::complex<double> >& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}

	threevector(const threevector<long double>& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}
	
	template<typename _Xp> 
	explicit threevector(const threevector<_Xp>& x)
	  : _M_x(x.X()), _M_y(x.Y()), _M_z(x.Z()) {}
	
	//    assignments
	threevector& operator=(const threevector& a) 
	{
	  if(this != &a) {
		_M_x = a._M_x; _M_y = a._M_y; _M_z = a._M_z;
	  }
	  return *this;
	}
	
	template<typename _Xp> 
	threevector<value_type>& operator=(const threevector<_Xp>& x) {
	  _M_x = x.X(); _M_y = x.Y(); _M_z = x.Z();
	  return *this;
	}
	
	//  elements access
	const value_type& X() const { return _M_x;}
	const value_type& Y() const { return _M_y;}
	const value_type& Z() const { return _M_z;} 
	
	value_type& X() { return _M_x;}
	value_type& Y() { return _M_y;}
	value_type& Z() { return _M_z;} 
	
	//   computed assignments
	threevector& operator+=(const threevector& a) {
	  _M_x += a._M_x; _M_y += a._M_y; _M_z += a._M_z;
	  return *this;
	}
	
	threevector& operator-=(const threevector& a) {
	  _M_x -= a._M_x; _M_y -= a._M_y; _M_z -= a._M_z;
	  return *this;
	}
	
	threevector& operator*=(const value_type& a) {
	  _M_x *= a; _M_y *= a; _M_z *= a;
	  return *this;
	}
	
	threevector& operator/=(const value_type& a) {
	  _M_x /= a; _M_y /= a; _M_z /= a; 
	  return *this ;
	}
	
	threevector& operator*=(const long double& a) {
	  _M_x *= a; _M_y *= a; _M_z *= a;
	  return *this;
	}
	
	threevector& operator/=(const long double& a) {
	  _M_x /= a; _M_y /= a; _M_z /= a; 
	  return *this ;
	}
	
	template<typename _Xp> 
	threevector<value_type>& operator+=(const threevector<_Xp>& x) {
	  _M_x += x.X(); _M_y += x.Y(); _M_z += x.Z();
	  return *this;
	}
	
	template<typename _Xp> 
	threevector<value_type>& operator-=(const threevector<_Xp>& x) {
	  _M_x -= x.X(); _M_y -= x.Y(); _M_z -= x.Z();
	  return *this;
	}
	
  protected:
	//        data member
	value_type _M_x, _M_y, _M_z; 
  };
  
  //     other operators
  template<typename _Tp> 
  inline threevector<std::complex<_Tp> > real(const threevector<std::complex<_Tp> >& a) {
	return threevector<std::complex<_Tp> >(std::real(a.X()), std::real(a.Y()), std::real(a.Z()));
  }
  
  template<typename _Tp> 
  inline threevector<std::complex<_Tp> > imag(const threevector<std::complex<_Tp> >& a) {
	return threevector<std::complex<_Tp> >(std::imag(a.X()), std::imag(a.Y()), std::imag(a.Z()));
  }
  
  template<typename _Tp> 
  inline threevector<std::complex<_Tp> > conj(const threevector<std::complex<_Tp> >& a) {
	return threevector<std::complex<_Tp> >(std::conj(a.X()), std::conj(a.Y()), std::conj(a.Z()));
  }
  
  //     aditional operetors
  template<typename _Tp> inline threevector<std::complex<_Tp> >
  operator+(const threevector<std::complex<_Tp> >& a, const threevector<_Tp>& b) {
	return threevector<std::complex<_Tp> >(a) += b;
  }

  template<typename _Tp> inline threevector<std::complex<_Tp> >
  operator+(const threevector<_Tp>& b, const threevector<std::complex<_Tp> >& a) {
	return threevector<std::complex<_Tp> >(a) += b;
  }
  
  template<typename _Tp> inline threevector<std::complex<_Tp> >
  operator-(const threevector<std::complex<_Tp> >& a, const threevector<_Tp>& b) {
	return threevector<std::complex<_Tp> >(a) -= b;
  }
  
  template<typename _Tp> inline threevector<std::complex<_Tp> >
  operator-(const threevector<_Tp>& b, const threevector<std::complex<_Tp> >& a) {
	return threevector<std::complex<_Tp> >(a) -= b;
  }
  
  template<typename _Tp> inline threevector<std::complex<_Tp> >
  operator*(const threevector<std::complex<_Tp> >& a, const _Tp & b) {
	return threevector<std::complex<_Tp> >(a) *= b;
  }
  
  template<typename _Tp> inline threevector<std::complex<_Tp> >
  operator*(const _Tp & b, const threevector<std::complex<_Tp> >& a) {
	return threevector<std::complex<_Tp> >(a) *= b;
  }
  
  template<typename _Tp> inline threevector<std::complex<_Tp> >
  operator/(const threevector<std::complex<_Tp> >& a, const _Tp & b) {
	return threevector<std::complex<_Tp> >(a) /= b;
  }
  
  template<typename _Tp> inline std::complex<_Tp> 
  operator*(const threevector<std::complex<_Tp> >& a, const threevector<_Tp>& b) {
	return a.X()*b.X() + a.Y()*b.Y() + a.Z()*b.Z();
  }
  
  template<typename _Tp> inline std::complex<_Tp> 
  operator*(const threevector<_Tp>& b, const threevector<std::complex<_Tp> >& a) {
	return a.X()*b.X() + a.Y()*b.Y() + a.Z()*b.Z();
  }
  
  template<typename _Tp> inline std::complex<_Tp> 
  dot(const threevector<std::complex<_Tp> >& a, const threevector<_Tp>& b) {
	return a.X()*b.X() + a.Y()*b.Y() + a.Z()*b.Z();
  }
  
  template<typename _Tp> inline std::complex<_Tp> 
  dot(const threevector<_Tp>& b, const threevector<std::complex<_Tp> >& a) {
	return a.X()*b.X() + a.Y()*b.Y() + a.Z()*b.Z();
  }
  
  template<typename _Tp> inline bool 
  operator==(const threevector<std::complex<_Tp> >& a, const threevector<_Tp>& b) {
	return real(a) == b && imag(a) == threevector<_Tp>();
  }
  
  template<typename _Tp> inline bool 
  operator==(const threevector<_Tp>& b, const threevector<std::complex<_Tp> >& a) {
	return real(a) == b && imag(a) == threevector<_Tp>();
  }
  
  template<typename _Tp> inline bool 
  operator!=(const threevector<std::complex<_Tp> >& a, const threevector<_Tp>& b) {
	return real(a) != b || imag(a) != threevector<_Tp>();
  }
  
  template<typename _Tp> inline bool 
  operator!=(const threevector<_Tp>& b, const threevector<std::complex<_Tp> >& a) {
	return real(a) != b || imag(a) != threevector<_Tp>();
  }
  
  template<typename _Tp> inline threevector<std::complex<_Tp> >
  cross(const threevector<std::complex<_Tp> >& a, const threevector<_Tp>& b) 
  {
	return threevector<_Tp>(a.Y()*b.Z() - a.Z()*b.Y(),
							a.Z()*b.X() - a.X()*b.Z(),
							a.X()*b.Y() - a.Y()*b.X());
  }
  
  template<typename _Tp> inline threevector<std::complex<_Tp> >
  cross(const threevector<_Tp>& b, const threevector<std::complex<_Tp> >& a) 
  {
	return threevector<_Tp>(a.Y()*b.Z() - a.Z()*b.Y(),
							a.Z()*b.X() - a.X()*b.Z(),
							a.X()*b.Y() - a.Y()*b.X());
  }
}

#endif

