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
#ifndef __NLO_HEP_THREEVECTOR_H__
#define __NLO_HEP_THREEVECTOR_H__ 1


//   Standard includes
#include <cmath>
#include <iostream>


namespace nlo {
  
  
  template<typename _Tp>
  class threevector 
  {
  public:
    //        types
    typedef _Tp value_type;
    
    //  constructors
    threevector();
    threevector(const _Tp&, const _Tp&, const _Tp&);
    threevector(const threevector&);
    template<typename _Xp> explicit threevector(const threevector<_Xp>&);
    
    //   assignments
    threevector& operator=(const threevector&);
    template<typename _Xp> threevector<_Tp>& operator=(const threevector<_Xp>&);
    
    threevector& operator+=(const threevector&);
    threevector& operator-=(const threevector&);
    threevector& operator*=(const value_type&);
    threevector& operator/=(const value_type&);
    
    template<typename _Xp> threevector<_Tp>& operator+=(const threevector<_Xp>&);
    template<typename _Xp> threevector<_Tp>& operator-=(const threevector<_Xp>&);
    
    //   elements access
    const value_type& X() const;
    const value_type& Y() const;
    const value_type& Z() const;
    
    value_type& X();
    value_type& Y();
    value_type& Z();
  };


  //    Specializations
  
  template<>
  class threevector<float> 
  {
  public:
    //        types
    typedef float value_type;
    
    //  constructors
    threevector() 
      : _M_x(0.0f), _M_y(0.0f), _M_z(0.0f) {}
  
    threevector(const value_type& x, const value_type& y, const value_type& z)
	  : _M_x(x), _M_y(y), _M_z(z) {}
    
    threevector(const threevector& x)
      : _M_x(x._M_x), _M_y(x._M_y), _M_z(x._M_z) {}
    
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

    //  the unit vector parallel to this
    threevector unit() const {
      value_type tot = this -> mag2();
      threevector p(*this);
      return tot > 0.0f ? p *= (1.0f/std::sqrt(tot)) : p;
    }
    
    //        magnitude and the transvers component
    value_type mag () const { return std::sqrt(this -> mag2());}
    value_type perp() const { return std::sqrt(this -> perp2());}
    
    //        magnitude and the transvers component squared
    value_type mag2 () const { return _M_x*_M_x + _M_y*_M_y + _M_z*_M_z;}
    value_type perp2() const { return _M_x*_M_x + _M_y*_M_y;}
    
    //        azimuth and polar angles
    value_type phi() const {
      return _M_x == 0.0f && _M_y == 0.0f ? 0.0f : std::atan2(_M_y,_M_x);
    }
    
    value_type theta() const {
      value_type p = this -> perp();
      return p == 0.0f && _M_z == 0.0f ? 0.0f : std::atan2(p, _M_z);
    }
    
    value_type cosTheta() const {
      value_type ptot = this -> mag();
      return ptot == 0.0f ? 1.0f : _M_z/ptot;
    }
    
    //        rotations   
    void rotateX(const value_type&);
    void rotateY(const value_type&);
    void rotateZ(const value_type&);
    void rotate (const value_type&, const threevector &);
    
  protected:
    //        data member
    value_type _M_x, _M_y, _M_z; 
  };
  
  
  template<>
  class threevector<double> 
  {
  public:
	//        types
	typedef double value_type;
	
	//  constructors
	threevector() 
	  : _M_x(0.0), _M_y(0.0), _M_z(0.0) {}
	
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
	
	//  the unit vector parallel to this
	threevector unit() const {
	  value_type tot = this -> mag2();
	  threevector p(*this);
	  return tot > 0.0 ? p *= (1.0/std::sqrt(tot)) : p;
	}
	
	//        magnitude and the transvers component
	value_type mag () const { return std::sqrt(this -> mag2());}
	value_type perp() const { return std::sqrt(this -> perp2());}
	
	//        magnitude and the transvers component squared
	value_type mag2 () const { return _M_x*_M_x + _M_y*_M_y + _M_z*_M_z;}
	value_type perp2() const { return _M_x*_M_x + _M_y*_M_y;}
	
	//        azimuth and polar angles
	value_type phi() const {
	  return _M_x == 0.0 && _M_y == 0.0 ? 0.0 : std::atan2(_M_y,_M_x);
	}
	
	value_type theta() const {
	  value_type p = this -> perp();
	  return p == 0.0 && _M_z == 0.0 ? 0.0 : std::atan2(p, _M_z);
	}
	
	value_type cosTheta() const {
	  value_type ptot = this -> mag();
	  return ptot == 0.0 ? 1.0 : _M_z/ptot;
	}
	
	//        rotations   
	void rotateX(const value_type&);
	void rotateY(const value_type&);
	void rotateZ(const value_type&);
	void rotate (const value_type&, const threevector &);
	
  protected:
	//        data member
	value_type _M_x, _M_y, _M_z; 
  };
  

  template<>
  class threevector<long double> 
  {
  public:
	//        types
	typedef long double value_type;
	
	//  constructors
	threevector() 
	  : _M_x(0.0L), _M_y(0.0L), _M_z(0.0L) {}
	
	threevector(const value_type& x, const value_type& y, const value_type& z)
	  : _M_x(x), _M_y(y), _M_z(z) {}
	
	threevector(const threevector& x)
	  : _M_x(x._M_x), _M_y(x._M_y), _M_z(x._M_z) {}
	
	threevector(const threevector<float>& x)
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
	
	//  the unit vector parallel to this
	threevector unit() const {
	  value_type tot = this -> mag2();
	  threevector p(*this);
	  return tot > 0.0L ? p *= (1.0L/std::sqrt(tot)) : p;
	}
	
	//        magnitude and the transvers component
	value_type mag () const { return std::sqrt(this -> mag2());}
	value_type perp() const { return std::sqrt(this -> perp2());}
	
	//        magnitude and the transvers component squared
	value_type mag2 () const { return _M_x*_M_x + _M_y*_M_y + _M_z*_M_z;}
	value_type perp2() const { return _M_x*_M_x + _M_y*_M_y;}
	
	//        azimuth and polar angles
	value_type phi() const {
	  return _M_x == 0.0L && _M_y == 0.0L ? 0.0L : std::atan2(_M_y,_M_x);
	}
	
	value_type theta() const {
	  value_type p = this -> perp();
	  return p == 0.0L && _M_z == 0.0L ? 0.0L : std::atan2(p, _M_z);
	}
	
	value_type cosTheta() const {
	  value_type ptot = this -> mag();
	  return ptot == 0.0L ? 1.0L : _M_z/ptot;
	}
	
	//        rotations   
	void rotateX(const value_type&);
	void rotateY(const value_type&);
	void rotateZ(const value_type&);
	void rotate (const value_type&, const threevector &);
	
  protected:
	//        data member
	value_type   _M_x, _M_y, _M_z; 
  };
  
  //  unary operations
  template<typename _Tp> 
  inline threevector<_Tp> operator+(const threevector<_Tp>& x) { 
    return x;
  }
  
  template<typename _Tp> 
  inline threevector<_Tp> operator-(const threevector<_Tp>& x) {
    return threevector<_Tp>(-x.X(), -x.Y(), -x.Z());
  }
  
  //     other operators
  template<typename _Tp> inline threevector<_Tp> 
  operator+(const threevector<_Tp>& a, const threevector<_Tp>& b) {
    return threevector<_Tp>(a) += b;
  }
  
  template<typename _Tp> inline threevector<_Tp> 
  operator-(const threevector<_Tp>& a, const threevector<_Tp>& b) {
    return threevector<_Tp>(a) -= b;
  }
  
  template<typename _Tp> inline threevector<_Tp> 
  operator*(const threevector<_Tp>& a, const _Tp & b) {
    return threevector<_Tp>(a) *= b;
  }

  template<typename _Tp> inline threevector<_Tp> 
  operator*(const _Tp & b, const threevector<_Tp>& a) {
    return threevector<_Tp>(a) *= b;
  }
  
  template<typename _Tp> inline threevector<_Tp> 
  operator/(const threevector<_Tp>& a, const _Tp & b) {
    return threevector<_Tp>(a) /= b;
  }

  template<typename _Tp> inline _Tp 
  operator*(const threevector<_Tp>& a, const threevector<_Tp>& b) {
    return a.X()*b.X() + a.Y()*b.Y() + a.Z()*b.Z();
  }

  template<typename _Tp> inline bool 
  operator==(const threevector<_Tp>& a, const threevector<_Tp>& b) {
    return a.X() == b.X() && a.Y() == b.Y() && a.Z() == b.Z();
  }

  template<typename _Tp> inline bool 
  operator!=(const threevector<_Tp>& a, const threevector<_Tp>& b) {
    return a.X()!= b.X() || a.Y() != b.Y() || a.Z() != b.Z();
  }

  template<typename _Tp> inline _Tp 
  dot(const threevector<_Tp>& a, const threevector<_Tp>& b) {
    return a.X()*b.X() + a.Y()*b.Y() + a.Z()*b.Z();
  }

  template<typename _Tp> inline threevector<_Tp>
  cross(const threevector<_Tp>& a, const threevector<_Tp>& b) 
  {
    return threevector<_Tp>(a.Y()*b.Z() - a.Z()*b.Y(),
						    a.Z()*b.X() - a.X()*b.Z(),
			    			a.X()*b.Y() - a.Y()*b.X());
  }
  
  //        Specializations
  inline float 
  cosAngle(const threevector<float>& a, const threevector<float>& b) 
  {
    float ptot2 = a.mag2()*b.mag2();
    return ptot2 <= 0.0 ? 1.0 : dot(a, b)/std::sqrt(ptot2);
  }
  
  inline double 
  cosAngle(const threevector<double>& a, const threevector<double>& b) 
  {
    double ptot2 = a.mag2()*b.mag2();
    return ptot2 <= 0.0 ? 1.0 : dot(a, b)/std::sqrt(ptot2);
  }

  inline long double 
  cosAngle(const threevector<long double>& a, const threevector<long double>& b)
  {
    long double ptot2 = a.mag2()*b.mag2();
    return ptot2 <= 0.0L ? 1.0L : dot(a, b)/std::sqrt(ptot2);
  }

  
  inline float
  angle(const threevector<float>& a, const threevector<float>& b) 
  {
    double ptot2 = a.mag2()*b.mag2();
    return ptot2 <= 0.0f ? 0.0f : std::acos(dot(a, b)/std::sqrt(ptot2));
  }

  inline double
  angle(const threevector<double>& a, const threevector<double>& b) 
  {
    double ptot2 = a.mag2()*b.mag2();
    return ptot2 <= 0.0 ? 0.0 : std::acos(dot(a, b)/std::sqrt(ptot2));
  }

  inline long double
  angle(const threevector<long double>& a, const threevector<long double>& b) 
  {
    long double ptot2 = a.mag2()*b.mag2();
    return ptot2 <= 0.0L ? 0.0L : std::acos(dot(a, b)/std::sqrt(ptot2));
  }

  //    I/O operations
  template<typename _Tp>
  inline std::istream& operator>>(std::istream& is, threevector<_Tp>& a) {
	char delin;
	is>>delin>>a.X()>>delin>>a.Y()>>delin>>a.Z()>>delin;
	return is;
  }
  
  template<typename _Tp>
  inline std::ostream& operator<<(std::ostream& os, const threevector<_Tp>& q) {
	return os<<"("<<q.X()<<","<<q.Y()<<","<<q.Z()<<")";
  }
}  //  namespace nlo

#endif
