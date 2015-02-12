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

//   nlo includes
#include "threevector.h"


namespace nlo {

  
#define __THREEVECTOR_MEMBERS(_Tp)					\
  void threevector<_Tp>::                                               \
  rotate(const value_type& psi, const threevector<_Tp>& axis)   	\
  {                                                                     \
    value_type ph = axis.phi(), th = axis.theta();      		\
									\
    rotateZ(-ph);							\
    rotateY(-th);							\
    rotateZ(psi);							\
    rotateY(th);							\
    rotateZ(ph);							\
  }									\
									\
  void threevector<_Tp>::rotateX(const value_type& psi)                 \
  {                                                     		\
    value_type sp = std::sin(psi), cp = std::cos(psi);			\
    value_type py = _M_y, pz = _M_z;		        		\
									\
    _M_y = cp*py - sp*pz;         					\
    _M_z = sp*py + cp*pz;        					\
  }									\
									\
  void threevector<_Tp>::rotateY(const value_type& th)                  \
  {                                                     		\
    value_type st = std::sin(th), ct = std::cos(th);       		\
    value_type px = _M_x, pz = _M_z;	        			\
									\
    _M_x =  ct*px + st*pz;				         	\
    _M_z = -st*px + ct*pz;			        		\
  }									\
									\
  void threevector<_Tp>::rotateZ(const value_type& ph)                  \
  {                                                          		\
    value_type cp = std::cos(ph), sp = std::sin(ph);			\
    value_type px = _M_x, py = _M_y;    				\
									\
    _M_x = cp*px - sp*py;       					\
    _M_y = sp*px + cp*py;	            				\
  }

__THREEVECTOR_MEMBERS(float)
__THREEVECTOR_MEMBERS(double)
__THREEVECTOR_MEMBERS(long double)

}   //  namespace nlo
