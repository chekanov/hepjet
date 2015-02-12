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
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

//   nlo includes
#include "lorentzvector.h"


namespace nlo {


#define __LORENTZVECTOR_MEMBERS_SPINOR(_Tp)				       \
  void lorentzvector<_Tp>::boost(const _Tp& bx, const _Tp& by, const _Tp& bz)  \
  {									       \
    _Tp b2 = bx*bx + by*by + bz*bz;					       \
    _Tp gamma = 1.0/std::sqrt(1.0 - b2);				       \
    _Tp bp = bx*_M_x + by*_M_y + bz*_M_z;				       \
    _Tp gamma2 = (b2 > 0.0 ? (gamma - 1.0)/b2 : 0.0);			       \
									       \
    _M_x = _M_x + gamma2*bp*bx + gamma*bx*_M_t;				       \
    _M_y = _M_y + gamma2*bp*by + gamma*by*_M_t;				       \
    _M_z = _M_z + gamma2*bp*bz + gamma*bz*_M_t;				       \
    _M_t = gamma*(_M_t + bp);						       \
  }									       \
									       \
  std::complex<_Tp> scalar_mp(lorentzvector<_Tp> a, lorentzvector<_Tp> b)      \
  {									       \
    std::complex<_Tp> I(0.0,1.0), ea(1.0), eb(1.0), ph(1.0);		       \
    _Tp pa = a.perp(), pb = b.perp();					       \
									       \
    if(a.T() < 0.0) a *= -1.0, ph *= I;					       \
    if(b.T() < 0.0) b *= -1.0, ph *= I;					       \
    									       \
    if(pa != 0.0) ea = std::complex<_Tp>(a.X(), a.Y())/pa;		       \
    if(pb != 0.0) eb = std::complex<_Tp>(b.X(), b.Y())/pb;		       \
    _Tp mp = a.minus()*b.plus();					       \
    _Tp pm = a.plus()*b.minus();					       \
    mp = (mp <= 0.0 ? 0.0 : mp);					       \
    pm = (pm <= 0.0 ? 0.0 : pm);					       \
        								       \
    return ph*(std::sqrt(mp)*ea - std::sqrt(pm)*eb);			       \
  } \
									       \
  lorentzvector<std::complex<_Tp> > vector_pp(lorentzvector<_Tp> a, lorentzvector<_Tp> b) \
  {																						  \
    typedef lorentzvector<std::complex<_Tp> > return_type;\
    std::complex<_Tp> I(0.0,1.0), ph(1.0), ea(1.0), eb(1.0);  \
    _Tp pa = a.perp(), pb = b.perp(); \
	 \
    if(a.T() < 0.0) a *= -1.0, ph *= I; \
    if(b.T() < 0.0) b *= -1.0, ph *= I; \
     \
	if(pa != 0.0) ea = std::complex<_Tp>(a.X(), a.Y())/pa; \
	if(pb != 0.0) eb = std::complex<_Tp>(b.X(), b.Y())/pb; \
 \
    _Tp ap = a.plus(), am = a.minus(), bp = b.plus(), bm = b.minus(); \
	ap = (ap <= 0.0 ? 0.0 : ap); bp = (bp <= 0.0 ? 0.0 : bp); \
	am = (am <= 0.0 ? 0.0 : am); bm = (bm <= 0.0 ? 0.0 : bm); \
 \
    std::complex<_Tp> t1 = ph*std::sqrt(bp*ap); \
    std::complex<_Tp> t2 = ph*std::sqrt(bm*am)*eb*conj(ea); \
    std::complex<_Tp> t3 = ph*std::sqrt(bp*am)*conj(ea); \
    std::complex<_Tp> t4 = ph*std::sqrt(bm*ap)*eb; \
     \
    return return_type(t3 + t4, I*(t3 - t4), t1 - t2, t1 + t2); \
  } 
  
  
  __LORENTZVECTOR_MEMBERS_SPINOR(float)
  __LORENTZVECTOR_MEMBERS_SPINOR(double)
  __LORENTZVECTOR_MEMBERS_SPINOR(long double)

}   //  namespace nlo
