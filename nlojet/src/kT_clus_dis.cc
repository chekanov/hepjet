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
#include "kT_clus.h"



namespace nlo {

  
  lorentzvector<double> kT_clus_dis::_M_ktmom(unsigned int i) const {
    return _M_p[i];
  }
  
  void kT_clus_dis::_M_ktpmove(unsigned int j, unsigned int n) const {
    _M_p[j] = _M_p[n];
  } 
  
  void kT_clus_dis::
  _M_ktmerg(unsigned int n, unsigned int i, unsigned int j) const
  {
    _M_p[i] += _M_p[j];
    _M_ktp[i][i] = this -> _M_ktsing(i);
    
    double tmp;
    unsigned int ii, kk;
    for(unsigned int k = 1; k <= n; k++)
      if(k != i && k != j) {
	if((ii = i) > (kk = k)) std::swap(ii, kk);
	_M_ktp[ii][kk] = this -> _M_ktpair(i, k, tmp);
      }  
  }
  
  void kT_clus_dis::_M_ktpmerg(unsigned int i, unsigned int j) const {
    _M_p[i] += _M_p[j];
  }
  
  void kT_clus_dis::_M_ktcopy(const bounded_vector<_Lv>& p) const {
    _M_p = p;
  }
  
  double kT_clus_dis::_M_ktsing(unsigned int i) const 
  {
    double E = _M_p[i].T();
    return 2.0*E*E*(1.0 - cosAngle(_M_p[i], _M_n));
  }
  
  double kT_clus_dis::
  _M_ktpair(unsigned int i, unsigned int j, double&) const
  {
    double E = (_M_p[i].T() < _M_p[j].T() ? _M_p[i].T() : _M_p[j].T());
    double angle = 1.0 - cosAngle(_M_p[i], _M_p[j]);
    return 2.0*E*E*angle;
  }
}
