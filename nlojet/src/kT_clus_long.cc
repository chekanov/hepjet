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

//using namespace std;


namespace nlo {
  
     
  void kT_clus_long::_M_ktpmove(unsigned int j, unsigned int n) const {
    _M_p[j] = _M_p[n];
  } 
  
  lorentzvector<double> kT_clus_long::_M_ktmom(unsigned int i) const 
  {
    if(_M_reco == 1) return _M_p[i].p;
    else {
      double pT = _M_p[i].pt, ei = _M_p[i].eta, fi = _M_p[i].phi;
      return pT*_Lv(std::cos(fi), std::sin(fi), std::sinh(ei), std::cosh(ei)); 
    }
  }
  
  void kT_clus_long::_M_ktcopy(const bounded_vector<_Lv>& p) const 
  {
    int nt = p.upper();
    
    _M_p.resize(1, nt);
    for(int i = 1; i <= nt; i++) {
      if(_M_reco == 1) _M_p[i].p = p[i];
      _M_p[i].pt = p[i].perp();
      _M_p[i].eta = p[i].rapidity();
      _M_p[i].phi = p[i].phi();
    }
  }
  
  double kT_clus_long::_M_ktsing(unsigned int i) const 
  {

    double pi;
    // kT
   if (_M_angle !=-1) {
        pi = _M_p[i].pt;
   }
    // antikT
    if (_M_angle==-1) {
        pi = 1.0/_M_p[i].pt;  
    }

    if (_M_angle==0) return pi; // C/A case
 
    return pi*pi;
  }
  
  double kT_clus_long::
  _M_ktpair(unsigned int i, unsigned int j, double& ang) const
  {

    double pi,pj,pT;


    // kT case
    if (_M_angle !=-1) { 
         pi = _M_p[i].pt, pj = _M_p[j].pt;
         pT = (pi < pj ? pi : pj);
    }


    // anti-kT case
    if (_M_angle==-1) {
        pi = 1.0/_M_p[i].pt, pj = 1.0/_M_p[j].pt;
        pT = (pi < pj ? pi : pj);
   }
    
    double deta = _M_p[i].eta - _M_p[j].eta;
    double dphi = _M_ktdphi(_M_p[i].phi - _M_p[j].phi);
    
    if(_M_angle == 1) ang = deta*deta+dphi*dphi;
    else ang = 2.0*(std::cosh(deta) - std::cos(dphi));
 
    // C/A case 
    if (_M_angle==0) return pT*ang;
 
    // either kT and anti-kT
    return pT*pT*ang;

  }
  
  double kT_clus_long::_M_ktdphi(double phi) const
  {
    static const double pi = 3.14159265358979323846;
    static const double twopi = 6.28318530717958647692;
    
    if(phi >= pi) phi = std::fmod(pi+phi, twopi) - pi;
    else if(phi < -pi) phi = -std::fmod(pi-phi, twopi) + pi;
    
    return phi;
  }

  void kT_clus_long::_M_ktpmerg(unsigned int i, unsigned int j) const 
  {
    //--- combine the two momenta ---
    switch(_M_reco) {
      //--- E recombination scheme ---
    case 1: _M_p[i].p += _M_p[j].p; break;
    case 2: case 3: 
      {
        //--- pT or pT^2 weighted schemes ---
        double pi = _M_p[i].pt, pj = _M_p[j].pt;
        double ei = _M_p[i].eta, ej = _M_p[j].eta; 
        double fi = _M_p[i].phi, fj = _M_p[j].phi;
        
        //--- weighted sum ---
        _M_p[i].pt = pi + pj;
        if(_M_reco == 3) { pi *= pi; pj *= pj;} 
        _M_p[i].eta = (pi*ei + pj*ej)/(pi+pj);
	_M_p[i].phi = _M_ktdphi(fi + pj*_M_ktdphi(fj-fi)/(pi+pj));
      }
      break;
    }
  }


  void kT_clus_long::
  _M_ktmerg(unsigned int n, unsigned int i, unsigned int j) const
  {
    //--- if the scheme is monotonic then combine the relative angles ---
    if(_M_mono == true) {
      double pi, pj;
      unsigned int ii, jj, ik, jk;
      
      for(unsigned int k = 1; k <= n; k++) {
        if(_M_reco == 1) { pi = _M_p[i].p.T(); pj = _M_p[j].p.T();} 
	else { pi = _M_p[i].pt; pj = _M_p[j].pt;}
	
	if(_M_reco == 3) { pi *= pi; pj *= pj;}
        if(pi == 0.0 && pj == 0.0) pi = pj = 1.0;
	
        if((ii = i) < (ik = k)) std::swap(ii, ik);
        if((jj = j) < (jk = k)) std::swap(jj, jk);
        _M_ktp[ii][ik] = (pi*_M_ktp[ii][ik] + pj*_M_ktp[jj][jk])/(pi+pj);
      }
    }
    
    //--- combine the two momenta ---
    this -> _M_ktpmerg(i, j);
    if(_M_reco == 1) {
      _M_p[i].pt = _M_p[i].p.perp();
      _M_p[i].eta = _M_p[i].p.rapidity();
      _M_p[i].phi = _M_p[i].p.phi();
    }
    
    //--- recalculate the corresponding kT values ---
    _M_ktp[i][i] = this -> _M_ktsing(i);
    
    double tmp, pT = _M_p[i].pt;
    unsigned int ii, ik;
    for(unsigned int k = 1; k <= n; k++)
      if(k != i && k != j) {
        if((ii = i) > (ik = k)) std::swap(ii, ik);
        if(_M_mono == true) {
          tmp = _M_p[k].pt;
          _M_ktp[ii][ik] = (pT < tmp ? pT*pT : tmp*tmp)*_M_ktp[ik][ii];
        } else _M_ktp[ii][ik] = _M_ktpair(i, k, tmp);
      }
  }
}
