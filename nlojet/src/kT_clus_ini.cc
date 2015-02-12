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

#include "bits/phys-kT_clus_ini.h"


using namespace std;

namespace nlo {

  
  void kT_clus_ini::_M_ktclus(double ecut)
  {

    cout << "# Distance parameter R=" << sqrt(_M_rsq) << endl;

    //  number of particles
    unsigned int nt = _M_pp.upper();
    
    //  initialize
    this -> _M_ktcopy(_M_pp);
    
    _M_y.resize(1, (int) nt + 1);
    _M_kt.resize(1, (int) nt + 1);
    _M_ktl.resize(1, (int) nt);
    _M_ktp.resize(1, (int) nt);
    _M_hist.resize(1, (int) nt);
    _M_injet.resize(1, (int) nt);
    _M_kt[nt+1] = _M_y[nt+1] = 0.0;
    
    //----- calculate all single and pair kt's -----
    unsigned int i, j;
    for(i = 1; i <= nt; i++) {
      _M_ktp[i][i] = this -> _M_ktsing(i);
      for(j = i+1; j <= nt; j++) {
	_M_ktp[j][i] = -1.0;
	_M_ktp[i][j] = this -> _M_ktpair(i, j, _M_ktp[j][i]);
      }
    }
    
    //----- main loop -----
    bool merge;
    unsigned int imin = 1, jmin = 2, kmin = 1;
    double ktp, ktp_min, kts_min, ktmin, ktmax = 0.0;
    double etsq = ecut*ecut;
    const double eps = 1.0e-10;
    
    for(unsigned n = nt; n > 0; n--) {
      //----- find minimum member of ktp and kts -----
      ktp_min = kts_min = 9.9e123*etsq;
      for(i = 1; i <= n; i++) {
	if((ktp = _M_ktp[i][i]) < kts_min) {
	  kts_min = ktp;
	  kmin = i;
	}
	
	for(j = i+1; j <= n; j++)
	  if((ktp = _M_ktp[i][j]) < ktp_min) {
	    ktp_min = ktp;
	    imin = i; jmin = j;
	  }
      }
      
      kts_min *= _M_rsq; ktmin = ktp_min; merge = true;
      if(kts_min < eps || kts_min <= ktp_min) {
	ktmin = kts_min;
	merge = false;
      }
      
      if(ktmin > ktmax) ktmax = ktmin;
      
      //----- smallest resolution variables ----- 
      _M_y[n] = ktmin/etsq;
      _M_kt[n] = ktmin;
      
      if(!merge) {
	_M_ktmove(kmin, n);
	_M_hist[n] = _M_injet[n] = kmin;
	
	for(i = n; i <= nt; i++)
	  if(_M_injet[i] == kmin) {
	    _M_ktl[i] = ktmax;
	    _M_injet[i] = 0;
	  } else if(_M_injet[i] == n) _M_injet[i] = kmin;
      } else {
	_M_ktmerg(n, imin, jmin);
	_M_ktmove(jmin, n);
	_M_hist[n] = imin*(nt+1) + jmin;
	_M_injet[n] = imin;
	
	for(i = n; i <= nt; i++)
	  if(_M_injet[i] == jmin) _M_injet[i] = imin;
	  else if(_M_injet[i] == n) _M_injet[i] = jmin;
      }
    }
  }
  
  void kT_clus_ini::_M_ktmove(unsigned int j, unsigned int n) const
  {
    unsigned int i;
    this -> _M_ktpmove(j, n);
    _M_ktp[j][j] = _M_ktp[n][n];
    
    for(i = 1; i < j; i++) {
      _M_ktp[i][j] = _M_ktp[i][n];
      _M_ktp[j][i] = _M_ktp[n][i];
    }
    
    for(i = j+1; i < n; i++) {
      _M_ktp[j][i] = _M_ktp[i][n];
      _M_ktp[i][j] = _M_ktp[n][i];
    } 
  } 
  
  void  kT_clus_ini::
  ycut(double ecut, unsigned int ny, double *yc, unsigned int *njet) const
  {
    unsigned int i, j, nn = (unsigned int) _M_pp.upper();
    double etsq = 1.0/(ecut*ecut);
    
    for(i = 0; i < ny; i++) 
      njet[i] = 0;
    
    for(i = nn; i >= 1; i--)
      for(j = 0; j < ny; j++)
	if(njet[j] == 0 && _M_kt[i]*etsq >= yc[j]*0.99999) njet[j] = i;
  }
  
  unsigned int kT_clus_ini::ycut(double ecut, double yc) const
  {
    unsigned int nn = (unsigned int) _M_pp.upper();
    double etsq = 1.0/(ecut*ecut);
    
    for(unsigned int i = nn; i >= 1; i--)
      if(_M_kt[i]*etsq >= yc*0.99999) return i;
    
    return 0;
  }
  
  void kT_clus_ini::
  ysub(double ecut, unsigned int ny, double *ycut, double ymac, unsigned int *nsub) const
  {
    unsigned int i, j, nn = (unsigned int) _M_pp.upper();
    double etsq = 1.0/(ecut*ecut);
    
    for(i= 0; i < ny; i++) 
      nsub[i] = 0;
    
    for(i = nn; i >= 1; i--)
      for(j = 0; j < ny; j++) {
	if(nsub[j] == 0 && _M_kt[i]*etsq >= ycut[j]*0.99999) nsub[j] = i;
	if(nsub[j] != 0 && _M_ktl[i]*etsq < ymac*0.99999) --nsub[j];
      }
  }
  
  unsigned int kT_clus_ini::ysub(double ecut, double ycut, double ymac) const
  {
    unsigned int nsub = 0, nn = (unsigned int) _M_pp.upper();
    double etsq = 1.0/(ecut*ecut);
    
    for(unsigned int i = nn; i >= 1; i--) {
      if(nsub == 0 && _M_kt[i]*etsq >= ycut*0.99999) nsub = i;
      if(nsub != 0 && _M_ktl[i]*etsq < ymac*0.99999) --nsub;
    }
    
    return nsub;
  }
  
  void kT_clus_ini::beam(double ecut, bounded_vector<double>& y) const
  {
    unsigned int i, j = 1, nn = (unsigned int) _M_pp.upper();
    double etsq = 1.0/(ecut*ecut);
    y.resize(1, (int) nn);
    
    for(i = 1; i <= nn; i++)
      if(_M_hist[i] <= nn + 1) 
	y[j++] = etsq*_M_kt[i];
    
    for(i = j; i <= nn; i++) 
      y[i] = 0.0;
  }

  void kT_clus_ini::join(double ecut, double ymac, bounded_vector<double>& y) const
  {
    unsigned int i, j = 1, nn = (unsigned int) _M_pp.upper();
    double etsq = 1.0/(ecut*ecut);
    y.resize(1, (int) nn);
    
    for(i = 1; i <= nn; i++)
      if(_M_hist[i] > nn + 1 && etsq*_M_ktl[i] >= ymac*0.99999)
	y[j++] = etsq*_M_kt[i];
    
    for(i = j; i <= nn; i++) 
      y[i] = 0.0;
  }
  
  void kT_clus_ini::
  reco(double ecut, double ycut, double ymac, 
       bounded_vector<_Lv>& pjet, bounded_vector<unsigned int>& jet, 
       unsigned int& njet, unsigned int& nsub) const
  {
    unsigned int i, imin, jmin, nn = (unsigned int) _M_pp.upper();
    double etsq = 1.0/(ecut*ecut);
    
    //----- copy the momenta -----
    this -> _M_ktcopy(_M_pp);
    
    //----- keep merging until ycut ----
    unsigned int n = nn;
    while(etsq*_M_kt[n] < ycut*0.99999 && n > 0) {
      if(_M_hist[n] <= nn+1) 
	this -> _M_ktpmove(_M_hist[n], n);
      else {
	imin = _M_hist[n]/(nn+1);
	jmin = _M_hist[n] - (nn+1)*imin;
	this -> _M_ktpmerg(imin, jmin);
	this -> _M_ktpmove(jmin, n);
      }
      --n;
    }
    
    //----- if ycut is too large there are no jets -----
    njet = nsub = n;
    if(n == 0) return;
    
    //----- set up output momenta -----
    jet.resize(1, (int) njet);
    pjet.resize(1, (int) njet);
    for(i = 1; i <= njet; i++) {
      jet[i] = i;
      pjet[i] = this -> _M_ktmom(i);
    }
    
    //----- keep merging until ymac to find the fate of each jet-----
    while(etsq*_M_kt[n] < ymac*0.99999 && n > 0) {
      if(_M_hist[n] <= nn+1) {
	imin = 0;
	jmin = _M_hist[n];
	--nsub;
      } else {
	imin = _M_hist[n]/(nn+1);
	jmin = _M_hist[n] - (nn+1)*imin;
	if(etsq*_M_ktl[n] < ymac*0.99999) --nsub;
      }
      
      for(i = 1; i <= njet; i++) {
	if(jet[i] == jmin) jet[i] = imin;
	if(jet[i] == n) jet[i] = jmin;
      }
      --n;
    }
  }
  
  void kT_clus_ini::_M_ktreco(unsigned int njet) const
  {
    unsigned int imin, jmin, nt = (unsigned int) _M_pp.upper();
    unsigned int n = nt;
    
    //----- copy the momenta -----
    this -> _M_ktcopy(_M_pp);

    while(n > njet) {
      if(_M_hist[n] <= nt+1) 
	this -> _M_ktpmove(_M_hist[n], n);
      else {
	imin = _M_hist[n]/(nt+1);
	jmin = _M_hist[n] - (nt+1)*imin;
	this -> _M_ktpmerg(imin, jmin);
	this -> _M_ktpmove(jmin, n);
      }
      --n;
    }
  }

  void kT_clus_ini::
  _M_ktreco(unsigned int njet, bounded_vector<unsigned int>& jet) const
  {
    unsigned int imin, jmin, nt = (unsigned int) _M_pp.upper();
    unsigned int i, n = nt;
    
    //----- copy the momenta -----
    this -> _M_ktcopy(_M_pp);

    //----- set up output momenta -----
    jet.resize(1, (int) nt);
    for(i = 1; i <= nt; i++)
      jet[i] = i;

    while(n > njet) {
      if(_M_hist[n] <= nt+1) {
	imin = 0;
	jmin = _M_hist[n];
	this -> _M_ktpmove(_M_hist[n], n);
      } else {
	imin = _M_hist[n]/(nt+1);
	jmin = _M_hist[n] - (nt+1)*imin;
	this -> _M_ktpmerg(imin, jmin);
	this -> _M_ktpmove(jmin, n);
      }

      for(i = 1; i <= nt; i++) {
	if(jet[i] == jmin) jet[i] = imin;
	if(jet[i] == n) jet[i] = jmin;
      }

      --n;
    }
  }

  void kT_clus_ini::
  incl(bounded_vector<_Lv>& pjet, bounded_vector<unsigned int>& jet) const
  {
    unsigned int i, imin, jmin, nn = (unsigned int) _M_pp.upper();
    
    //----- copy the momenta -----
    this -> _M_ktcopy(_M_pp);
    
    //----- set up output momenta -----
    jet.resize(1, (int) nn);
    for(i = 1; i <= nn; i++)
      jet[i] = i;
    
    //----- keep merging to the bitter end -----
    unsigned int njet = 0U;
    pjet.resize(1, 0);
    
    for(unsigned int n = nn; n > 0U; n--){ 
      if(_M_hist[n] <= nn+1) {
	imin = 0;
	jmin = _M_hist[n];
	njet++;
	pjet.push_back(this -> _M_ktmom(jmin));
	this -> _M_ktpmove(jmin, n);
      } else {
	imin = _M_hist[n]/(nn+1);
	jmin = _M_hist[n] - imin*(nn+1);
	this -> _M_ktpmerg(imin, jmin);
	this -> _M_ktpmove(jmin, n);
      }
      
      for(i = 1; i <= nn; i++) {
	if(jet[i] == jmin) jet[i] = imin;
	if(jet[i] == n) jet[i] = jmin;
	if(jet[i] == 0) jet[i] = nn+1 + njet;
      }
    }
    
    //----- finally every particle must be in an inclusive jet -----
    for(i = 1; i <= nn; i++) jet[i] -= nn+1;
  }
  
  void kT_clus_ini::
  isub(unsigned int n, unsigned int ny, double *ycut, unsigned int *nsub) const
  {
    unsigned int i, nn = (unsigned int) _M_pp.upper();

    for(i = 0; i < ny; i++) 
      nsub[i] = 0U;
 
    //----- find which merging corresponds to the nth inclusive jet -----
    unsigned int nm = 0, j = 0;
    for(i = nn; i > 0; i--) {
      if(_M_hist[i] <= nn+1) j++;
      if(j == n) {
	nm = i;
	break;
      }
    }

    //----- give up if there are less than n inclusive jets
    if(nm == 0) {
      std::cerr<<"kT_clus_ini::isub : there are less than"
	       <<n<<"inclusive jets"<<std::endl;
      return;
    }
    
    for(i = nn; i > 0; i--)
      for(j = 0; j < ny; i++) {
	if(nsub[j] == 0 && _M_rsq*_M_kt[i] >= 0.99999*ycut[j]*_M_kt[nm]) nsub[j] = i;
	if(nsub[j] == 0 && std::abs(_M_ktl[i] - _M_ktl[nm]) > 1e-6) nsub[j] = nsub[j]-1;
      }
  }



}
















