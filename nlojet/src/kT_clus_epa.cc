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
#include "bits/phys-kT_clus_epa.h"


namespace nlo {
  

  void kT_clus_epa::_M_ktclus(double ecut)
  {
    //  number of particles
    unsigned int nt = _M_pp.upper();
    
    //  initialize
    _M_kt.resize(2, (int) nt + 1);
    _M_hist.resize(2, (int) nt);
    _M_ktp.resize(1, nt);
    
    _M_p = _M_pp;
    _M_kt[nt+1] = 0.0;  
    
    //----- calculate all single and pair kt's -----
    unsigned int i, j; 
    
    for(i = 1; i < nt; i++)
      for(j = i+1; j <= nt; j++)
	_M_ktp[i][j] = _M_ktpair(i, j);
    
    // ----- main loop -----
    unsigned int imin = 1, jmin = 2;
    double ktp, ktmin, etsq = ecut*ecut;
    
    for(unsigned n = nt; n > 1; n--) {
      //----- find minimum member of ktp and kts -----
      ktmin = 9.9e123*etsq;
      for(i = 1; i < n; i++) 
	for(j = i+1; j <= n; j++) {
	  ktp = _M_ktp[i][j];
	  if(ktp < ktmin) {
	    ktmin = ktp;
	    imin = i; jmin = j;
	  }
	}
      
      //----- smallest resolution variables ----- 
      _M_kt[n] = ktmin/etsq;
      
      //----- do merge -----
      _M_ktmerg(n, imin, jmin);
      _M_ktmove(n, jmin);
      _M_hist[n] = imin*(nt + 1) + jmin;
    }
  }
  
  void kT_clus_epa::_M_ktreco(unsigned int njet) const
  {
    unsigned int imin, jmin, nt = (unsigned int) _M_pp.upper();
        
    //----- copy the momenta -----
    _M_p = _M_pp;
    
    //----- keep merging until njet ----
    unsigned int n = nt;
    while(n > njet) {
      imin = _M_hist[n]/(nt + 1);
      jmin = _M_hist[n] - (nt + 1)*imin;
      _M_p[imin] += _M_p[jmin];
      _M_p[jmin]  = _M_p[n];
      --n;
    }
  }

  void kT_clus_epa::_M_ktmove(unsigned int n, unsigned int j) const
  {
    unsigned int i;
    _M_p[j] = _M_p[n];
    
    for(i = 1; i < j; i++)
      _M_ktp[i][j] = _M_ktp[i][n];
    
    for(i = j+1; i < n; i++)
      _M_ktp[j][i] = _M_ktp[i][n];
  } 
  
  
  void kT_clus_epa::
  _M_ktmerg(unsigned int n, unsigned int i, unsigned int j) const
  {
    _M_p[i] += _M_p[j];
    
    for(unsigned int k = 1; k <= n; k++)
      if(k != i && k != j)
	_M_ktp[i][k] = _M_ktpair(i, k);
  }
}   //  namespace nlo
