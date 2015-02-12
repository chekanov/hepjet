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
#ifndef __NLO_PHYS_KT_CLUS_EPA_H__
#define __NLO_PHYS_KT_CLUS_EPA_H__ 1


//   Standard includes
#include <iostream>

#include <bits/epa-event.h>
#include <bits/hep-bounded_vector.h>
#include <bits/hep-bounded_matrix.h>


namespace nlo {

  
  class kT_clus_epa
  {
    //   private types
    typedef lorentzvector<double> _Lv;
    
  public:
    //  Do cluster analysis.
    //  Eg.:  
    //     ktclus_dis durham;
    //     event_dis event(10);
    //     const bounded_vector<double>& y;
    //     
    //     y = durham(p.begin()+2, p.end());       //  ecut = etot 
    //     y = durham(p.begin()+2, p.end(), ecut);
    //
    //  The first argument refers to the momentum of the incoming parton.
    //
    //  or
    //     y = durham(p);                          //  ecut = etot 
    //     y = durham(p, ecut);
    //     
    //    ecut = input  : denominator of kt measure
    //    y[j] = output : value of y for which event changes from being
    //                    j jet to j-1 jet
    template<class _InputIterator> const bounded_vector<double>& 
    operator()(_InputIterator, _InputIterator);
   
    template<class _InputIterator> const bounded_vector<double>& 
    operator()(_InputIterator, _InputIterator, double);
    
    const bounded_vector<double>& operator()(const event_epa&);
    const bounded_vector<double>& operator()(const event_epa&, double);
    
    //   output scales at which jets were merged 
    double operator[](unsigned int i) const { return _M_kt[i];}
    
    // Reconstruct kinematics of n-jet system, which has already been 
    // analysed by operator()(..). 
    //
    //     bounded_vector<_Lv> pj(1, njet);
    //     durham.reco(pj.begin(), pj.begin() + njet);
    // or
    //     durham.reco(njet, pj);
    //
    // njet  = input  : the number of the jet
    // pj[j] = output : 4-momentum of jth jet at scale ycut
    template<class _OutputIterator> 
    void reco(_OutputIterator, _OutputIterator) const;
    
    void reco(unsigned int njet, bounded_vector<_Lv>& pj) const {
      this -> reco(pj.begin(), pj.begin()+njet);
    }

  private:
    //    private members
    void   _M_ktclus(double);
    void   _M_ktreco(unsigned int) const;
    void   _M_ktmerg(unsigned int, unsigned int, unsigned int) const;
    void   _M_ktmove(unsigned int, unsigned int) const;
    double _M_ktpair(unsigned int, unsigned int) const;
    
    //   momenta of the particles
    bounded_vector<_Lv> _M_pp;
    
    //   jet event shapes and the smallest resolution variables
    bounded_vector<double> _M_kt;
    
    //   record the merging history
    bounded_vector<unsigned int> _M_hist;
    
    //     temporary variables
    //   temporary member to strore the momenta
    mutable bounded_vector<_Lv> _M_p;
    
    //   temporary variables to store the pair resolution parameters 
    mutable bounded_matrix<double, symmetric> _M_ktp;
  };

  
  template<class _InputIterator> inline 
  const bounded_vector<double>& kT_clus_epa::
  operator()(_InputIterator first, _InputIterator last)
  {
    double ecut = 0.0;
    _M_pp.assign(first, last, 1);
    
    for(_InputIterator iter = first; iter < last; iter++)
      ecut += iter -> T();
    
    _M_ktclus(ecut);
    return _M_kt;
  }
  
  template<class _InputIterator> inline 
  const bounded_vector<double>& kT_clus_epa::
  operator()(_InputIterator first, _InputIterator last, double ecut)
  {
    _M_pp.assign(first, last, 1);
    _M_ktclus(ecut);
    return _M_kt;
  }
    
  inline const bounded_vector<double>& 
  kT_clus_epa::operator()(const event_epa& p) {
    return this -> operator()(p.begin()+2, p.end());
  }
 
  inline const bounded_vector<double>& 
  kT_clus_epa::operator()(const event_epa& p, double ecut) {
    return this -> operator()(p.begin()+2, p.end(), ecut);
  }
  
  template<class _OutputIterator> void kT_clus_epa::
  reco(_OutputIterator first, _OutputIterator last) const
  {
    unsigned int njet = 0U;
    for(_OutputIterator iter = first; iter < last; iter++) ++njet;
    
    if(_M_pp.upper() < (int) njet) {
      std::cerr<<"unable to resolve "<<njet<<" jets"<<std::endl;
      return;
    }
    
    this -> _M_ktreco(njet);
  
    unsigned int i = 1U;
    for(_OutputIterator iter = first; iter < last; iter++) 
      *iter = _M_p[i++];	
  }

  inline double kT_clus_epa::_M_ktpair(unsigned int i, unsigned int j) const
  {
    double E = (_M_p[i].T() < _M_p[j].T() ? _M_p[i].T() : _M_p[j].T());
    return 2.0*E*E*(1.0 - cosAngle(_M_p[i], _M_p[j]));
  }
} //  namespace nlo


#endif
