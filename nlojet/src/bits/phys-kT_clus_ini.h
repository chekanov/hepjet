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
#ifndef __NLO_PHYS_KT_CLUS_INI_H__
#define __NLO_PHYS_KT_CLUS_INI_H__ 1

//   Standard includes
#include <iostream>
#include <lorentzvector.h>
#include <bits/hep-bounded_vector.h>
#include <bits/hep-bounded_matrix.h>


namespace nlo {
  
  /**
   *   Abstract class for implementing the longitudinary boost invariant
   *   \f$k_\perp\f$ algorithm.
   *
   *
   *
   */
  class kT_clus_ini
  {
    //   private types
    typedef lorentzvector<double> _Lv;
    
  public:
      
    ///  Do cluster analysis.
    ///  Eg.:  
    ///  ktclus_hhc durham;
    ///  event_hhc event(10);
    ///  const bounded_vector<double>& y;
    ///
    ///  y = durham(p.begin()+2, p.end());       //  ecut = etot, R = 1.0
    ///  y = durham(ecut, p.begin()+2, p.end(), R);
    ///
    /// ecut = input  : denominator of kt measure
    /// R    = input  : cone size related parameter (Ellis-Soper)
    /// y[j] = output : value of y for which event changes from being
    ///                 j jet to j-1 jet
    /// \param first ecut = input  : denominator of kt measure
    /// \param last ecut = input  : denominator of kt measure
    /// \return value of y for which event changes from being
    ///                 j jet to j-1 jet
    template<class _InputIterator> const bounded_vector<double>& 
    operator()(_InputIterator first, _InputIterator last, double ecut = 1.0);
    
    template<class _InputIterator> const bounded_vector<double>& 
    operator()(double, _InputIterator, _InputIterator, double = 1.0);
    
    ///   output scales at which jets were merged 
    double operator[](unsigned int i) const { return _M_y[i];}
    
    ///  Count the number of jets at each value of ycut, for event which has
    ///  already been analysed by operator()(..). 
    ///
    ///      njet[0] = durham.ysub(ecut, ycut[0]); 
    ///  or
    ///      durham.ycut(ecut, ny, ycut, njet);
    ///       
    ///  ecut    = input  : denominator of kt measure
    ///  ny      = input  : number of ycut values
    ///  ycut[j] = input  : y values at which numbers of jets are counted
    ///  njet[j] = output : number of jets at ycut[j]
    unsigned int ycut(double, double) const;
    void ycut(double, unsigned int, double *, unsigned int *) const;
    
    ///  Count the number of sub-jets at each value of ycut, for event which
    ///  has already been analysed by operator()(..). Remember that a sub-jet
    ///  is defined as a jet at y=ycut which has not yet been merged with the 
    ///  beam at y=ymac.
    ///
    ///      nsub[0] = durham.ysub(ecut, ycut[0], ymac); 
    ///  or
    ///      durham.ysub(ecut, ny, ycut, ymac, nsub);
    ///
    ///    ecut    = input : denominator of kt measure
    ///    ny      = input : number of ycut values
    ///    ycut[j] = input : y values at which numbers of sub-jets are counted
    ///    ymac    = input : y value used to define macro-jets, to determine
    ///                      which jets are sub-jets
    ///    nsub[j] = output : number of sub-jets at ycut[j]
    unsigned int ysub(double, double, double) const;
    void ysub(double, unsigned int, double *, double, unsigned int *) const;
    
    //  Give same information as last call to operator()(..) except that only
    //  transitions where a jet was merged with the beam jet are recorded.
    //
    //     durham.beam(ecut, y); 
    //
    //   ecut = input  : denominator of kt measure
    //   y[j] = output : y value where jth hardest jet was merged with beam
    void beam(double, bounded_vector<double>&) const;
    
    //  Give same information as last call to operator()(..) except    
    //  that only transitions where two sub-jets were joined are recorded
    //  remember that a sub-jet is defined as a jet at y=ycut which has not
    //  yet been merged with the beam at y=ymac. 
    //
    //    durham.join(ecut, ymac, y);
    //
    //    ecut = input  : denominator of kt measure
    //    ymac = input  : value of y used to define macro-jets
    //    y[j] = output : y value where event changed from having
    //                    n+j sub-jets to having n+j-1, where n is
    //                    the number of macro-jets at scale ymac
    void join(double, double, bounded_vector<double>&) const;

    // Reconstruct kinematics of jet system, which has already been 
    // analysed by operator()(..). 
    //
    //    durham.reco(ecut, ycut, ymac, pjet, jet, njet, nsub);
    //
    //  ecut = input : denominator of kt measure
    //  ycut = input : y value at which to reconstruct jet momenta
    //  ymac = input : y value used to define macro-jets, to determine
    //                   which jets are sub-jets
    //  pjet[j] = output : 4-momentum of jth jet at scale ycut
    //   jet[j] = output : the macro-jet which contains the jth jet,
    //                     set to zero if jet is not a sub-jet
    //  njet    = output : the number of jets
    //  nsub    = output : the number of sub-jets (equal to the number of
    //                     non-zero entries in jet[])
    void reco(double, double, double, bounded_vector<_Lv>&, 
	      bounded_vector<unsigned int>&, unsigned int&, 
	      unsigned int&) const;

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

    template<class _OutputIterator> 
    void reco(_OutputIterator, _OutputIterator, bounded_vector<unsigned int>&) const;
    
    void reco(unsigned int njet, bounded_vector<_Lv>& pj) const {
      this -> reco(pj.begin(), pj.begin()+njet);
    }

    void reco(unsigned int njet, bounded_vector<_Lv>& pj,
	      bounded_vector<unsigned int>& jet) const {
      this -> reco(pj.begin(), pj.begin()+njet, jet);
    }
    
    // Reconstruct kinematics of jet system, which has already been
    // analysed by operator()(..) according to the inclusive jet definition. 
    //
    //    durham.incl(pjet, jet);
    //
    // pjet[j] = output : 4-momentum of jth jet at scale ycut
    //  jet[j] = output : the jet which contains the jth particle
    void incl(bounded_vector<_Lv>&, bounded_vector<unsigned int>&) const;
    
    //  Count the number of sub-jets in the nth inclusive jet of an event
    //  that has already been analysed operator()(..).
    //
    //    durham.isub(n, ny, ycut, nsub);
    //    
    //  n       = input : which inclusive jet to use
    //  ny      = input : number of ycut values
    //  ycut[j] = input : y values at which numbers of sub-jets are counted
    //  nsub[j] =output : number of sub-jets at ycut[j]
    void isub(unsigned int, unsigned int, double *, unsigned int *) const;
    
    //    destructor
    virtual ~kT_clus_ini() {}
    
  private:
    //    private members
    void _M_ktclus(double);
    void _M_ktreco(unsigned int) const;
    void _M_ktreco(unsigned int, bounded_vector<unsigned int>&) const;
    void _M_ktmove(unsigned int, unsigned int) const;

    //    virtual members
    virtual double _M_ktsing(unsigned int) const = 0;    
    virtual double _M_ktpair(unsigned int, unsigned int, double&) const = 0;

    virtual _Lv  _M_ktmom(unsigned int) const = 0;
    virtual void _M_ktcopy(const bounded_vector<_Lv>&) const = 0;
    virtual void _M_ktpmerg(unsigned int, unsigned int) const = 0;
    virtual void _M_ktpmove(unsigned int, unsigned int) const = 0;
    virtual void _M_ktmerg(unsigned int, unsigned int, unsigned int) const = 0;

  protected: 
    //   cone size related parameter (Ellis-Soper)
    double _M_rsq;

    //   momenta of the particles
    bounded_vector<_Lv> _M_pp;
    
    //   jet event shapes and the smallest resolution variables
    bounded_vector<double> _M_kt, _M_y;
    
    //  ktlast records for each merging, the highest Ecut^2 for 
    //  which the result is not merged with the beam
    bounded_vector<double> _M_ktl;

    //   record the merging history
    bounded_vector<unsigned int> _M_hist;
    
    //          temporary variables
    //   temporary variables to store the single
    //   and pair resolution parameters 
    mutable bounded_matrix<double> _M_ktp;
    
    //   temporary variable recording the jet structure 
    mutable bounded_vector<unsigned int> _M_injet;
  };
  
  template<class _InputIterator> inline 
  const bounded_vector<double>& kT_clus_ini::
  operator()(_InputIterator first, _InputIterator last, double R)
  {
    double ecut = 0.0;
    _M_pp.assign(first, last, 1);
    _M_rsq = R*R;

    for(int i = 1; i <= _M_pp.upper(); i++)
      ecut += _M_pp[i].T();
    
    this -> _M_ktclus(ecut);
    return _M_y;
  }
  
  template<class _InputIterator> inline 
  const bounded_vector<double>& kT_clus_ini::
  operator()(double ecut, _InputIterator first, _InputIterator last, double R)
  {
    _M_rsq = R*R;
    _M_pp.assign(first, last, 1);
    this -> _M_ktclus(ecut);
    return _M_y;
  }

  template<class _OutputIterator> void kT_clus_ini::
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
      *iter = this -> _M_ktmom(i++);	
  }

  template<class _OutputIterator> void kT_clus_ini::
  reco(_OutputIterator first, _OutputIterator last,
       bounded_vector<unsigned int>& jet) const
  {
    unsigned int njet = 0U;
    for(_OutputIterator iter = first; iter < last; iter++) ++njet;
    
    if(_M_pp.upper() < (int) njet) {
      std::cerr<<"unable to resolve "<<njet<<" jets"<<std::endl;
      return;
    }
    
    this -> _M_ktreco(njet, jet);
    
    unsigned int i = 1U;
    for(_OutputIterator iter = first; iter < last; iter++) 
      *iter = this -> _M_ktmom(i++);	
  }
} //  namespace nlo


#endif
