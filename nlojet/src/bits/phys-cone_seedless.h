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
#ifndef __NLO_PHYS_CONE_SEEDLESS_H__
#define __NLO_PHYS_CONE_SEEDLESS_H__ 1

//  Standard C++ includes
#include <list>

//   nlo++ includes
#include <lorentzvector.h>

#include <bits/dis-event.h>
#include <bits/hhc-event.h>
#include <bits/hhc2ph-event.h>
#include <bits/hep-bounded_vector.h>


namespace nlo {

  class cone_seedless
  {
    //   private types
    typedef lorentzvector<double> _Lv;
    
    struct _Vector {
      _Vector() 
	: p(0.0,0.0,0.0,0.0), pt(0.0), eta(0.0), phi(0.0) {}
      
      _Vector(const _Lv& q) 
	: p(q), pt(q.perp()), eta(q.rapidity()), phi(q.phi()) {}
      
      lorentzvector<double> p;
      double pt, eta, phi;
    };

    struct _Proto {
      _Vector psum;
      std::list<unsigned int> pars;
    };
    
  public:
    explicit cone_seedless(int reco = 0, int midpoint = 0) 
      : _M_reco(reco), _M_midpoint(midpoint) {}
    
    template<class _InputIterator> const bounded_vector<_Lv>& 
    operator()(_InputIterator first, _InputIterator last, double R=0.7, double f=0.5)
    {
      _M_pp.assign(first, last, 1);
      this -> _M_do_clustering(R, f);
      return _M_pj;
    }

    const bounded_vector<_Lv>& 
    operator()(const event_dis& p, double R=0.7, double f=0.5) {
      return this -> operator()(p.begin()+3, p.end(), R, f);
    }
    
    const bounded_vector<_Lv>& 
    operator()(const event_hhc& p, double R=0.7, double f=0.5) {
      return this -> operator()(p.begin()+2, p.end(), R, f);
    }
    
    const bounded_vector<_Lv>& 
    operator()(const event_hhc2ph& p, double R=0.7, double f=0.5) {
      return this -> operator()(p.begin()+4, p.end(), R, f);
    }
    
  private:
    void _M_do_clustering(double, double);
    void _M_check_trial_cone(unsigned int, unsigned int *, double);
    unsigned int _M_check_shared_towers(const _Proto&, const _Proto&);
    void _M_split_merge(double, std::list<_Proto>::iterator, std::list<_Proto>::iterator);
    
    void _M_assadd(_Vector&, const _Vector&) const;
    void _M_iterate_cone(double, const _Vector&);
    void _M_was_it_already_found(const _Proto *);

    //       private data members
    //   recombination scheme
    int _M_reco, _M_midpoint;
    
    //   momenta of the particles
    bounded_vector<_Vector> _M_pp;
    
    //   list of the protojets
    std::list<_Proto> _M_proto; 
    
    //   the list of the overlapping partons
    std::list<unsigned int> _M_overlap; 
    
    //   temporary variable to store the jet momenta
    mutable bounded_vector<_Lv> _M_pj; 
    
    
    //   private static members
    static double _S_dphi(double);
    
    struct pT_sort {
      bool operator()(const _Proto& p1, const _Proto& p2) const {
	return p1.psum.pt > p2.psum.pt;
      }
    };
  };

}


#endif
