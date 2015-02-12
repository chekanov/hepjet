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
#ifndef __NLO_KT_CLUS_H__
#define __NLO_KT_CLUS_H__ 1

//   nlo includes
#include <bits/dis-event.h>
#include <bits/hhc-event.h>
#include <bits/hhc2ph-event.h>

#include <bits/phys-kT_clus_ini.h>
#include <bits/phys-kT_clus_epa.h>



namespace nlo {

  
  class kT_clus_dis : public kT_clus_ini
  {
  public:
    //   constructor
    explicit kT_clus_dis(bool is_Pz_positive = true) 
      : _M_n(0.0, 0.0, (is_Pz_positive ? 1.0 : -1.0)) {}
    
    explicit kT_clus_dis(const threevector<double>& pa) 
      : _M_n(pa) {}
    
    const bounded_vector<double>& 
    operator()(const event_dis& p, double R = 1.0) {
      _M_n = p[0];
      return this -> kT_clus_ini::operator()(p.begin()+3, p.end(), R);
    }
    
    const bounded_vector<double>& 
    operator()(double ecut, const event_dis& p, double R = 1.0) {
      _M_n = p[0];
      return this -> kT_clus_ini::operator()(ecut, p.begin()+3, p.end(), R);
    }
    
  private: 
    //   private members
    typedef lorentzvector<double> _Lv;

    double _M_ktsing(unsigned int) const;    
    double _M_ktpair(unsigned int, unsigned int, double&) const;
    
    _Lv  _M_ktmom(unsigned int) const;
    void _M_ktcopy(const bounded_vector<_Lv>&) const;
    void _M_ktpmerg(unsigned int, unsigned int) const;
    void _M_ktpmove(unsigned int, unsigned int) const;
    void _M_ktmerg(unsigned int, unsigned int, unsigned int) const;
    
    //  direction of the incoming parton
    threevector<double> _M_n;
    
    //   temporary variable to store the jet momenta
    mutable bounded_vector<_Lv> _M_p; 
  };
  
  
  class kT_clus_long : public kT_clus_ini
  {
  public:
    //    constructor
    kT_clus_long()
      : _M_mono(false), _M_reco(1), _M_angle(1) {}
    
    // angle: 1     => DeltaR
    //        2     => f(DeltaEta,DeltaPhi) where f()=2(cosh(eta)-cos(phi))
    //                 is the QCD emission metric
    //        0     => DeltaR Cambridge/Achim case 
    //       -1     => DeltaR (anti-kT case)
    // mono:  false => derive relative pseudoparticle angles from jets
    //        true  => monotonic definitions of relative angles
    // reco:  1     => E recombination scheme
    //        2     => pt scheme
    //        3     => pt**2 scheme
    kT_clus_long(unsigned int angle, bool mono, unsigned int reco)
      : _M_mono(mono), _M_reco(reco), _M_angle(angle) {

     int mode=_M_angle;

     if (mode == 1)
        {
                std::cout << "# NLOjet: Longitudinally invariant kt algorithm" << std::endl;
        }
        else if (mode == 0)
        {
                std::cout << "# Nlojet: Cambridge/Aachen algorithm" << std::endl;
        }
        else if (mode == 2)
        {
                std::cout << "# NLOjet: f(DeltaEta,DeltaPhi) where f()=2(cosh(eta)-cos(phi))" << std::endl;
        }
        else if (mode == -1)
        {
                std::cout << "# NLOjet: Longitudinally invariant anti-kt algorithm" << std::endl;
        }
        else
        {
                std::cout << "# NLOjet: Not correct mode:  Fallback to the inclusive kT algorithm using E-scheme " << std::endl;
        }


    }
    
    void set_up(unsigned int angle, bool mono, unsigned int reco) {
      _M_mono = mono; _M_reco = reco; _M_angle = angle;
    }
    
    //   operators to analize the event
    const bounded_vector<double>& 
    operator()(const event_dis& p, double R = 1.0) {
      return this -> kT_clus_ini::operator()(p.begin()+3, p.end(), R);
    }
    
    const bounded_vector<double>& 
    operator()(double ecut, const event_dis& p, double R = 1.0) {
      return this -> kT_clus_ini::operator()(ecut, p.begin()+3, p.end(), R);
    }

    const bounded_vector<double>& 
    operator()(const event_hhc& p, double R = 1.0) {
      return this -> kT_clus_ini::operator()(p.begin()+2, p.end(), R);
    }
    
    const bounded_vector<double>& 
    operator()(double ecut, const event_hhc& p, double R = 1.0) {
      return this -> kT_clus_ini::operator()(ecut, p.begin()+2, p.end(), R);
    }
    
    const bounded_vector<double>& 
    operator()(const event_hhc2ph& p, double R = 1.0) {
      return this -> kT_clus_ini::operator()(p.begin()+4, p.end(), R);
    }
    
    const bounded_vector<double>& 
    operator()(double ecut, const event_hhc2ph& p, double R = 1.0) {
      std::cout <<  "-> Distance R=" << R << std::endl;
      return this -> kT_clus_ini::operator()(ecut, p.begin()+4, p.end(), R);
    }
    
  private: 
    //   private members
    typedef lorentzvector<double> _Lv;

    double _M_ktdphi(double) const;
    double _M_ktsing(unsigned int) const;    
    double _M_ktpair(unsigned int, unsigned int, double&) const;
    
    _Lv  _M_ktmom(unsigned int) const;
    void _M_ktcopy(const bounded_vector<_Lv>&) const;
    
    void _M_ktpmerg(unsigned int, unsigned int) const;
    void _M_ktpmove(unsigned int, unsigned int) const;
    void _M_ktmerg(unsigned int, unsigned int, unsigned int) const;
    

    struct _Vector {
      lorentzvector<double> p;
      double pt, eta, phi;
    };
    
    bool _M_mono;
    int _M_reco, _M_angle;
    
    //   temporary variable to store the jet momenta
    mutable bounded_vector<_Vector> _M_p; 
  };


}   //   namespace nlo

#endif
