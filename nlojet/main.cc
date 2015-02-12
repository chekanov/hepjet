// A kt-jet code taken from NLOjet++ program
// modified to use as anti-kT
// use it as: main < input.txt
// S.Chekanov (ANL)

#include <iostream> 
#include <cstdio>   
#include <kT_clus.h>
#include <vector>
#include <algorithm>
#include "Timer.h"

using namespace std;
using namespace nlo;

typedef lorentzvector<double> Lv;

const double PI2 = 6.28318530716;
const double PI  = 3.14159265358;

struct pT_sort {
    bool operator()(const Lv& p1, const Lv& p2) const {
      return p1.perp2() > p2.perp2();
    }
  };


int main(){
  
  //  private types
  // the jet structore in breit & lab. frame
  bounded_vector<Lv> inclusive_jets;
  bounded_vector<unsigned int> jet;

/*
   p[-1] = lorentzvector<double>(0.0, 0.0, 1.0,1.0);  
   p[0] = lorentzvector<double>(0.0, 0.0, 1.0,1.0);
   p[1] = lorentzvector<double>(0.0, 0.0, 1.0,1.0);
   p[2] = lorentzvector<double>(0.0, 0.0, 10.0,10.0);
   p[3] = lorentzvector<double>(0.0, 0.0, 11.0,11.0);
   p[4] = lorentzvector<double>(0.0, 10.0, 10.0,20.0);
   p[5] = lorentzvector<double>(0.0, 15.0, 11.0,29.0);
*/

  double px, py , pz, E;
  event_hhc p(355);

  // inclomin hadron
  p[hadron(-1)] = lorentzvector<double>(0.0, 0.0, 900.0,900.0);
  p[hadron(0)] = lorentzvector<double>(0.0, 0.0, 900.0,900.0);

  
  int n=1; 
  while (cin >> px >> py >> pz >> E) {
     p[n] = lorentzvector<double>(px,py,pz,E);
     n++;
  }
  cout << "Reading " << n << " particles" << endl; 

/*
  event_hhc::const_iterator p1 = p.begin()+2;
  event_hhc::const_iterator pn = p.end();
  for(event_hhc::const_iterator idir=p1; idir != pn ; ++idir ) {
    cout << idir->X() << " " << idir->Y() << " " << idir->Z() << " " << idir->T() << endl; 
  }
*/

  Timer tm;
  tm.start();

  double ptmin=5.0;

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

// kt-case
//  kT_clus_long jetclus(1,false,1);

// anti-KT case
   kT_clus_long jetclus(-1,false,1);

  jetclus(p,0.6);

  // run 
  jetclus.incl(inclusive_jets, jet);
  std::sort(inclusive_jets.begin(), inclusive_jets.end(), pT_sort());
  tm.stop();

  // label the columns
  printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
  for (unsigned int i = 1; i < inclusive_jets.size(); i++) {
    if (inclusive_jets[i].perp()<ptmin) continue;
    double phi=inclusive_jets[i].phi();
    if (phi<0) phi=PI2+phi;
    printf("%5u %15.8f %15.8f %15.8f\n",
           i, inclusive_jets[i].rapidity(), phi, inclusive_jets[i].perp());
  }

  std::cout << "Final nlojet calculation :" << tm.duration() << " ms " << std::endl;

  return 0;
}
