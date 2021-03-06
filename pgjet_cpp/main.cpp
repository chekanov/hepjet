// An example program to show how to run scjet (C++)
// using an input 
// use it as: main < input.txt
// S.Chekanov (ANL)


#include <iostream>
#include <iomanip>
#include <list>
#include <set>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include "ParticleD.h"  
#include "KTpg.h"   
#include "Timer.h"


using namespace std;

int main(){
  
  // read in input particles
  vector<ParticleD> input_particles;

  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
      ParticleD pp(px,py,pz,E);
      input_particles.push_back(pp);
  }

   cout << "Reading " << input_particles.size() << " particles" << endl;

  Timer tm;
  tm.start();

  // build anti- jets with R=0.6
  KTpg* kt= new KTpg(0.6, 1, -1, 5.0);
  kt->setDebug(false);
  kt->buildJets(input_particles);
  vector<ParticleD> cjets=kt->getJetsSorted();
  tm.stop();

  cout << "Final output:" << endl;
  printf("%5s %15s %15s %15s %8s\n","jet #", "rapidity", "phi", "pt", " const");
  for (unsigned int i = 0; i < cjets.size(); i++) {
    ParticleD lp = cjets[i];
    double phi = lp.phi();
   if (phi < 0) phi = PI2 + phi; 
    std::vector<ParticleD> con = lp.getConstituents();
    printf("%5u %15.8f %15.8f %15.8f %8d\n",
           i, cjets[i].rapidity(), phi,
           cjets[i].perp(), (int)con.size());
  }

   std::cout << "Final PGjet calculation :" << tm.duration() << " ms " << std::endl;

  return 0;
}
