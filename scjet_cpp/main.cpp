#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include "ParticleD.h"   // needed for io
#include "KT.h"   // needed for io

using namespace std;

/// an example program showing how to use fastjet
int main(){
  
  // read in input particles
  //----------------------------------------------------------

    // read in input particles
  //----------------------------------------------------------
  vector<ParticleD*> input_particles;

  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
      ParticleD *pp = new ParticleD(px,py,pz,E);
      input_particles.push_back(pp);
  }

   cout << "Reading " << input_particles.size() << " particles" << endl; 
   KT* jet= new KT(0.6, 1, -1, 5.0);
   jet->setDebug(true);
   jet->buildJets(input_particles);

  return 0;
}
