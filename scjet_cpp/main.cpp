// An example program to show how to run scjet (C++)
// using an input 
// use it as: main < input.txt
// S.Chekanov (ANL)


#include <iostream> 
#include <cstdio>   
#include "ParticleD.h"  
#include "KT.h"   

using namespace std;

int main(){
  
  // read in input particles
  vector<ParticleD*> input_particles;

  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
      ParticleD *pp = new ParticleD(px,py,pz,E);
      input_particles.push_back(pp);
  }

   cout << "Reading " << input_particles.size() << " particles" << endl; 
   // build anti-kt jets with R=0.6
   KT* jet= new KT(0.6, 1, -1, 5.0);
   jet->setDebug(true);
   jet->buildJets(input_particles);

  return 0;
}
