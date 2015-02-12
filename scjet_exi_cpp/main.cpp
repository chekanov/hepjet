// An example program to show how to run scjet (C++)
// using an input 
// use it as: main < input.txt
// S.Chekanov (ANL)


#include <iostream> 
#include <cstdio>   
#include "ParticleD.h"  
#include "KT.h"   
#include "Timer.h"

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

  Timer tm;
  tm.start();

  // build anti- jets with R=0.6
  KT* jet= new KT(0.6, 1, -1, 5.0);
  jet->setDebug(false);
  jet->buildJets(input_particles);
  vector<ParticleD*> cjets=jet->getJetsSorted();
  tm.stop();

  cout << "Final output:" << endl;
  printf("%5s %15s %15s %15s %8s\n","jet #", "rapidity", "phi", "pt", " const");
  for (unsigned int i = 0; i < cjets.size(); i++) {
     ParticleD *lp = cjets[i];
     //double phi = lp->phi();
    std::vector<int> con = lp->getConstituents();
    printf("%5u %15.8f %15.8f %15.8f %8d\n",
           i, cjets[i]->rapidity(), cjets[i]->phi(),
           cjets[i]->et(), (int)con.size());
  }

   std::cout << "Final SCjet calculation :" << tm.duration() << " ms " << std::endl;

  return 0;
}
