//----------------------------------------------------------------------
// Benchmark 
// Goal: Compare fastJet and SCJet algorithms for the same input files.
// S.Chekanov (ANL)
//----------------------------------------------------------------------
#include "fjcore.hh"
#include "Timer.h"
#include <iostream> 
#include <cstdio>   
#include "ParticleD.h"
#include "KT.h"

using namespace std;

/// an example program showing how to use fastjet
int main(){

  int n=0; 

  // set max number of input particles
  int maxPart=500;

  // read in input particles
  vector<fjcore::PseudoJet> input_particles1;
  vector<ParticleD*> input_particles2;
 
  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
    // create a fjcore::PseudoJet with these components and put it onto
    // back of the input_particles vector
    if (n>maxPart) continue; 
    input_particles1.push_back(fjcore::PseudoJet(px,py,pz,E)); 
    ParticleD *pp = new ParticleD(px,py,pz,E);
    input_particles2.push_back(pp);
    n++;
  }
  
   cout << "Reading " << input_particles1.size() << " particles" << endl;

   Timer tm;
   tm.start();

  // create a jet definition: 
  // a jet algorithm with a given radius parameter
  //----------------------------------------------------------
  double R = 0.6;
  fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, R);
// fjcore::JetDefinition jet_def(fjcore::kt_algorithm, R);

  // run the jet clustering with the above jet definition
  fjcore::ClusterSequence clust_seq(input_particles1, jet_def);


  // get the resulting jets ordered in pt
  double ptmin = 5.0;
  vector<fjcore::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));


  // tell the user what was done
  //  - the description of the algorithm used
  //  - extract the inclusive jets with pt > 5 GeV
  //  show the output as 
  //  {index, rap, phi, pt}
  //----------------------------------------------------------
  cout << "Run " << jet_def.description() << endl;

  tm.stop();
  std::cout << "Final fastjet calculations (ms):" << tm.duration() << std::endl;

  // label the columns
  printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    printf("%5u %15.8f %15.8f %15.8f\n",
	   i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
	   inclusive_jets[i].perp());
  }

   // SCjet implementation
   // build anti-kt jets with R=0.6
   cout << "\nBuilding SCJet for comparison.." << endl; 
   KT* jet= new KT(R, 1, -1, ptmin);
   jet->setDebug(false);
   jet->buildJets(input_particles2);
   vector<ParticleD*> cjets=jet->getJetsSorted(); 
   printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
   for (unsigned int i = 0; i < cjets.size(); i++) {
    printf("%5u %15.8f %15.8f %15.8f\n",
           i, cjets[i]->rapidity(), cjets[i]->phi(),
           cjets[i]->et());
  }

  int xmax=inclusive_jets.size();
  if (cjets.size()<xmax) xmax=cjets.size();

  cout << "\nRelative difference between FastJet and CJet (%)" << endl;
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
      double x1=(double)(100*(inclusive_jets[i].rap()-cjets[i]->rapidity())/inclusive_jets[i].rap());
      double x2=(double)(100*(inclusive_jets[i].phi()-cjets[i]->phi())/inclusive_jets[i].phi());
      double x3=(double)(100*(inclusive_jets[i].perp()-cjets[i]->et())/inclusive_jets[i].perp());
      printf("%5u %15.8f %15.8f  %15.8f \n", i, x1,x2,x3); 
  }
   

  return 0;
}
