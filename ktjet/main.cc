// An example program to show how to run scjet (C++)
// using an input 
// use it as: main < input.txt
// S.Chekanov (ANL)


#include <iostream> 
#include <cstdio>   

/** Need to include these KtJet Headers */
#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"
#include "Timer.h"


using namespace std;
using namespace KtJet;
const double PI2=6.2831853071;
const double PI=3.1415926535;


int main(){
  
  // read in input particles
  std::vector<KtJet::KtLorentzVector> avec;

  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
    KtJet::KtLorentzVector p(px,py,pz,E);
    avec.push_back(p);
  }

   Timer tm;
   tm.start();


   cout << "Reading " << avec.size() << " particles" << endl; 

   /** set KtEvent flags */
  int type  = 4; // PP
//  int angle = 2; // deltaR 
  int angle = -2; // deltaR, but anti-kT 
  int recom = 1; // E
  double rparameter = 0.6;
  double ETmin=5.0;


  /** Construct the KtEvent object */
  KtJet::KtEvent ev(avec,type,angle,recom,rparameter);

  /** Print out the number of final state jets */
  std::cout << "Number of final state jets: " << ev.getNJets() << std::endl;

  /** Retrieve the final state jets from KtEvent sorted by Et*/
  std::vector<KtLorentzVector> jets = ev.getJetsEt();

  /** Print out jets 4-momentum and Pt */
  std::vector<KtLorentzVector>::const_iterator itr = jets.begin();

/*
  for( ; itr != jets.end() ; ++itr) {
      double et=(*itr).perp();
      if (et<ETmin) continue;
    std::cout << "Jets 4 vector: " 
	      << (*itr).px() << " " 
	      << (*itr).py() << " " 
	      << (*itr).pz() << " " 
	      << (*itr).e() << std::endl;
    std::cout << "Jets Pt: " << (*itr).perp() << std::endl; 
  }
*/

  tm.stop();

  // label the columns
  printf("%5s %15s %15s %15s %7s\n","jet #", "rapidity", "phi", "pt", "const");

  // print out the details for each jet
  for (unsigned int i = 0; i < jets.size(); i++) {
    double et=jets[i].perp();
    int cons=jets[i].getNConstituents();
    if (et<ETmin) continue;
    double phi=jets[i].phi();
    if (phi<0) phi=PI2+phi;
    printf("%5u %15.8f %15.8f %15.8f %7d\n",
           i, jets[i].rapidity(), phi,
           jets[i].perp(),cons);
  }

  std::cout << "Final KtJet++ calculation :" << tm.duration() << " ms" << std::endl;

  return 0;
}
