//----------------------------------------------------------------------
/// \file
/// \page Example01 01 - basic usage example
///
/// fastjet basic example program:
///   simplest illustration of the usage of the basic classes:
///   fjcore::PseudoJet, fjcore::JetDefinition and 
///   fjcore::ClusterSequence
///
/// run it with    : ./01-basic < single-event.dat
///
/// Source code: 01-basic.cc
//----------------------------------------------------------------------

//STARTHEADER
// $Id: 01-basic.cc 2684 2011-11-14 07:41:44Z soyez $
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
#include "fjcore.hh"
#include "Timer.h"


#include <iostream> // needed for io
#include <cstdio>   // needed for io

using namespace std;

/// an example program showing how to use fastjet
int main(){
 

  int n=0; 
  // read in input particles
  //----------------------------------------------------------
  vector<fjcore::PseudoJet> input_particles;
  
  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
    // create a fjcore::PseudoJet with these components and put it onto
    // back of the input_particles vector
    // cout << n++ << endl;
    input_particles.push_back(fjcore::PseudoJet(px,py,pz,E)); 
  }
  



   Timer tm;
   tm.start();

  // create a jet definition: 
  // a jet algorithm with a given radius parameter
  //----------------------------------------------------------
  double R = 0.6;
 fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, R);
// fjcore::JetDefinition jet_def(fjcore::kt_algorithm, R);


  // run the jet clustering with the above jet definition
  //----------------------------------------------------------
  fjcore::ClusterSequence clust_seq(input_particles, jet_def);


  // get the resulting jets ordered in pt
  //----------------------------------------------------------
  double ptmin = 5.0;
  vector<fjcore::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));


  // tell the user what was done
  //  - the description of the algorithm used
  //  - extract the inclusive jets with pt > 5 GeV
  //    show the output as 
  //      {index, rap, phi, pt}
  //----------------------------------------------------------
  cout << "Ran " << jet_def.description() << endl;

  tm.stop();

  // label the columns
  printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt",  "n const");
 
  // print out the details for each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    vector<fjcore::PseudoJet> constituents = inclusive_jets[i].constituents();
    printf("%5u %15.8f %15.8f %15.8f %8u\n",
	   i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
	   inclusive_jets[i].perp(),constituents.size());
  }

/*
  vector<fjcore::PseudoJet> cons = inclusive_jets[4].constituents();

  cout << endl;
  for (unsigned int i = 0; i < cons.size(); i++) {
  //  printf("%5u %15.8f %15.8f %15.8f\n",
  //         i, cons[i].rap(), cons[i].phi(),
  //         cons[i].perp());
      printf("%15.8f %15.8f %15.8f %15.8f\n",
           i, cons[i].px(), cons[i].py(),
           cons[i].pz(),cons[i].e());

  }
*/

  std::cout << "Final fastjet calculation :" << tm.duration() << " ms " << std::endl;
  return 0;
}
