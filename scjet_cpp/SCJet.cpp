/**
This is class to build kT-type jets using the CSJet clustering algorithm.
The algorithm is released under the GNU public license.
S.Chekanov (ANL) 
**/


#include <iostream>
#include <cstdio>
#include "SCJet.h"
#include "Timer.h"
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace std;
const double PI2 = 2*M_PI;
const double PI  = M_PI;

bool comp(ParticleD* rhs1, ParticleD* rhs2) {
	return rhs1->getPt2() > rhs2->getPt2();
};


/** @brief Initialize calculations of the KT-type of jets
    @param R distance parameter
    @param recom recombination mode. Only recom=1 is supported (E-scheme, p=p1+p2) 
    @param mode defines the algorithm. 1 means KT, -1 is anti-KT, 0 is Cambridge/Aachen
    @param isfast. if true, use N*N algorithm for anti-kT. If false, use the traditional that scales as N^3. 
    @author S.Chekanov 
    @version 1.0 02/02/14
*/
SCJet::SCJet(double R, int recom, int mode, double minpt, bool isfast)
{
	m_R = R;
	m_R2 = (R * R);
	m_recom = recom;
	m_debug = false;
	m_minpt = minpt;
	m_mode = mode;
	m_fast=isfast;

	if (mode == 1)
	{
		std::cout << "SCJet: Longitudinally invariant kT algorithm" << std::endl;
	}
	else if (mode == 0)
	{
		std::cout << "SCjet: Cambridge/Aachen algorithm" << std::endl;
	}
	else if (mode == -1)
	{
		std::cout << "SCJet: Longitudinally invariant anti-kT algorithm" << std::endl;
	}
	else
	{
		std::cout << "SCjet: Not correct mode:  Fallback to the inclusive kT algorithm using E-scheme and R=" << R << std::endl;
	}

	if (m_recom !=1) {
		std::cout << "SCJet: Error. Currently only the E-mode is supported (p1+p2)" <<  std::endl;
		exit(0);
	}

	if (m_fast ==true && m_mode==-1){
		std::cout <<  "SCjet: Fast mode for anti-kT is enabled." << endl;
	}

	if (m_fast ==true && m_mode>=0) {
		std::cout << "SCjet: Currently, the fast mode is enabled for anti-KT jets. Exit." <<  std::endl;
		exit(0);
	}


}


/** @brief Main method to build  KT-type of jets
    @param list list with particles 
    @return list with unsorted jets 
    @author S.Chekanov 
    @version 1.0 02/02/14
*/
std::vector<ParticleD*> SCJet::buildJets(std::vector<ParticleD*> &list)
{

	Timer tm;
	jets = std::vector<ParticleD*>();
	int size = list.size();

	if (m_debug)
		tm.start();

	double min12 = std::numeric_limits<double>::max(); // distance in a pair d_{12}
	double min1  = std::numeric_limits<double>::max(); // distance to the beam d_{iB}

	// *****************************************
	// build a  cache and find first closest pair
	// ******************************************
	double ktdistance1[size];
	int    is_consider[size]; // 0-ignore, 1-original,  n>1 - merged, -1 - jet
	for (int m = 0; m < size; m++) {
		is_consider[m] = 1;
		ParticleD *p1 = list[m];
		ktdistance1[m] = getKtDistance1(p1);
	}


	double ktdistance12[size][size];
	for (int i = 0; i < size - 1; i++) {
		ParticleD *p1 = list[i];
		for (int j = i + 1; j < size; j++) {
			ParticleD *p2 = list[j];
			ktdistance12[i][j] = getKtDistance12(p1, p2);
		}
	}


	if (m_debug) {
		tm.stop();
		std::cout << "After creating a cache of distances (ms):" << tm.duration() << std::endl;
		tm.start();
	}


	int Nstep = size;
	//  cout << "Merge " << Nstep << " j1=" << j1 << " j2=" << j2<< endl;
	bool merged=false;
	int j1 = -1;
	int j2 = -1;
	int km=-1;
	int iter=0;
	int i,j;

	// start loop over all objects
	while (Nstep > 0) {

		min12 = std::numeric_limits<double>::max(); // distance in a pair d_{12}
		min1 =  std::numeric_limits<double>::max(); // distance to the beam d_{iB}

		// this is fast antiKT jet algorithm
		// build pseudo-jet around particles with large pT
		if (m_fast ==true) {

			// find smallest d12.
			// this is after reseting to a new jet
			if (!merged) {

				for (i=0; i < size-1; i++) {
					if (is_consider[i]<=0) continue;
					for (j=i+1; j < size; j++) {
						if (is_consider[j]<=0) continue;
						if (ktdistance12[i][j] < min12) {
							min12 = ktdistance12[i][j];
							j1 = i;
							j2 = j;
						}
					}
				}
			} else {
				// find another minimum around this jet  when j1>0
				for (j=0; j < size; j++) {
					if (is_consider[j]<=0 || j==j1) continue;
					if (ktdistance12[j1][j] < min12) {
						min12 = ktdistance12[j1][j];
						j1 = j1;
						j2 = j;
					}
				}
			} // end of min finding


			// find min distance to the beam
			min1 = ktdistance1[j1];
			if (ktdistance1[j2]<min1) {min1 = ktdistance1[j2];};

			// protect against -1
			if (merged==false && Nstep==1) break;

			// make the decision about this particle
			merged = false;
			if (min12 < min1)  merged = true;



		} else  {  // end fast antiKT
			// start the usual kT algorithm..
			// -----------------------------//


			// find min distance to the beam
			km =-1;
			for (j = 0; j < size; j++) {
				if (is_consider[j]<=0) continue;
				if (ktdistance1[j] < min1) {
					min1 = ktdistance1[j];
					km = j;
				}
			}


			j1=-1;
			j2=-1;
			// find smallest dij distances
			for (i=0; i < size-1; i++) {
				if (is_consider[i]<=0) continue;
				for (j=i+1; j < size; j++) {
					if (is_consider[j]<=0) continue;
					if (ktdistance12[i][j]<min1) {
						min1=ktdistance12[i][j];
						j1 = i;
						j2 = j;
					}
				}
			}


			// make the decision about this particle
			merged=false;
			if (j1>-1 && j2>-1) merged=true;


		} // end standard kT




		// merging ..
		if (merged && j1 != j2) {
			//if (Nstep==1) cout << "Merge " << Nstep << " j1=" << j1 << " j2=" << j2<< endl;
			ParticleD *p1 = list[j1];
			ParticleD *p2 = list[j2];
			p1->add(p2,j2); // p1=p1+p2. Also keeps an index j2
			Nstep--;
			list[j1] = p1;   // replaced with p1=p1+p2
			//list[j2] = NULL; //
			is_consider[j2]=0; // p2, but keep in the list
			is_consider[j1]=is_consider[j1]+1;
			// recalculate distance for this particle
			ktdistance1[j1] = getKtDistance1(p1);
			for (int i = 0; i < size; i++) {
				if (is_consider[i]<=0) continue;
				ParticleD *pp1 = list[i];
				ktdistance12[j1][i] = getKtDistance12(p1, pp1);
				if (m_mode <0) ktdistance12[i][j1] = getKtDistance12(p1, pp1);
			}

		}


		// create a jet
		if (!merged) {   // add this to the jet
			if (!m_fast) j1=km; // this is for KT and C/A
			is_consider[j1] = -1;
			ParticleD *pj = list[j1];
			Nstep--;
			if (pj->getPt() > m_minpt)
			{
				jets.push_back(pj); // fill jets
			}
		}


		if (m_debug) {
			cout << "## Iteration:" << iter++ << endl;
			for (int i = 0; i < size; i++) {
				ParticleD *p1 = list[i];
				std::string mess="original";
				if (is_consider[i]==-1) mess="!final-jet!";
				if (is_consider[i]>1) mess="(proto-jet)";
				if (is_consider[i]==0) mess="(removed)";
				cout << i << "  E=" << p1->e() << " " << mess << endl;
			}
		}



	} // end loop


	if (m_debug) cout << "Final Nr of iterations=" << iter << endl;

	// attempt to deal with un-mearged particle
	int ins=-1;
	for (int i = 0; i < size; i++)
	if (is_consider[i]==1) {ins=i;};

	if (ins>-1) {
		//cout << "working on the last un-merged particle " << endl;
		ParticleD *p2 = list[ins];
		if (m_debug) cout << "Unmerged particle id=" << ins << " y=" << p2->getRapidity() << " phi=" <<  p2->getPhi() << " pt=" <<  p2->getPt() << endl;
		min12 = std::numeric_limits<double>::max();
		for (unsigned int j = 0; j < jets.size(); j++) {
			ParticleD *lp = jets[j];
			double d=getDistance(p2, lp);
			if (d<min12) { j1=j; min12 =d; };
		}
		if (m_debug) cout << "Distance R to closest jet=" << min12 << endl;
		if (min12<m_R) {
			if (m_debug) cout << " --> Particle merged" << endl;
			ParticleD *lp = jets[j1];
			lp->add(p2,j1);
			is_consider[ins] = 0;
		}
	}

	// sanity test. All particles were merged?
	if (m_debug) {
		int nn=0; int ins=-1;
		for (int i = 0; i < size; i++)
		if (is_consider[i]==1) {nn++; ins=i;};
		if (nn != 0)   std::cout << "--> WARNING: particle with ID=" << ins<< " un-merged" << endl;
	}


	if (m_debug) {
		tm.stop();
		std::cout << "  --> Final time for calculation (ms):" << tm.duration() << std::endl;
		printJets();
	}

	return jets;

}



/** @brief Get sorted jets
    @return list with sorted jets 
    @author S.Chekanov 
    @version 1.0 02/02/14
*/
std::vector<ParticleD*> SCJet::getJetsSorted() {
	sort(jets.begin(), jets.end(), comp);
	return jets;
}


// get deltaR distance
double SCJet::getDistance(ParticleD *a, ParticleD *b)
{
        double rsq, deltaEta, deltaPhi;
        deltaEta = b->getRapidity() - a->getRapidity();
        double phi1 = a->getPhi();
        double phi2 = b->getPhi();
        deltaPhi = phi2 - phi1;
        if (deltaPhi>PI) deltaPhi=PI2-deltaPhi;
        if (deltaPhi<-PI) deltaPhi=PI2+deltaPhi;
        rsq = (deltaEta*deltaEta + deltaPhi*deltaPhi);
        return sqrt(rsq);
}


/** Print sorted jets
*/
void SCJet::printJets() {

	std::vector<ParticleD*> sjets = getJetsSorted();

	std::cout << "# Nr of jets=" << sjets.size() << std::endl;
	printf("%5s %15s %15s %15s %8s\n","jet #", "rapidity", "phi", "pt", " const");
	for (unsigned int i = 0; i < sjets.size(); i++)
	{
		ParticleD *lp = sjets[i];
		double phi = lp->phi();
		std::vector<int> con = lp->getConstituents();
		if (phi < 0) phi = PI + phi;
		double srap=lp->getRapidity();
		double pt=lp->getPt();
		printf("%5u %15.8f %15.8f %15.8f %8d\n",
		       i, srap, phi,  pt, (int)con.size());
	}


}


// calculate KT distance between any 2 particles
double SCJet::getKtDistance12(ParticleD *a, ParticleD *b)
{
	double rsq, esq, deltaEta, deltaPhi;
	deltaEta = b->getRapidity() - a->getRapidity();
	double phi1 = a->getPhi();
	double phi2 = b->getPhi();
	deltaPhi = phi2 - phi1;
	if (deltaPhi>PI) deltaPhi=PI2-deltaPhi;
	if (deltaPhi<-PI) deltaPhi=PI2+deltaPhi;

	rsq = (deltaEta*deltaEta + deltaPhi*deltaPhi);
	esq = 0;
	if (m_mode == 1)
	{
		esq = std::min(a->getPt2(), b->getPt2()); // kT
	}
	else if (m_mode == 0)
	{
		esq = 1.0 ;  // C-A
	}
	else if (m_mode == -1)
	{
		esq = std::min(1.0 / a->getPt2(), 1.0 / b->getPt2()); // anti-KT
	}
	else
	{
		esq = std::min(a->getPt2(), b->getPt2()); // kT (fallback)
	}

	return (esq * rsq / m_R2);
}



// calculate distance to the beam
double SCJet::getKtDistance1(ParticleD *a)
{
	if (m_mode == 1)
	{
		return a->getPt2();
	}
	else if (m_mode == 0)
	{
		return 1.0;
	}
	else if (m_mode == -1)
	{
		return (1.0 / a->getPt2());
	}
	return a->getPt2();
}


// set debugging mode
void SCJet::setDebug(bool debug)
{

	m_debug = debug;
}

