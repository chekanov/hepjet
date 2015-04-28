/**
This is class to build kT-type jets using the CSJet clustering algorithm.
The algorithm is released under the GNU public license.
authors: Ivan Pogrebnyak and Sergei Chekanov  
**/


#include <iostream>
#include <cstdio>
#include "KTpg.h"
#include "Timer.h"
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace std;
//const double PI2 = 2*M_PI;
//const double PI  = M_PI;

bool comp(ParticleD rhs1, ParticleD rhs2) {
	return rhs1.getPt2() > rhs2.getPt2();
};


/** @brief Initialize calculations of the KT-type of jets
    @param R distance parameter
    @param recom recombination mode. Only recom=1 is supported (E-scheme, p=p1+p2) 
    @param mode defines the algorithm. 1 means KT, -1 is anti-KT, 0 is Cambridge/Aachen
    @author S.Chekanov 
    @version 1.0 02/02/14
*/
KTpg::KTpg(double R, int recom, int mode, double minpt)
{
	m_R = R;
	m_R2 = (R * R);
	m_recom = recom;
	m_debug = false;
	m_minpt = minpt;
	m_mode = mode;
	if (mode == 1)
	{
		std::cout << "Longitudinally invariant kt algorithm" << std::endl;
	}
	else if (mode == 0)
	{
		std::cout << "Cambridge/Aachen algorithm" << std::endl;
	}
	else if (mode == -1)
	{
		std::cout << "Longitudinally invariant anti-kt algorithm" << std::endl;
	}
	else
	{
		std::cout << "Not correct mode:  Fallback to the inclusive kT algorithm using E-scheme and R=" << R << std::endl;
	}



}


/** @brief Main method to build  KT-type of jets
    @param list list with particles 
    @return list with unsorted jets 
    @author S.Chekanov 
    @version 1.0 02/02/14
*/

std::vector<ParticleD> KTpg::buildJets(const std::vector<ParticleD> &list) { 


         jets = std::vector<ParticleD>();

         set<ParticleD> particles;
         for (unsigned int i = 0; i < list.size(); i++)
        {
                ParticleD pp = list[i];
                ParticleD pj(pp.px(),pp.py(),pp.pz(),pp.e(),m_mode); 
                particles.insert(pj);
         }       

    if (m_debug) {
       cout << "Input particles and distances after sorting" << endl;
       for(set<ParticleD>::iterator it = particles.begin(); it != particles.end(); it++)
       {
        cout << it->getD() << endl;
       }
    }



   while (particles.size()) cluster(particles, jets);

 //   cout << jets.size() << endl;

   return jets;

    }



void KTpg::cluster(set<ParticleD>& particles, std::vector<ParticleD>& jets) { 

  typedef set<ParticleD>::iterator iter_t;

  iter_t it = particles.begin();
  iter_t end = particles.end(), it1 = end, it2 = end;
  bool merged = false;
  double dist = it->getD()*m_R2;
  double merged_R = 0.;
  for (;it!=end;++it) {
    for (iter_t jt=next(it);jt!=end;++jt) {
        double d = it->distance(*jt);
        if (d < dist) {
          dist = d;
          merged_R = d;
          it1 = it;
          it2 = jt;
          if (!merged) merged = true;
        }
    }
  }


   if (merged) {

    ParticleD isum=(*it1)+(*it2);
    /*
    std::vector<ParticleD> p1=it1->getConstituents(); 
    std::vector<ParticleD> p2=it2->getConstituents();

    isum.setConstituent((*it1));
    isum.setConstituent((*it2));
    isum.setConstituents(p1);
    isum.setConstituents(p2);
    */
    particles.insert(isum);
    particles.erase(it1);
    particles.erase(it2);

    if (m_debug) {
    cout << "Merged "
         << setw(2) << it1->getID() << " & "
         << setw(2) << it2->getID() << "  d = ";
    const ios_base::fmtflags flags = cout.flags();
    cout << fixed << setprecision(8) << dist
         << " R = " << merged_R << endl;
    cout.flags(flags);
   }

  } else {
    iter_t itj = particles.begin();
    double pt = itj->getPt();
    if (pt> m_minpt) jets.push_back(*particles.begin()); 
    particles.erase(particles.begin());
    }

}



/** @brief Get sorted jets
    @return list with sorted jets 
    @author S.Chekanov 
    @version 1.0 02/02/14
*/
std::vector<ParticleD> KTpg::getJetsSorted() {
	sort(jets.begin(), jets.end(), comp);
	return jets;
}


/** Print sorted jets
*/
void KTpg::printJets() {

	std::vector<ParticleD> sjets = getJetsSorted();

	std::cout << "# Nr of jets=" << sjets.size() << std::endl;
	printf("%5s %15s %15s %15s %8s\n","jet #", "rapidity", "phi", "pt", " const");
	for (unsigned int i = 0; i < sjets.size(); i++)
	{
		ParticleD lp = sjets[i];
		double phi = lp.phi();
		std::vector<ParticleD> con = lp.getConstituents();
		if (phi < 0) phi = PI + phi;
		double srap=lp.getRapidity();
		double pt=lp.getPt();
		printf("%5u %15.8f %15.8f %15.8f %8d\n",
		       i, srap, phi,  pt, (int)con.size());
	}


}



// set debugging mode
void KTpg::setDebug(bool debug)
{

	m_debug = debug;
}

