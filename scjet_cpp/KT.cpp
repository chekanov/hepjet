// KT is class to build KT-type of jets in inclusive mode
// S.Chekanov (ANL)

#include "KT.h"
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
    @author S.Chekanov 
    @version 1.0 02/02/14
*/
KT::KT(double R, int recom, int mode, double minpt)
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
std::vector<ParticleD*> KT::buildJets(std::vector<ParticleD*> &list)
{

        Timer tm;
        jets = std::vector<ParticleD*>();
	int size = list.size();

	if (m_debug)
	     tm.start();	

        int j1 = -1;
        int j2 = -1;
        double min12 = std::numeric_limits<double>::max();

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


        // start loop over all objects
        while (Nstep > 0) { 

                min12 = std::numeric_limits<double>::max();

                /*
                cout << "Nstep=" << Nstep << endl;
                for (int i = 0; i < size; i++) { 
                        ParticleD *p1 = list[i];
                        cout << i << " is_consider[i]=" << is_consider[i] << " E=" << p1->e() << endl;

                }
                */



                // find smallest d12.
                // this is after reseting to a new jet
                if (!merged) { 

                for (int i = 0; i < size-1; i++) { 
                        if (is_consider[i]<=0) continue; 
                        for (int j = i + 1; j < size; j++) {  
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
                  for (int j = 0; j < size; j++) {
                                if (is_consider[j]<=0 || j==j1) continue;
                                if (ktdistance12[j1][j] < min12) {
                                        min12 = ktdistance12[j1][j];
                                        j1 = j1;
                                        j2 = j;
                                }
                        }
                } // end of min finding 

                  /*
                  if (Nstep>2) {
                    cout << "working on the last unmerged particle " << endl;
                    cout << "current pair=" << j1 << " " << j2 << endl;
                  }
                  */

                if (merged==false && Nstep==1) break;


                // find min distance to the beam
                double min1 = ktdistance1[j1];
                if (ktdistance1[j2]<min1) {min1 = ktdistance1[j2];};
                

                // make the decision about this particle 
                merged=false;
                if (min12<min1) merged=true;
         
                 if (merged) {   // merge particles  
                        //if (Nstep==1) cout << "Merge " << Nstep << " j1=" << j1 << " j2=" << j2<< endl;
                        ParticleD *p1 = list[j1];
                        ParticleD *p2 = list[j2];
                        if (j1 != j2)  p1->add(p2,j2); // p1=p1+p2. Also keeps an index j2 
                        Nstep--;
                        list[j1] = p1; // replace with p1=p1+p2
                        is_consider[j2]=0; // p2, but keep in the list 
                        is_consider[j1]=is_consider[j1]+1;
                        // recalculate distance for this particle
                        ktdistance1[j1] = getKtDistance1(p1);
                        for (int i = 0; i < size; i++) { 
                                if (is_consider[i]<=0) continue;  
                                ParticleD *pp1 = list[i];
                                ktdistance12[j1][i] = getKtDistance12(p1, pp1);
                        }

                }


              
               if (!merged) {   // add this to the jet 
                        is_consider[j1] = -1;
                        ParticleD *pj = list[j1];
                        Nstep--;
                        if (pj->getPt() > m_minpt)
                        {
                                jets.push_back(pj); // fill jets
                        }

                } 

        } // end loop 



               // attempt to deal with unmeargable particle
               int ins=-1;
               for (int i = 0; i < size; i++)
                        if (is_consider[i]==1) {ins=i;};
                     //cout << "working on the last unmerged particle " << endl;
                     ParticleD *p2 = list[ins];
                     if (m_debug) cout << "Unmergable particle=" << ins << " rap=" << p2->getRapidity() << " phi=" <<  p2->getPhi() << endl;
                     min12 = std::numeric_limits<double>::max();
                     for (unsigned int j = 0; j < jets.size(); j++) {
                                ParticleD *lp = jets[j];
                                double d=getDistance(p2, lp);
                                if (d<min12) { j1=j; min12 =d; }; 
                                }
                    if (m_debug) cout << "Unmergable distance R=" << min12 << endl;
                    if (min12<m_R) {
                           if (m_debug) cout << " --> Particle merged" << endl;
                           ParticleD *lp = jets[j1];
                           lp->add(p2,j1);
                           is_consider[ins] = 0;
                     }




        if (m_debug) { 
               tm.stop();
               std::cout << "  --> Final time for calculation (ms):" << tm.duration() << std::endl;
               tm.start();
               std::cout << "  --> Nr of jets : " << jets.size() << std::endl;
               printJets();
               // sanity test. All particles were merged
               int nn=0; int ins=-1;
               for (int i = 0; i < size; i++) 
                        if (is_consider[i]==1) {nn++; ins=i;};
               if (nn != 0)  { std::cout << "--> ERROR: =" << nn << " particles were not assigned.." << std::endl; 
                               std::cout << "--> ERROR: particle with ID=" << ins << " is lost!" << std::endl; 
                             }

        }

	return jets;

}



/** @brief Get sorted jets
    @return list with sorted jets 
    @author S.Chekanov 
    @version 1.0 02/02/14
*/
std::vector<ParticleD*> KT::getJetsSorted() { 
       sort(jets.begin(), jets.end(), comp);
       return jets;
}


/** Print sorted jets
*/
void KT::printJets() { 

        std::vector<ParticleD*> sjets = getJetsSorted();

        std::cout << "# Nr of jets=" << sjets.size() << std::endl;
        for (unsigned int i = 0; i < sjets.size(); i++)
        {
                ParticleD *lp = sjets[i];
                double phi = lp->phi();
                std::vector<int> con = lp->getConstituents();
                if (phi < 0) phi = PI + phi; 
                double srap=lp->getRapidity();
                double pt=lp->getPt();
                std::cout << "n=" << i << " y=" << srap << " phi=" << phi << " pt=" << pt << " const=" << con.size() << std::endl;

        }


}

double KT::phiAngle(double phi)
{
        if (phi > PI2)
        {
                phi -= (PI2);
        }
        if (phi < -PI2)
        {
                phi += (PI2);
        }
        return phi;
}


// calculate KT distance between any 2 particles
double KT::getKtDistance12(ParticleD *a, ParticleD *b)
{
	double rsq, esq, deltaEta, deltaPhi;
	deltaEta = b->getRapidity() - a->getRapidity();
	double phi1 = a->getPhi();
	double phi2 = b->getPhi();
	//deltaPhi = phiAngle(phi2 - phi1);
        deltaPhi = phi2 - phi1;
        if (deltaPhi>PI) deltaPhi=PI2-deltaPhi;
        if (deltaPhi<-PI) deltaPhi=PI2+deltaPhi;

        //deltaPhi = phi2 - phi1;
        //if(deltaPhi >= PI) deltaPhi = std::fmod(PI+deltaPhi, PI2) - PI;
        //else if(deltaPhi < -PI) deltaPhi = -std::fmod(PI-deltaPhi, PI2) + PI;


	rsq = (deltaEta*deltaEta + deltaPhi*deltaPhi);
	esq = 0;
	if (m_mode == 1)
	{
		esq = std::min(a->getPt2(), b->getPt2()); // kT
	}
	else if (m_mode == 0)
	{
		esq = std::min(a->getPt(), b->getPt()); // C-A
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


// get deltaR distance
double KT::getDistance(ParticleD *a, ParticleD *b)
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

// calculate distance to the beam
double KT::getKtDistance1(ParticleD *a)
{
	if (m_mode == 1)
	{
		return a->getPt2();
	}
	else if (m_mode == 0)
	{
		return a->getPt();
	}
	else if (m_mode == -1)
	{
		return (1.0 / a->getPt2());
	}
	return a->getPt2();
}


// set debugging mode
void KT::setDebug(bool debug)
{

	m_debug = debug;
}

