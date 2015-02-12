// KT is class to build KT-type of jets in inclusive mode
// S.Chekanov (ANL)

#include "KT.h"
#include "Timer.h"
#include <vector>
#include <algorithm>

using namespace std;
const double PI2 = 6.28318530716;
const double PI  = 3.14159265358;

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

	double ktdistance1[size];
	int    is_consider[size]; // 0-ignore, 1-original,  n>1 - merged, -1 - jet
	for (int m = 0; m < size; m++) { 
		is_consider[m] = 1;
		ParticleD *p1 = static_cast<ParticleD*>(list[m]);
		ktdistance1[m] = getKtDistance1(p1);
	}




        double ktdistance12[size][size];
        for (int i = 0; i < size - 1; i++)
        {
                ParticleD *p1 = static_cast<ParticleD*>(list[i]);
                for (int j = i + 1; j < size; j++)
                {
                        ParticleD *p2 = static_cast<ParticleD*>(list[j]);
                        ktdistance12[i][j] = getKtDistance12(p1, p2);
                        //System.out.println(ktdistance12[i][j]);
                }
        }

        if (m_debug) { 
               tm.stop();
               std::cout << "After creating a cache of distances (ms):" << tm.duration() << std::endl;
               tm.start();
        }


         int Nstep = size;
        // start loop over all objects
        while (Nstep > 0) { 
                int j1 = 0;
                int j2 = 1;
                double min12 = 1e+64;
                //
                // find smallest d12.
                //
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

                // find min distance to the beam
                int km = 0;
                double min1 = 1e+64;
/*
                for (int i = 0; i < size; i++) { 
                        if (is_consider[i]<=0) continue; 
                        if (ktdistance1[i] < min1)
                        {
                                min1 = ktdistance1[i];
                                km = i;
                        }
                }
*/

                km = j1;
                min1 = ktdistance1[j1];
                if (ktdistance1[j2]<min1) {min1 = ktdistance1[j2];
                                          km = j2; } 

                // make the decision about this particle 
                bool merged=false;
                if (min12<min1) merged=true;

         
                 if (merged) {   // merge particles  
                        ParticleD *p1 = static_cast<ParticleD*>(list[j1]);
                        ParticleD *p2 = static_cast<ParticleD*>(list[j2]);
                        if (j1 != j2)  p1->add(p2,j2); // p1=p1+p2. Also keeps an index j2 
                        Nstep--;
                        list[j1] = p1; // replace with p1=p1+p2
                        is_consider[j2] = 0; // p2, but keep in the list 
                        is_consider[j1]=is_consider[j1]+1;
                        // recalculate distance for this particle
                        ktdistance1[j1] = getKtDistance1(p1);
                        for (int i = 0; i < size; i++) { 
                                if (is_consider[i]<=0) continue;  
                                ParticleD *pp1 = static_cast<ParticleD*>(list[i]);
                                ktdistance12[j1][i] = getKtDistance12(p1, pp1);
                        }

                }


               if (!merged) {   // add this to the jet 
                        is_consider[km] = -1; // make as a jet
                        Nstep--;
                        ParticleD *pj = static_cast<ParticleD*>(list[km]);
                        if (pj->getPt() > m_minpt)
                        {
                                jets.push_back(pj); // fill jets
                                // System.out.println(pj.getPt() + " " + minpt);
                        }

                } 



                // end loop
        }


        if (m_debug) { 
               tm.stop();
               std::cout << "  --> Final time for calculation (ms):" << tm.duration() << std::endl;
               tm.start();
               std::cout << "  --> Nr of jets : " << jets.size() << std::endl;
               printJets();
               // sanity test. All particles were merged
               int nn=0;
               for (int i = 0; i < size; i++) 
                        if (!is_consider[i]) nn++; 
                if (nn != (int)list.size())  std::cout << "!!!! Error! Not all particles were assigned" << list.size() << " done=" << nn << std::endl; 


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
                if (phi < 0) phi = PI2 + phi; 
                double srap=lp->getRapidity();
                double et=lp->getPt();
                std::cout << "n=" << i << " y=" << srap << " phi=" << phi << " pt=" << et << " const=" << con.size() << std::endl;

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
        //if (deltaPhi>PI) deltaPhi=PI2-deltaPhi;
        //if (deltaPhi<-PI) deltaPhi=PI2+deltaPhi;
        deltaPhi = phi2 - phi1;
        if(deltaPhi >= PI) deltaPhi = std::fmod(PI+deltaPhi, PI2) - PI;
        else if(deltaPhi < -PI) deltaPhi = -std::fmod(PI-deltaPhi, PI2) + PI;


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

