#ifndef _KT_H_
#define _KT_H_


#include "ParticleD.h"
#include <string>
#include <vector>
#include <iostream>
#include <limits>

/**
* SCJet is an implementation of the longitudinally-invariant kT, anti-KT and  Cambridge/Aachen clustering algorithms. 
* The algorithm uses rapidity-phi for the distance parameter and double values for merging. 
* This class uses double values for calculations, caching and requires more memory.
* This algorithm is similar to the FastJet http://fastjet.fr/ implementation that uses rapidity.  
* The algorithm with disabled fast calculations  has identical outputs as FastJet, but slower by a factor 15.
* When fast calculations are enabled for the anti-kT mode, the algorithm has small difference for very soft proto-jets,
* but in this mode it is a factor 4 faster compared to the same algorithm with disabled fast calculations.   
* @author S.Chekanov
**/

class SCJet
{

private:
	int m_recom;
	double m_R;
	double m_R2;
	std::vector<ParticleD*> jets;
	bool m_debug;
	double m_minpt;
	int m_mode;
	bool m_fast;


	/**
	* Initialize calculations of the longitudinally invariant kT algorithm in inclusive mode. 
	* Jet can be clustered using Cambridge/Aachen or anti-kT approaches, depending on the "mode" parameter. 
	* The distance parameters are rapidity and phi. 
	* @param R
	*          distance measure defined in rapidity-phi. 
	* @param recom
	*            recombination scheme.
	*            1: The E-scheme Simple 4-vector addition.
	*            2: The pT-scheme. 
	*            3: The pT^2 scheme.
	*            Currently only E-scheme is implemented.
	* @param mode 
	*          clustering mode dij=min(kT_i^{2* mode},kT_j^{2* mode})).
	*          mode=1 means inclusive KT jet algorithm  
	*          mode=0 means Cambridge/Aachen jet algorithm 
	*          mode=-1 means anti-KT jet algorithm 
	* @param minpt
	*            min pT for final jets.
	*
	* @param isfast
	*           when true, enable fast calculations for anti-kT jets. They scale as N^2, unlike the traditional anti-kT N^3.
	*           Some small difference for very soft jets may exist compared to the traditional anti-kT. 
	**/
public: SCJet(double R=0.6, int recom=1, int mode=-1, double minpt=5.0, bool isfast=true);

	/** Run the jet algorithm using the list of particles.
	* @param list
	*            list with particles
	* @returns final jets without sorting.
	       **/ 
	virtual std::vector<ParticleD*> buildJets(std::vector<ParticleD*> &list);

	/**
	*Get jets after  sorting in jet pT. Run  buildJets before calling this method. 
	*@return list with sorted jets
	**/	
	virtual std::vector<ParticleD*> getJetsSorted();

	/**
	*Print the kT jets for debugging.
	**/
	virtual void printJets();

	/**
	 *This is the kT-distance to  another particle.
	 * @param a
	 *           particle
	 * @param b 
	 *           particle
	 * @return kT distance 
	 */
	virtual double getKtDistance12(ParticleD *a, ParticleD *b);

	/**
	*This is the kT-distance to the beam (assuming Z=Y=0). 
	*The distance measure depends on the mode parameter.
	* @param a
	*           particle
	* @return kT distance 
	*/
	virtual double getKtDistance1(ParticleD *a);


	/** Print debugging information.
	**/
	virtual void setDebug(bool debug);


        /** Get distance in y-phi between 2 particles. 
        **/
        virtual double getDistance(ParticleD *a, ParticleD *b); 

};

#endif
