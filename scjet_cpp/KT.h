#ifndef _KT_H_
#define _KT_H_


#include "ParticleD.h"
#include <string>
#include <vector>
#include <iostream>
#include <limits>

/**
 Longitudinally-invariant kT, anti-KT and  Cambridge/Aachen clustering algorithms. 
 The algorithm uses rapidity-phi for the distance parameter and double values for merging. 
 This class uses double values for calculations, caching and requires more memory.
 This algorithm is similar to the FastJet http://fastjet.fr/ implementation that uses rapidity.  
 @author S.Chekanov
**/

class KT
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
	*            Currently only E-scheme is implemented. </param>
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
        *           Some small difference for very soft jets may exist compared to the traditonal anti-kT. 
        **/
	public: KT(double R=0.6, int recom=1, int mode=-1, double minpt=5.0, bool isfast=true);



	/** Run the jet algorithm using the list of particles. 
	* @param list
	*            list with particles
	* @returns final jets without sorting.
        **/ 
	virtual std::vector<ParticleD*> buildJets(std::vector<ParticleD*> &list);




	/// <summary>
	/// Get jets after  sorting in jet pT. Run  buildJets before calling this method. 
	/// </summary>
	/// <returns> list with sorted jets </returns>

	virtual std::vector<ParticleD*> getJetsSorted();

	/// <summary>
	/// Print the kT jets for debugging.
	/// </summary>
	virtual void printJets();

private:
	virtual double phiAngle(double phi);

	/// <summary>
	/// Calculate delta R distance.
	/// </summary>
	/// <param name="a">
	///            input particle </param>
	/// <param name="b">
	///            input particle </param>
	/// <param name="p">
	///            power parameter </param>
	/// <returns> Kt distance </returns>
public:
	virtual double getKtDistance12(ParticleD *a, ParticleD *b);

	/// <summary>
	/// This is the KT distance to the beam (assuming Z=Y=0). 
	/// The distance measure depends on the mode parameter. </summary>
	/// <param name="a">
	///            particle </param>
	/// <returns> kT distance </returns>
	virtual double getKtDistance1(ParticleD *a);


        /// <summary>
        /// This is the distance in Rapidity-Phi. 
        virtual double getDistance(ParticleD *a, ParticleD *b);


	/// <summary>
	/// Print debugging information. It shows how much time spend to make jets in ms. </summary>
	/// <param name="debug"> true if printing benchmark information.  </param>
	virtual void setDebug(bool debug);


};

#endif
