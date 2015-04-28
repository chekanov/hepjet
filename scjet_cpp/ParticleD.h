#ifndef _ParticleD_H_
#define _ParticleD_H_


#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>


using namespace std;

/**
 A class representing a jet or particle with pre-computed px,py,pz,e (double precision). 
 It uses double types to keep information on 4-meomenta. The merging is done in the 4-momentum space. 
 The class has a minimum dynamic computation to minimize CPU. 
 
 @author S.Chekanov 
 
**/



class ParticleD {

private:
	double m_px, m_py, m_pz;
	double energy;
	static const double PI2 = 6.28318530716;
	double m_phi, m_pt2, m_rapidity;
	std::vector<int> consts;
        bool operator<(ParticleD& score) const;

public:
	ParticleD();

        ~ParticleD();

	ParticleD(double m_px, double m_py, double m_pz, double energy);

	void setPxPyPzE(double m_px, double m_py, double m_pz, double energy);

        /**	
	* Compute pseudorapidity.
        **/	
	double eta();

        /**	
	* Compute transverse energy squared.
        **/	
	double et2();

	/**
	*Compute transverse energy.
	**/
	double et();

        /**
        *Compute mass.
        **/
        double m();

        /**
        *Compute mass.
        **/
        double mass();

        /**	
	* Compute rapidity. 0.5*log( (m+z)/(m-z) );
        **/	
	virtual double rapidity();

	/**
	* Compute magnitude sqrt(px**2+py**2+pz**2)
        **/	
	double mag();

	/**
	* Compute pT**2  
	**/
	double perp2();

	/**
	*Compute Transverse momentum (pT)   
        **/	
	double perp();

        /**
        * Set transverse momentum (pT)   
        **/
	void setEnergy(double energy);

	/**
	* Comparator. using perp2  for comparison (in increasing order)
        **/
	int compareTo(ParticleD *o);

	ParticleD *copy(ParticleD *o);

	// Get  px. 
	double px();

	// Get  py. 
	double py();

	// Get  pz. 
	double pz();

	// Get  cached rapidity 
	double getRapidity();
	
	// Get  cached pt**2. 
        double getPt2();

	// Get  cached  pT. 
	double getPt();

	// Get  cached phi 
	double getPhi();

	// Get indexes of constituents, iff filled with the "add" method
	std::vector<int> getConstituents();

	// Compute Phi
	double phi();

	// Get energy. 
	double e();

        void setPx(double m_px);

	void setPy(double m_py);

	void setPz(double m_pz);

	void setRapidity(double m_rapidity);

	void setPhi(double m_phi);

	void setPt2(double m_pt2);

	/** The method precomputers Phi, Rapidity and Et2 and store them. 
	* Such caching makes faster computations. Use getRapidity(), getPhi(), getEt2() methods
	* to return such values. 
        */
 	void cachePhiRapidity();

	/** Add to this particle another particle. The method also precomputes rapidity, phi and et2
	* for fast retrival using getPhi, getRapidity methods. </summary>
        **/
	void add(ParticleD *a);

	/** Add to this particle another particle. Also add index of the added particle, which will be stored as an array. 
	*   The method also precomputes rapidity, phi and et2
	*   for fast retrival using getPhi, getRapidity methods. </summary>
        */
	void add(ParticleD *a, int index);


};

#endif
