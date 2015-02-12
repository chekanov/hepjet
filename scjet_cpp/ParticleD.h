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
	Compute transverse energy.
	**/
	double et();

        /**	
	 Compute rapidity. 0.5*log( (m+z)/(m-z) );
        **/	
	virtual double rapidity();

	/**
	 Compute magnitude sqrt(px**2+py**2+pz**2)
        **/	
	double mag();

	/**
	 Compute pT**2  
	**/
	double perp2();


	/**
	Compute Transverse momentum (pT)   
        **/	
	double perp();



	void setEnergy(double energy);


	/**
	 Comparator. using perp2  for comparison (in increasing order)
        **/
	int compareTo(ParticleD *o);

	ParticleD *copy(ParticleD *o);

	/// <summary>
	/// Get  px. 
	/// </summary>
	/// <returns> px component </returns>
	double px();

	/// <summary>
	/// Get  py. 
	/// </summary>
	/// <returns> py component </returns>
	double py();

	/// <summary>
	/// Get  pz. 
	/// </summary>
	/// <returns> pz component </returns>
	double pz();


	/// <summary>
	/// Get  cached rapidity 
	/// </summary>
	/// <returns> rapidity  </returns>
	double getRapidity();


	/// <summary>
	/// Get  cached pt**2. 
	/// </summary>
	/// <returns> pt**2  </returns>
        double getPt2();

	/// <summary>
	/// Get  cached  pT. 
	/// </summary>
	/// <returns> pt transverse momentum  </returns>
	double getPt();


	/// <summary>
	/// Get  cached phi 
	/// </summary>
	/// <returns> cached phi  </returns>
	double getPhi();

	/// <summary>
	/// Get indexes of constituents, iff filled with the "add" method
	/// </summary>
	/// <returns> list of constituents.  </returns>
	std::vector<int> getConstituents();

	/// Compute Phi
	/// @return
	/// </summary>
	double phi();


	/// <summary>
	/// Get energy. 
	/// </summary>
	/// <returns> energy component </returns>
	double e();


        void setPx(double m_px);


	void setPy(double m_py);


	void setPz(double m_pz);


	void setRapidity(double m_rapidity);


	void setPhi(double m_phi);


	void setPt2(double m_pt2);


	/// <summary>
	/// The method precomputers Phi, Rapidity and Et2 and store them. 
	/// Such caching makes faster computations. Use getRapidity(), getPhi(), getEt2() methods
	/// to return such values. 
	/// </summary>
	void cachePhiRapidity();

	/// <summary>
	/// Add to this particle another particle. The method also precomputes rapidity, phi and et2
	/// for fast retrival using getPhi, getRapidity methods. </summary>
	/// <param name="a"> </param>
	void add(ParticleD *a);

	/// <summary>
	/// Add to this particle another particle. Also add index of the added particle, which will be stored as an array. 
	/// The method also precomputes rapidity, phi and et2
	/// for fast retrival using getPhi, getRapidity methods. </summary>
	/// <param name="a"> </param>
	/// <param name="index"> index of the particle to be stored </param>
	void add(ParticleD *a, int index);


};

#endif
