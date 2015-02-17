#ifndef _ParticleD_H_
#define _ParticleD_H_


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>
#include <list>
#include <set>
#include <cstdlib>
#include <cstdio>

using namespace std;


const double PI2 = 2*M_PI;
const double PI  = M_PI;


template<class ForwardIt> ForwardIt next(ForwardIt it) { return ++it; }


/**
 A class representing a jet or particle with pre-computed px,py,pz,e (double precision). 
 It uses double types to keep information on 4-meomenta. The merging is done in the 4-momentum space. 
 The class has a minimum dynamic computation to minimize CPU. 
 @authors Ivan Pogrebnyak and Sergei Chekanov 
**/

class ParticleD {

private:
	double m_px, m_py, m_pz;
	double energy, m_d;
	double m_phi, m_pt2, m_rapidity;
	bool operator<(ParticleD& score) const;
	int  m_mode;
	int m_id;
	static int num;
	std::vector<ParticleD> consts;

public:

	ParticleD();

	~ParticleD(){};

	// also include distance mode
	ParticleD(double m_px, double m_py, double m_pz, double energy, int mode=1);

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

	ParticleD *copy(ParticleD *o);

	/// <summary>
	/// Get  px.
	/// </summary>
	/// <returns> px component </returns>
	double px() const { return m_px; };

	/// <summary>
	/// Get  py.
	/// </summary>
	/// <returns> py component </returns>
	double py() const { return m_py; };

	/// <summary>
	/// Get  pz.
	/// </summary>
	/// <returns> pz component </returns>
	double pz() const { return m_pz; };

	/// <summary>
	/// Get  cached rapidity
	/// </summary>
	/// <returns> rapidity  </returns>
	double getRapidity() const { return m_rapidity; };

	// return distance to the beam
	double getD() const { return m_d; }

	double getID() const { return m_id; }

	/// <summary>
	/// Get  cached pt**2.
	/// </summary>
	/// <returns> pt**2  </returns>
	double getPt2() const { return m_pt2; };

	/// <summary>
	/// Get  cached  pT.
	/// </summary>
	/// <returns> pt transverse momentum  </returns>
	double getPt()  const { return sqrt(m_pt2); };


	/// <summary>
	/// Get  cached phi
	/// </summary>
	/// <returns> cached phi  </returns>
	double getPhi()   const { return sqrt(m_phi); };


	/// <summary>
	/// Get indexes of constituents, iff filled with the "add" method
	/// </summary>
	/// <returns> list of constituents.  </returns>
	std::vector<ParticleD> getConstituents() const {return consts;};


        void  setConstituents(const std::vector<ParticleD> p) {
                 consts.insert(consts.end(), p.begin(), p.end());

        } 

       void  setConstituent(const ParticleD p) {
                 consts.push_back(p);

        }

	/// Compute Phi
	/// @return
	/// </summary>
	double phi();


	/// <summary>
	/// Get energy.
	/// </summary>
	/// <returns> energy component </returns>
	double e() const { return energy; };


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


	double distanceB();

	// calculated distance to another particle;
	double distance(const ParticleD &p) const {
		double deltaPhi=m_phi-p.m_phi;
		if (deltaPhi>PI)  deltaPhi=PI2-deltaPhi;
		if (deltaPhi<-PI) deltaPhi=PI2+deltaPhi;
		double deltaY=m_rapidity - p.m_rapidity;
		double sd;

		if (m_mode==-1) {
			sd=std::min(1.0/m_pt2,1.0/p.m_pt2); // anti-KT
		} else if (m_mode==0) {
			sd=1.0; // C/A
			// kT
		}else if (m_mode==1) {
			sd=std::min(m_pt2,p.m_pt2); // KT
		}

		return sd*(deltaY*deltaY+deltaPhi*deltaPhi);
	}


	// add 2 particles
	ParticleD operator+(const ParticleD &p) const {
        	ParticleD pp(m_px+p.m_px,m_py+p.m_py,m_pz+p.m_pz,energy+p.energy,m_mode);
		return pp;
	}


	bool operator<(const ParticleD &p) const {
		if (m_d==p.m_d) return ( m_id < p.m_id );
		else return ( m_d < p.m_d );
	}


	ostream& prt(ostream &out, int depth=0) const {
		for (int i=0;i<depth;++i) out << "  |";
		out << fixed
		<< setw(3) << m_id << ": ("
		<< setw(9) << m_px << ", "
		<< setw(9) << m_py << ", "
		<< setw(9) << m_pz << ", "
		<< setw(9) << energy << " | "
		<< setw(9) << m_pt2 << ", "
		<< setw(9) << m_d
		<< " )" << endl;
		return out;
	}


	/*
	ostream& operator<<(ostream &out, const ParticleD& p) {
	  const ios_base::fmtflags flags = out.flags();
	  const streamsize prec = out.precision(5);
	  p.prt(out);
	  out.flags(flags);
	  out.precision(prec);
	  return out;
}


	istream& operator>>(istream &in, set<ParticleD>& p) {
	  static double px, py, pz, E;
	  in >> px >> py >> pz >> E;
	  if (in) cout << *p.insert(ParticleD(px,py,pz,E)).first;
	  return in;
}
	*/



};

#endif
