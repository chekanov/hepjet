#ifndef _KT_H_
#define _KT_H_


#include "ParticleD.h"
#include <string>
#include <vector>
#include <iostream>
#include <limits>

/// <summary>
/// Longitudinally-invariant kT, anti-KT and  Cambridge/Aachen clustering algorithms. 
/// The algorithm uses rapidity-phi for the distance parameter and double values for merging. The input and output for this algorithm  is  <seealso cref="hephysics.jet.ParticleD"/> class. <para> 
/// </para>
/// This class uses double values for calculations, caching and requires more memory, compared to the light-weight <seealso cref="hephysics.jet.KTjet"/> class that uses floats and pseudo-rapidity to define the distance parameter. This implementation can access jet constituents.  <para></para> 
/// This algorithm is similar to the FastJet http://fastjet.fr/ implementation that uses rapidity.  Use light-weight <seealso cref="hephysics.jet.KTjet"/> class when using pseudo-rapidity and phi to define distance parameters. The method uses E-scheme to combine particles (p1+p2). More details is in http://arxiv.org/pdf/hep-ph/0210022v1.pdf. 
/// 
/// @author S.Chekanov
/// 
/// </summary>
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


	/// <summary>
	/// Initialize calculations of the longitudinally invariant kT algorithm in inclusive mode. 
	/// Jet can be clustered using Cambridge/Aachen or anti-kT approaches, depending on the "mode" parameter. 
	/// The distance parameters are rapidity and phi. 
	/// </summary>
	/// <param name="R">
	///            distance measure </param>
	/// <param name="recom">
	///            recombination scheme.<br>
	///            1: The E-scheme Simple 4-vector addition. <br>
	///            2: The pT-scheme. <br>
	///            3: The pT^2 scheme. <br>
	///            Currently only E-scheme is implemented. </param>
	/// <param name="mode"> 
	///          clustering mode dij=min(kT_i^{2* mode},kT_j^{2* mode})). <br>
	///          mode=1 means inclusive KT jet algorithm <br> 
	///          mode=0 means Cambridge/Aachen jet algorithm <br> 
	///          mode=-1 means anti-KT jet algorithm <br> </param>
	/// <param name="minpt">
	///            min pT for final jets. </param>
public:
	KT(double R, int recom, int mode, double minpt);


	/// <summary>
	/// Initialize calculations of the KT algorithm. Meaningful values are R=0.2- 1.
	/// Jets are clustered in rapidity and phi space. 
	/// </summary>
	/// <param name="R">
	///            distance measure </param>
	/// <param name="recom">
	///            recombination scheme.<br>
	///            1: The E-scheme Simple 4-vector addition. <br>
	///            2: The pT-scheme. <br>
	///            3: The pT^2 scheme. <br>
	///            Currently only E-scheme is implemented. </param>
	/// <param name="minpt">
	///            min pT for final jets. </param>
	KT(double R, int recom, double minpt);


	/// <summary>
	/// Initialize calculations of the KT algorithm. Meaningful values are R=0.2- 1.
	/// Jets are clustered in rapidity and phi space. The The E-scheme with 4-vector addition is used. 
	/// </summary>
	/// <param name="R">
	///            distance measure </param>
	/// <param name="minpt">
	///            min pT for final jets. </param>
	KT(double R, double minpt);




	/// <summary>
	/// Run the jet algorithm using the list of particles. 
	/// </summary>
	/// <param name="list">
	///            list with particles </param>
	/// <returns> final jets without sorting. </returns>
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
	/// Print debugging information. It shows how much time spend to make jets in ms. </summary>
	/// <param name="debug"> true if printing benchmark information.  </param>
	virtual void setDebug(bool debug);


};

#endif
