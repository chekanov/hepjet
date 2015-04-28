// -*- C++ -*-
// $Id: LorentzVectorB.cc,v 1.1 2006/12/19 19:16:09 madgraph Exp $
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the HepLorentzVector class:
// Those methods originating in ZOOM dealing with simple boosts and rotations.
// Use of one of these methods will not force loading of the HepRotation or
// HepLorentzRotation class.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/defs.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ZMxpv.h"

namespace CLHEP  {

//-*********
// rotationOf()
//-*********

// Each of these is a shell over a rotate method.

HepLorentzVector rotationXOf
	(const HepLorentzVector & vec, double phi){
  HepLorentzVector vv (vec);
  return vv.rotateX (phi);
}

HepLorentzVector rotationYOf
	(const HepLorentzVector & vec, double phi){
  HepLorentzVector vv (vec);
  return vv.rotateY (phi);
}

HepLorentzVector rotationZOf
	(const HepLorentzVector & vec, double phi){
  HepLorentzVector vv (vec);
  return vv.rotateZ (phi);
}

//-********
// boost
//-********

HepLorentzVector & HepLorentzVector::boost 
			( const Hep3Vector & axis,  double beta ) {
  if (beta==0) {
    return *this; // do nothing for a 0 boost
  }
  double r2 = axis.mag2();
  if ( r2 == 0 ) {
    ZMthrowA (ZMxpvZeroVector(
      "A zero vector used as axis defining a boost -- no boost done"));
    return *this;
  } 
  double b2 = beta*beta;
  if (b2 >= 1) {
    ZMthrowA (ZMxpvTachyonic(
      "LorentzVector boosted with beta >= 1 (speed of light) -- \n"
      "no boost done"));
  } else {
    Hep3Vector u = axis.unit();
    register double gamma = sqrt(1./(1.-b2));
    register double betaDotV = u.dot(pp)*beta;
    register double tt = ee;

    ee = gamma * (tt + betaDotV);
    pp += ( ((gamma-1)/b2)*betaDotV*beta + gamma*beta*tt ) * u;
    // Note:  I have verified the behavior of this even when beta is very
    //        small -- (gamma-1)/b2 becomes inaccurate by O(1), but it is then
    //        multiplied by O(beta**2) and added to an O(beta) term, so the
    //        inaccuracy does not affect the final result.
  }
  return *this;
} /* boost (axis, beta) */

}  // namespace CLHEP
