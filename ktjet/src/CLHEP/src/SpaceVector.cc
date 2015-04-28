// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// SpaceVector
//
// This is the implementation of those methods of the Hep3Vector class which
// originated from the ZOOM SpaceVector class.  Several groups of these methods
// have been separated off into the following code units:
//
// SpaceVectorR.cc	All methods involving rotation
// SpaceVectorD.cc	All methods involving angle decomposition
// SpaceVectorP.cc	Intrinsic properties and methods involving second vector
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/defs.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/ZMxpv.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>

namespace CLHEP  {

//-*****************************
//           - 1 -
// set (multiple components)
// in various coordinate systems
//
//-*****************************

void Hep3Vector::setSpherical (
		double r,
                double theta,
                double phi) {
  if ( r < 0 ) {
    ZMthrowC (ZMxpvNegativeR(
      "Spherical coordinates set with negative   R"));
    // No special return needed if warning is ignored.
  }
  if ( (theta < 0) || (theta > CLHEP::pi) ) {
    ZMthrowC (ZMxpvUnusualTheta(
      "Spherical coordinates set with theta not in [0, PI]"));
	// No special return needed if warning is ignored.
  }
  dz = r * cos(theta);
  double rho ( r*sin(theta));
  dy = rho * sin (phi);
  dx = rho * cos (phi);
  return;
} /* setSpherical (r, theta, phi) */

void Hep3Vector::setCylindrical (
 		double rho,
                double phi,
                double z) {
  if ( rho < 0 ) {
    ZMthrowC (ZMxpvNegativeR(
      "Cylindrical coordinates supplied with negative Rho"));
    // No special return needed if warning is ignored.
  }
  dz = z;
  dy = rho * sin (phi);
  dx = rho * cos (phi);
  return;
} /* setCylindrical (r, phi, z) */

void Hep3Vector::setRhoPhiTheta (
 		double rho,
                double phi,
                double theta) {
  if (rho == 0) {
    ZMthrowC (ZMxpvZeroVector(
      "Attempt set vector components rho, phi, theta with zero rho -- "
      "zero vector is returned, ignoring theta and phi"));
    dx = 0; dy = 0; dz = 0;
    return;
  }
  if ( (theta == 0) || (theta == CLHEP::pi) ) {
    ZMthrowA (ZMxpvInfiniteVector(
      "Attempt set cylindrical vector vector with finite rho and "
      "theta along the Z axis:  infinite Z would be computed"));
  }
  if ( (theta < 0) || (theta > CLHEP::pi) ) {
    ZMthrowC (ZMxpvUnusualTheta(
      "Rho, phi, theta set with theta not in [0, PI]"));
	// No special return needed if warning is ignored.
  }
  dz = rho / tan (theta);
  dy = rho * sin (phi);
  dx = rho * cos (phi);
  return;
} /* setCyl (rho, phi, theta) */

void Hep3Vector::setRhoPhiEta (
 		double rho,
                double phi,
                double eta ) {
  if (rho == 0) {
    ZMthrowC (ZMxpvZeroVector(
      "Attempt set vector components rho, phi, eta with zero rho -- "
      "zero vector is returned, ignoring eta and phi"));
    dx = 0; dy = 0; dz = 0;
    return;
  }
  double theta (2 * atan ( exp (-eta) ));
  dz = rho / tan (theta);
  dy = rho * sin (phi);
  dx = rho * cos (phi);
  return;
} /* setCyl (rho, phi, eta) */


//************
//    - 3 - 
// Comparisons
//
//************

int Hep3Vector::compare (const Hep3Vector & v) const {
  if       ( dz > v.dz ) {
    return 1;
  } else if ( dz < v.dz ) {
    return -1;
  } else if ( dy > v.dy ) {
    return 1;
  } else if ( dy < v.dy ) {
    return -1;
  } else if ( dx > v.dx ) {
    return 1;
  } else if ( dx < v.dx ) {
    return -1;
  } else {
    return 0;
  }
} /* Compare */


bool Hep3Vector::operator > (const Hep3Vector & v) const {
	return (compare(v)  > 0);
}
bool Hep3Vector::operator < (const Hep3Vector & v) const {
	return (compare(v)  < 0);
}
bool Hep3Vector::operator>= (const Hep3Vector & v) const {
	return (compare(v) >= 0);
}
bool Hep3Vector::operator<= (const Hep3Vector & v) const {
	return (compare(v) <= 0);
}



//-********
// Nearness
//-********

// These methods all assume you can safely take mag2() of each vector.
// Absolutely safe but slower and much uglier alternatives were 
// provided as build-time options in ZOOM SpaceVectors.
// Also, much smaller codes were provided tht assume you can square
// mag2() of each vector; but those return bad answers without warning
// when components exceed 10**90. 
//
// IsNear, HowNear, and DeltaR are found in ThreeVector.cc

double Hep3Vector::howParallel (const Hep3Vector & v) const {
  // | V1 x V2 | / | V1 dot V2 |
  double v1v2 = fabs(dot(v));
  if ( v1v2 == 0 ) {
    // Zero is parallel to no other vector except for zero.
    return ( (mag2() == 0) && (v.mag2() == 0) ) ? 0 : 1;
  }
  Hep3Vector v1Xv2 ( cross(v) );
  double abscross = v1Xv2.mag();
  if ( abscross >= v1v2 ) {
    return 1;
  } else {
    return abscross/v1v2;
  }
} /* howParallel() */

bool Hep3Vector::isParallel (const Hep3Vector & v,
                              double epsilon) const {
  // | V1 x V2 | **2  <= epsilon **2 | V1 dot V2 | **2
  // V1 is *this, V2 is v

  static const double TOOBIG = pow(2.0,507);
  static const double SCALE  = pow(2.0,-507);
  double v1v2 = fabs(dot(v));
  if ( v1v2 == 0 ) {
    return ( (mag2() == 0) && (v.mag2() == 0) );
  }
  if ( v1v2 >= TOOBIG ) {
    Hep3Vector sv1 ( *this * SCALE );
    Hep3Vector sv2 ( v * SCALE );
    Hep3Vector sv1Xsv2 = sv1.cross(sv2);
    double x2 = sv1Xsv2.mag2();
    double limit = v1v2*SCALE*SCALE;
    limit = epsilon*epsilon*limit*limit;
    return ( x2 <= limit );
  }

  // At this point we know v1v2 can be squared.

  Hep3Vector v1Xv2 ( cross(v) );
  if (  (fabs (v1Xv2.dx) > TOOBIG) ||
        (fabs (v1Xv2.dy) > TOOBIG) ||
        (fabs (v1Xv2.dz) > TOOBIG) ) {
    return false;
  }
                    
  return ( (v1Xv2.mag2()) <= ((epsilon * v1v2) * (epsilon * v1v2)) );

} /* isParallel() */


double Hep3Vector::howOrthogonal (const Hep3Vector & v) const {
  // | V1 dot V2 | / | V1 x V2 | 

  double v1v2 = fabs(dot(v));
	//-| Safe because both v1 and v2 can be squared
  if ( v1v2 == 0 ) {
    return 0;	// Even if one or both are 0, they are considered orthogonal
  }
  Hep3Vector v1Xv2 ( cross(v) );
  double abscross = v1Xv2.mag();
  if ( v1v2 >= abscross ) {
    return 1;
  } else {
    return v1v2/abscross;
  }

} /* howOrthogonal() */

bool Hep3Vector::isOrthogonal (const Hep3Vector & v,
			     double epsilon) const {
// | V1 x V2 | **2  <= epsilon **2 | V1 dot V2 | **2
// V1 is *this, V2 is v

  static const double TOOBIG = pow(2.0,507);
  static const double SCALE = pow(2.0,-507);
  double v1v2 = fabs(dot(v));
        //-| Safe because both v1 and v2 can be squared
  if ( v1v2 >= TOOBIG ) {
    Hep3Vector sv1 ( *this * SCALE );
    Hep3Vector sv2 ( v * SCALE );
    Hep3Vector sv1Xsv2 = sv1.cross(sv2);
    double x2 = sv1Xsv2.mag2();
    double limit = epsilon*epsilon*x2;
    double y2 = v1v2*SCALE*SCALE;
    return ( y2*y2 <= limit );
  }

  // At this point we know v1v2 can be squared.

  Hep3Vector eps_v1Xv2 ( cross(epsilon*v) );
  if (  (fabs (eps_v1Xv2.dx) > TOOBIG) ||
        (fabs (eps_v1Xv2.dy) > TOOBIG) ||
        (fabs (eps_v1Xv2.dz) > TOOBIG) ) {
    return true;
  }

  // At this point we know all the math we need can be done.

  return ( v1v2*v1v2 <= eps_v1Xv2.mag2() );

} /* isOrthogonal() */

double Hep3Vector::setTolerance (double tol) {
// Set the tolerance for Hep3Vectors to be considered near one another
  double oldTolerance (tolerance);
  tolerance = tol;
  return oldTolerance;
}


//-***********************
// Helper Methods:
//	negativeInfinity()
//-***********************

double Hep3Vector::negativeInfinity() const {
  // A byte-order-independent way to return -Infinity
  struct Dib {
    union {
      double d;
      unsigned char i[8];
    } u;
  };
  Dib negOne;
  Dib posTwo;
  negOne.u.d = -1.0;
  posTwo.u.d =  2.0;
  Dib value;
  int k;
  for (k=0; k<8; k++) {
    value.u.i[k] = negOne.u.i[k] | posTwo.u.i[k];
  }
  return value.u.d;
}

}  // namespace CLHEP


