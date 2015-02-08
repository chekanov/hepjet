// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of methods of the HepRotation class which
// were introduced when ZOOM PhysicsVectors was merged in, and which involve
// Euler Angles representation.
//
// Apr 28, 2003  mf  Modified way of computing Euler angles to avoid flawed
//                   answers in the case where theta is near 0 of pi, and
//                   the matrix is not a perfect rotation (due to roundoff).

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/defs.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/EulerAngles.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>
#include <stdlib.h>

using std::abs;

namespace CLHEP  {

static inline double safe_acos (double x) {
  if (abs(x) <= 1.0) return acos(x);
  return ( (x>0) ? 0 : CLHEP::pi );
}

// ----------  Constructors and Assignment:

// Euler angles

HepRotation & HepRotation::set(double phi, double theta, double psi) {

  register double sinPhi   = sin( phi   ), cosPhi   = cos( phi   );
  register double sinTheta = sin( theta ), cosTheta = cos( theta );
  register double sinPsi   = sin( psi   ), cosPsi   = cos( psi   );

  rxx =   cosPsi * cosPhi - cosTheta * sinPhi * sinPsi;
  rxy =   cosPsi * sinPhi + cosTheta * cosPhi * sinPsi;
  rxz =   sinPsi * sinTheta;

  ryx = - sinPsi * cosPhi - cosTheta * sinPhi * cosPsi;
  ryy = - sinPsi * sinPhi + cosTheta * cosPhi * cosPsi;
  ryz =   cosPsi * sinTheta;

  rzx =   sinTheta * sinPhi;
  rzy = - sinTheta * cosPhi;
  rzz =   cosTheta;

  return  *this;

}  // Rotation::set(phi, theta, psi)

HepRotation::HepRotation( double phi, double theta, double psi ) 
{
  set (phi, theta, psi);
}
HepRotation & HepRotation::set( const HepEulerAngles & e ) {
  return set(e.phi(), e.theta(), e.psi());
}
HepRotation::HepRotation ( const HepEulerAngles & e ) 
{
  set(e.phi(), e.theta(), e.psi());
}


 
double HepRotation::phi  () const {

  double s2 =  1.0 - rzz*rzz;
  if (s2 < 0) {
    ZMthrowC ( ZMxpvImproperRotation (
        "HepRotation::phi() finds | rzz | > 1 "));
    s2 = 0;
  }
  const double sinTheta = sqrt( s2 );

  if (sinTheta < .01) { // For theta close to 0 or PI, use the more stable
  			// algorithm to get all three Euler angles
    HepEulerAngles ea = eulerAngles();
    return ea.phi();
  }
  
  const double cscTheta = 1/sinTheta;
  double cosabsphi =  - rzy * cscTheta;
  if ( fabs(cosabsphi) > 1 ) {	// NaN-proofing
    ZMthrowC ( ZMxpvImproperRotation (
      "HepRotation::phi() finds | cos phi | > 1 "));
    cosabsphi = 1;
  }
  const double absPhi = acos ( cosabsphi );
  if (rzx > 0) {
    return   absPhi;
  } else if (rzx < 0) {
    return  -absPhi;
  } else {
    return  (rzy < 0) ? 0 : CLHEP::pi;
  }

} // phi()

double HepRotation::theta() const {

  return  safe_acos( rzz );

} // theta()

double HepRotation::psi  () const {

  double sinTheta;
  if ( fabs(rzz) > 1 ) {	// NaN-proofing
    ZMthrowC ( ZMxpvImproperRotation (
      "HepRotation::psi() finds | rzz | > 1"));
    sinTheta = 0;
  } else { 
    sinTheta = sqrt( 1.0 - rzz*rzz );
  }
  
  if (sinTheta < .01) { // For theta close to 0 or PI, use the more stable
  			// algorithm to get all three Euler angles
    HepEulerAngles ea = eulerAngles();
    return ea.psi();
  }

  const double cscTheta = 1/sinTheta;
  double cosabspsi =  ryz * cscTheta;
  if ( fabs(cosabspsi) > 1 ) {	// NaN-proofing
    ZMthrowC ( ZMxpvImproperRotation (
      "HepRotation::psi() finds | cos psi | > 1"));
    cosabspsi = 1;
  }
  const double absPsi = acos ( cosabspsi );
  if (rxz > 0) {
    return   absPsi;
  } else if (rxz < 0) {
    return  -absPsi;
  } else {
    return  (ryz > 0) ? 0 : CLHEP::pi;
  }

} // psi()


// Helpers for eulerAngles():

static		     
void correctByPi ( double& psi, double& phi ) {
  if (psi > 0) {
    psi -= CLHEP::pi;
  } else {
    psi += CLHEP::pi;
  }
  if (phi > 0) {
    phi -= CLHEP::pi;
  } else {
    phi += CLHEP::pi;
  }  
}

static
void correctPsiPhi ( double rxz, double rzx, double ryz, double rzy, 
		     double& psi, double& phi ) {

  // set up quatities which would be positive if sin and cosine of
  // psi and phi were positive:
  double w[4];
  w[0] = rxz; w[1] = rzx; w[2] = ryz; w[3] = -rzy;

  // find biggest relevant term, which is the best one to use in correcting.
  double maxw = abs(w[0]); 
  int imax = 0;
  for (int i = 1; i < 4; ++i) {
    if (abs(w[i]) > maxw) {
      maxw = abs(w[i]);
      imax = i;
    }
  }
  // Determine if the correction needs to be applied:  The criteria are 
  // different depending on whether a sine or cosine was the determinor: 
  switch (imax) {
    case 0:
      if (w[0] > 0 && psi < 0)           correctByPi ( psi, phi );
      if (w[0] < 0 && psi > 0)           correctByPi ( psi, phi );
      break;
    case 1:
      if (w[1] > 0 && phi < 0)           correctByPi ( psi, phi );
      if (w[1] < 0 && phi > 0)           correctByPi ( psi, phi );
      break;
    case 2:
      if (w[2] > 0 && abs(psi) > CLHEP::halfpi) correctByPi ( psi, phi );    
      if (w[2] < 0 && abs(psi) < CLHEP::halfpi) correctByPi ( psi, phi );    
      break;
    case 3:
      if (w[3] > 0 && abs(phi) > CLHEP::halfpi) correctByPi ( psi, phi );    
      if (w[3] < 0 && abs(phi) < CLHEP::halfpi) correctByPi ( psi, phi );    
      break;
  }          
}


HepEulerAngles HepRotation::eulerAngles() const {

  // Please see the mathematical justification in eulerAngleComputations.ps

  double phi, theta, psi;
  double psiPlusPhi, psiMinusPhi;
  
  theta = safe_acos( rzz );
  
  if (rzz > 1 || rzz < -1) {
    ZMthrowC ( ZMxpvImproperRotation (
        "HepRotation::eulerAngles() finds | rzz | > 1 "));
  }
  
  double cosTheta = rzz;
  if (cosTheta > 1)  cosTheta = 1;
  if (cosTheta < -1) cosTheta = -1;

  if (cosTheta == 1) {
    psiPlusPhi = atan2 ( rxy - ryx, rxx + ryy );
    psiMinusPhi = 0;     

  } else if (cosTheta >= 0) {

    // In this realm, the atan2 expression for psi + phi is numerically stable
    psiPlusPhi = atan2 ( rxy - ryx, rxx + ryy );

    // psi - phi is potentially more subtle, but when unstable it is moot
    double s = -rxy - ryx; // sin (psi-phi) * (1 - cos theta)
    double c =  rxx - ryy; // cos (psi-phi) * (1 - cos theta)
    psiMinusPhi = atan2 ( s, c );
        
  } else if (cosTheta > -1) {

    // In this realm, the atan2 expression for psi - phi is numerically stable
    psiMinusPhi = atan2 ( -rxy - ryx, rxx - ryy );

   // psi + phi is potentially more subtle, but when unstable it is moot
    double s = rxy - ryx; // sin (psi+phi) * (1 + cos theta)
    double c = rxx + ryy; // cos (psi+phi) * (1 + cos theta)
    psiPlusPhi = atan2 ( s, c );

  } else { // cosTheta == -1

    psiMinusPhi = atan2 ( -rxy - ryx, rxx - ryy );
    psiPlusPhi = 0;

  }
  
  psi = .5 * (psiPlusPhi + psiMinusPhi); 
  phi = .5 * (psiPlusPhi - psiMinusPhi); 

  // Now correct by pi if we have managed to get a value of psiPlusPhi
  // or psiMinusPhi that was off by 2 pi:
  correctPsiPhi ( rxz, rzx, ryz, rzy, psi, phi );
  
  return  HepEulerAngles( phi, theta, psi );

} // eulerAngles()



void HepRotation::setPhi (double phi) {
  set ( phi, theta(), psi() );
}

void HepRotation::setTheta (double theta) {
  set ( phi(), theta, psi() );
}

void HepRotation::setPsi (double psi) {
  set ( phi(), theta(), psi );
}

}  // namespace CLHEP

