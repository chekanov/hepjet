// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is part of the implementation of the HepLorentzVector class:
// Those methods which originated from ZOOM and which deal with relativistic
// kinematic properties.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/defs.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ZMxpv.h"

#include <cmath>

namespace CLHEP  {

//-******************
// Metric flexibility
//-******************

ZMpvMetric_t HepLorentzVector::setMetric( ZMpvMetric_t m ) {
  ZMpvMetric_t oldMetric = (metric > 0) ? TimePositive : TimeNegative;
  if ( m == TimeNegative ) {
    metric = -1.0;
  } else {
    metric =  1.0;
  }
  return oldMetric;
}

ZMpvMetric_t HepLorentzVector::getMetric() {
  return ( (metric > 0) ? TimePositive : TimeNegative );
}

//-********
// plus
// minus
//-********

double HepLorentzVector::plus (const Hep3Vector & ref) const {
  double r = ref.mag();
  if (r == 0) {
    ZMthrowA (ZMxpvZeroVector(
      "A zero vector used as reference to LorentzVector plus-part"));
    return ee;
  }
  return ee + pp.dot(ref)/r;
} /* plus */

double HepLorentzVector::minus (const Hep3Vector & ref) const {
  double r = ref.mag();
  if (r == 0) {
    ZMthrowA (ZMxpvZeroVector(
      "A zero vector used as reference to LorentzVector minus-part"));
    return ee;
  }
  return ee - pp.dot(ref)/r;
} /* plus */

HepLorentzVector HepLorentzVector::rest4Vector() const {
  return HepLorentzVector (0, 0, 0, (t() < 0.0 ? -m() : m()));
}

//-********
// beta
// gamma
//-********

double HepLorentzVector::beta() const {
  if (ee == 0) {
    if (pp.mag2() == 0) {
      return 0;
    } else {
    ZMthrowA (ZMxpvInfiniteVector(
      "beta computed for HepLorentzVector with t=0 -- infinite result"));
    return 1./ee;
    }
  }
  if (restMass2() <= 0) {
    ZMthrowC (ZMxpvTachyonic(
      "beta computed for a non-timelike HepLorentzVector"));
        // result will make analytic sense but is physically meaningless
  }
  return sqrt (pp.mag2() / (ee*ee)) ;
} /* beta */

double HepLorentzVector::gamma() const {
  double v2 = pp.mag2();
  double t2 = ee*ee;
  if (ee == 0) {
    if (pp.mag2() == 0) {
      return 1;
    } else {
    ZMthrowC (ZMxpvInfiniteVector(
      "gamma computed for HepLorentzVector with t=0 -- zero result"));
    return 0;
    }
  }
  if (t2 < v2) {
    ZMthrowA (ZMxpvSpacelike(
      "gamma computed for a spacelike HepLorentzVector -- imaginary result"));
        // analytic result would be imaginary.
    return 0;
  } else if ( t2 == v2 ) {
    ZMthrowA (ZMxpvInfinity(
      "gamma computed for a lightlike HepLorentzVector -- infinite result"));
  }
  return 1./sqrt(1. - v2/t2 );
} /* gamma */


//-***************
// rapidity
// pseudorapidity
// eta
//-***************

double HepLorentzVector::rapidity() const {
  register double z = pp.getZ();
  if (fabs(ee) == fabs(z)) {
    ZMthrowA (ZMxpvInfinity(
      "rapidity for 4-vector with |E| = |Pz| -- infinite result"));
  }
  if (fabs(ee) < fabs(z)) {
    ZMthrowA (ZMxpvSpacelike(
      "rapidity for spacelike 4-vector with |E| < |Pz| -- undefined"));
    return 0;
  }
  double q = (ee + z) / (ee - z);
        //-| This cannot be negative now, since both numerator
        //-| and denominator have the same sign as ee.
  return .5 * log(q);
} /* rapidity */

double HepLorentzVector::rapidity(const Hep3Vector & ref) const {
  register double r = ref.mag2();
  if (r == 0) {
    ZMthrowA (ZMxpvZeroVector(
      "A zero vector used as reference to LorentzVector rapidity"));
    return 0;
  }
  register double vdotu = pp.dot(ref)/sqrt(r);
  if (fabs(ee) == fabs(vdotu)) {
    ZMthrowA (ZMxpvInfinity(
      "rapidity for 4-vector with |E| = |Pu| -- infinite result"));
  }
  if (fabs(ee) < fabs(vdotu)) {
    ZMthrowA (ZMxpvSpacelike(
      "rapidity for spacelike 4-vector with |E| < |P*ref| -- undefined "));
    return 0;
  }
  double q = (ee + vdotu) / (ee - vdotu);
  return .5 * log(q);
} /* rapidity(ref) */

double HepLorentzVector::coLinearRapidity() const {
  register double v = pp.mag();
  if (fabs(ee) == fabs(v)) {
    ZMthrowA (ZMxpvInfinity(
      "co-Linear rapidity for 4-vector with |E| = |P| -- infinite result"));
  }
  if (fabs(ee) < fabs(v)) {
    ZMthrowA (ZMxpvSpacelike(
      "co-linear rapidity for spacelike 4-vector -- undefined"));
    return 0;
  }
  double q = (ee + v) / (ee - v);
  return .5 * log(q);
} /* rapidity */

//-*************
// invariantMass
//-*************

double HepLorentzVector::invariantMass(const HepLorentzVector & w) const {
  double m2 = invariantMass2(w);
  if (m2 < 0) {
    // We should find out why:
    if ( ee * w.ee < 0 ) {
      ZMthrowA (ZMxpvNegativeMass(
        "invariant mass meaningless: \n"
        "a negative-mass input led to spacelike 4-vector sum" ));
      return 0;
    } else if ( (isSpacelike() && !isLightlike()) ||
                (w.isSpacelike() && !w.isLightlike()) ) {
      ZMthrowA (ZMxpvSpacelike(
        "invariant mass meaningless because of spacelike input"));
      return 0;
    } else {
      // Invariant mass squared for a pair of timelike or lightlike vectors
      // mathematically cannot be negative.  If the vectors are within the
      // tolerance of being lightlike or timelike, we can assume that prior
      // or current roundoffs have caused the negative result, and return 0
      // without comment.
      return 0;
    }
  }
  return (ee+w.ee >=0 ) ? sqrt(m2) : - sqrt(m2);
} /* invariantMass */

//-***************
// findBoostToCM
//-***************

Hep3Vector HepLorentzVector::findBoostToCM() const {
  return -boostVector();
} /* boostToCM() */

Hep3Vector HepLorentzVector::findBoostToCM (const HepLorentzVector & w) const {
  double t = ee + w.ee;
  Hep3Vector v = pp + w.pp;
  if (t == 0) {
    if (v.mag2() == 0) {
      return Hep3Vector(0,0,0);
    } else {
    ZMthrowA (ZMxpvInfiniteVector(
    "boostToCM computed for two 4-vectors with combined t=0 -- "
        "infinite result"));
    return Hep3Vector(v*(1./t)); // Yup, 1/0 -- that is how we return infinity
    }
  }
  if (t*t - v.mag2() <= 0) {
    ZMthrowC (ZMxpvTachyonic(
    "boostToCM  computed for pair of HepLorentzVectors with non-timelike sum"));
        // result will make analytic sense but is physically meaningless
  }
  return Hep3Vector(v * (-1./t));
} /* boostToCM(w) */

}  // namespace CLHEP

