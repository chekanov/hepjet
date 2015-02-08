#ifndef KTJET_KTUTIL_H
#define KTJET_KTUTIL_H

// Includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"


namespace KtJet{
  
#ifdef KTDOUBLEPRECISION
  typedef double KtFloat;
#else
  typedef float KtFloat;
#endif
  
  class KtLorentzVector;
  
  /** Phi angle forced into range -pi to +pi */
  KtFloat phiAngle(KtFloat testphi);
  
} //end of namespace

#endif
