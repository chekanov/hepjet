#include "KtJet/KtUtil.h"
#include "KtJet/KtLorentzVector.h"
#include <cmath>

namespace KtJet {
  /** Put phi in range [-pi,+pi]. No such function in CLHEP 1.7. (But is in 1.8.)
   */
KtFloat phiAngle(KtFloat testphi) {
  KtFloat phi = testphi;
  while (phi>M_PI) phi -= (2*M_PI);
  while (phi<-M_PI) phi += (2*M_PI);
  return phi;
}

}//end of namespace
