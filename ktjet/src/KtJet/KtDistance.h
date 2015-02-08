#ifndef KTJET_KTDISTANCE_H
#define KTJET_KTDISTANCE_H

#include "KtJet/KtDistanceInterface.h"
#include <string>
#include "KtJet/KtUtil.h"
#include "KtJet/KtLorentzVector.h"


namespace KtJet {
  /**
   *  Function object to calculate Kt for jets and pairs.

   @author J.Butterworth J.Couchman B.Cox B.Waugh
  */

  /** Get required KtDistance object given integer argument
   */
  KtDistance* getDistanceScheme(int dist, int collision_type);

  class KtDistanceAngle : public KtDistance {
  public:
    KtDistanceAngle(int collision_type=1);
    virtual ~KtDistanceAngle(){}
    /** Jet Kt */
    KtFloat operator()(const KtLorentzVector &) const;
    /** Pair Kt */
    KtFloat operator()(const KtLorentzVector &, const KtLorentzVector &) const;
    /** Name of scheme */
    std::string name() const;
  private:
    int m_type;
    std::string m_name;
  };


  class KtDistanceDeltaR : public KtDistance {
  public:
    KtDistanceDeltaR(int collision_type=1);
    virtual ~KtDistanceDeltaR(){};
    /** Jet Kt */
    KtFloat operator()(const KtLorentzVector &) const;
    /** Pair Kt */
    KtFloat operator()(const KtLorentzVector &, const KtLorentzVector &) const;
    /** Name of scheme */
    std::string name() const;
  private:
    int m_type;
    std::string m_name;
  };


  class KtDistanceQCD : public KtDistance {
  public:
    KtDistanceQCD(int collision_type=1);
    virtual ~KtDistanceQCD(){};
    /** Jet Kt */
    KtFloat operator()(const KtLorentzVector &) const;
    /** Pair Kt */
    KtFloat operator()(const KtLorentzVector &, const KtLorentzVector &) const;
    /** Name of scheme */
    std::string name() const;
  private:
    int m_type;
    std::string m_name;
  };

}//end of namespace

#endif
