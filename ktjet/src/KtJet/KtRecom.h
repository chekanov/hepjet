#ifndef KTJET_KTRECOM_H
#define KTJET_KTRECOM_H

#include <string>
#include "KtJet/KtUtil.h"
#include "KtJet/KtLorentzVector.h"
#include "KtJet/KtRecomInterface.h"


namespace KtJet {
  /**
   *  Function object to combine 4-momenta
   *  @author J.Butterworth J.Couchman B.Cox B.Waugh
  */

  /** Get required KtRecom object given integer argument
   */
  KtRecom* getRecomScheme(int recom);

  class KtRecomE : public KtRecom {
  public:
    KtRecomE();
    virtual ~KtRecomE(){};
    /** Return merged 4-momentum */
    CLHEP::HepLorentzVector operator()(const CLHEP::HepLorentzVector &, const CLHEP::HepLorentzVector &) const;
    /** Process input 4-momentum */
    KtLorentzVector operator()(const KtLorentzVector &) const;
    /** Name of scheme */
    std::string name() const;
  private:
    std::string m_name;
  };


  class KtRecomPt : public KtRecom {
  public:
    KtRecomPt();
    virtual ~KtRecomPt(){};
    /** Return merged 4-momentum */
    CLHEP::HepLorentzVector operator()(const CLHEP::HepLorentzVector &, const CLHEP::HepLorentzVector &) const;
    /** Process input 4-momentum */
    KtLorentzVector operator()(const KtLorentzVector &) const;
    /** Name of scheme */
    std::string name() const;
  private:
    std::string m_name;
  };


  class KtRecomPt2 : public KtRecom {
  public:
    KtRecomPt2();
    virtual ~KtRecomPt2(){};
    /** Return merged 4-momentum */
    CLHEP::HepLorentzVector operator()(const CLHEP::HepLorentzVector &, const CLHEP::HepLorentzVector &) const;
    /** Process input 4-momentum */
    KtLorentzVector operator()(const KtLorentzVector &) const;
    /** Name of scheme */
    std::string name() const;
  private:
    std::string m_name;
  };


  class KtRecomEt : public KtRecom {
  public:
    KtRecomEt();
    virtual ~KtRecomEt(){};
    /** Return merged 4-momentum */
    CLHEP::HepLorentzVector operator()(const CLHEP::HepLorentzVector &, const CLHEP::HepLorentzVector &) const;
    /** Process input 4-momentum */
    KtLorentzVector operator()(const KtLorentzVector &) const;
    /** Name of scheme */
    std::string name() const;
  private:
    std::string m_name;
  };


  class KtRecomEt2 : public KtRecom {
  public:
    KtRecomEt2();
    virtual ~KtRecomEt2(){};
    /** Return merged 4-momentum */
    CLHEP::HepLorentzVector operator()(const CLHEP::HepLorentzVector &, const CLHEP::HepLorentzVector &) const;
    /** Process input 4-momentum */
    KtLorentzVector operator()(const KtLorentzVector &) const;
    /** Name of scheme */
    std::string name() const;
  private:
    std::string m_name;
  };

}//end of namespace

#endif
