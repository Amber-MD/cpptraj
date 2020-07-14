#ifndef INC_MDOPTS_H
#define INC_MDOPTS_H
#include "Constraints.h"
class MdOpts {
  public:
    MdOpts();
    MdOpts(ArgList&);

    double ScaleEE() const { return scaleEE_; }
    double ScaleNB() const { return scaleNB_; }
    double CutEE()   const { return cutEE_; }
    double CutNB()   const { return cutNB_; }
    Constraints::ShakeType Shake() const { return shakeType_; }
  private:
    Constraints::ShakeType shakeType_;
    double scaleEE_; ///< Global electrostatic 1-4 scaling factor
    double scaleNB_; ///< Global Lennard-Jones (VDW) 1-4 scaling factor
    double cutEE_;   ///< Nonbond electrostatic cutoff
    double cutNB_;   ///< Nonbond Lennard-Jones (VDW) cutoff
};
#endif
