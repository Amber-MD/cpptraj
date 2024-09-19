#ifndef INC_ENERGY_VDW_LONGRANGE_CORRECTION_H
#define INC_ENERGY_VDW_LONGRANGE_CORRECTION_H
#include <vector>
#include "../Constants.h"
class AtomMask;
class Topology;
namespace Cpptraj {
namespace Energy {
class VDW_LongRange_Correction {
  public:
    VDW_LongRange_Correction();
    void SetDebug(int);
    int Setup_VDW_Correction(Topology const&, AtomMask const&);
    /// \return Full VDW long range correction from cutoff and volume
    double Vdw_Correction(double cutoff_, double volume) const {
      double prefac = Constants::TWOPI / (3.0*volume*cutoff_*cutoff_*cutoff_);
      double e_vdwr = -prefac * Vdw_Recip_term_;
      //if (debug_ > 0) mprintf("DEBUG: Vdw correction %20.10f\n", e_vdwr);
      return e_vdwr;
    }
  private:
    typedef std::vector<int> Iarray;
    double Vdw_Recip_term_;
    int debug_;
};
}
}
#endif
