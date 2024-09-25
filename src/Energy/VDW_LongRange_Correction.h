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
    typedef std::vector<double> Darray;
    double Vdw_Recip_term_;
    int debug_;
    // Below variables are needed for per-atom decomp
    Iarray vdw_type_; ///< VDW type for each seleccted atom (#atoms)
    Iarray N_vdw_type_; ///< Count of atoms that have each VDW type index (#types)
    Darray atype_vdw_recip_terms_; ///< the nonbond interaction for each atom type (#types)
};
}
}
#endif
