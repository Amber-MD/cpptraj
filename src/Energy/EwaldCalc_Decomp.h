#ifndef INC_ENERGY_EWALDCALC_DECOMP_H
#define INC_ENERGY_EWALDCALC_DECOMP_H
#include "EwaldCalc.h"
namespace Cpptraj {
namespace Energy {
/// Abstract base class for decomposable Ewald calcs.
class EwaldCalc_Decomp : public EwaldCalc {
  public:
    typedef std::vector<double> Darray;

    EwaldCalc_Decomp() {}
    // virtual since inherited
    ~EwaldCalc_Decomp() {}

    virtual int CalcDecomposedNonbondEnergy(Frame const&, AtomMask const&,
                                            PairList const& pairList, ExclusionArray const& Excluded,
                                            double&, double&, Darray&, Darray&) = 0;

    /// Just call the CalcDecomposedNonbondEnergy() routine with temp arrays.
    int CalcNonbondEnergy(Frame const& frameIn, AtomMask const& maskIn,
                          PairList const& pairList, ExclusionArray const& Excluded,
                         double& e_elec, double& e_vdw)
    {
      Darray atom_elec, atom_vdw;
      return CalcDecomposedNonbondEnergy(frameIn, maskIn, pairList, Excluded, e_elec, e_vdw, atom_elec, atom_vdw);
    }
};
}
}
#endif
