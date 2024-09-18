#ifndef INC_ENERGY_ENERGYDECOMP_EWALD_H
#define INC_ENERGY_ENERGYDECOMP_EWALD_H
#include <vector>
#include "EwaldParams.h"
#include "../PairList.h" // For AtmType
class ExclusionArray;
class Frame;
class NonbondParmType;
namespace Cpptraj {
namespace Energy {
class EnergyDecomp_Ewald {
  public:
    EnergyDecomp_Ewald();
  private:
    typedef std::vector<int> Iarray;
    typedef std::vector<double> Darray;

    double ERFC(double);
    double adjust(double, double, double);
    void calcAdjust(double&, double&, PairList::AtmType const&, PairList::AtmType const&,
                    double, double, double);
    void ene_nb(double&, double&, double&,
                double, double, double, PairList::AtmType const&, PairList::AtmType const&);
    void ene_ewald_direct(double&, double&, double&, Frame const&,
                          PairList const&, ExclusionArray const&);

    EwaldParams EW_;
    Darray Charge_;       ///< Array of charges
    Darray Cparam_;       ///< Array of C6 coefficients for LJPME
    Iarray TypeIndices_;          ///< Hold atom type indices for selected atoms
    NonbondParmType const* NB_;   ///< Pointer to nonbonded parameters

    Timer t_direct_;
    Timer t_adjust_;
    Timer t_erfc_;

};
}
}
#endif
