#ifndef INC_ENERGY_EWALDCALC_DECOMP_H
#define INC_ENERGY_EWALDCALC_DECOMP_H
#include "EwaldCalc.h"
#include <vector>
namespace Cpptraj {
namespace Energy {
/// Abstract base class for decomposable Ewald calcs.
class EwaldCalc_Decomp : public EwaldCalc {
  public:
    typedef std::vector<double> Darray;

    EwaldCalc_Decomp() {}
    // virtual since inherited
    ~EwaldCalc_Decomp() {}

    Darray const& Atom_Elec() const { return atom_elec_; }
    Darray const& Atom_VDW() const { return atom_vdw_; }
  protected:
    /// For DEBUG, \return sum of given array
    static inline double sumArray(Darray const& arrayIn) {
      double sum = 0;
      for (Darray::const_iterator it = arrayIn.begin(); it != arrayIn.end(); ++it)
        sum += *it;
      return sum;
    }

    Darray atom_elec_; ///< Hold electrostatic contribution for each atom
    Darray atom_vdw_;  ///< Hold VDW contribution for each atom
};
}
}
#endif
