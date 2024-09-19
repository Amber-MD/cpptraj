#ifndef INC_ENERGY_CALC_PME_H
#define INC_ENERGY_CALC_PME_H
#include "../ExclusionArray.h"
#include "../helpme_standalone.h"
#include "../PairList.h"
#include "../PairListEngine_Ewald_LJLR.h"
class AtomMask;
class Box;
class EwaldOptions;
class Frame;
class Topology;
namespace Cpptraj {
namespace Energy {
class Calc_PME {
  public:
    Calc_PME();
    /// Init with Box, EwaldOptions and debug level
    int Init(Box const&, EwaldOptions const&, int);
    int Setup(Topology const&, AtomMask const&);
    int CalcNonbondEnergy(Frame const&, AtomMask const&, double&, double&);

  private:
    PairListEngine_Ewald_LJLR<double> NBengine_;
    PMEInstanceD pme_object_;
    PairList pairList_;
    ExclusionArray Excluded_;
};
}
}
#endif
