#ifndef INC_ENERGY_EWALDCALC_H
#define INC_ENERGY_EWALDCALC_H
#include "../Timer.h"
class AtomMask;
class Box;
class EwaldOptions;
class ExclusionArray;
class Frame;
class PairList;
class Topology;
namespace Cpptraj {
namespace Energy {
/// Abstract base class for Ewald calcs
class EwaldCalc {
  public:
    EwaldCalc() {}
    // virtual since inherited
    virtual ~EwaldCalc() {}

    virtual int Init(Box const&, EwaldOptions const&, int) = 0;
    virtual int Setup(Topology const&, AtomMask const&) = 0;
    virtual int CalcNonbondEnergy(Frame const&, AtomMask const&,
                                  PairList const&, ExclusionArray const&,
                                  double&, double&) = 0;
    virtual void Timing(double) const = 0;
  protected:
    Timer t_total_;
    Timer t_direct_;
};
}
}
#endif
