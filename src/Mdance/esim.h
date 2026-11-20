#ifndef INC_CPPTRAJ_MDANCE_ESIM_H
#define INC_CPPTRAJ_MDANCE_ESIM_H
#ifdef HAS_EIGEN
#include "types.h"
namespace Cpptraj {
namespace Mdance {
Cpptraj::Mdance::MD::Indices genSimIdx(const ArrayXd& cTotal, int nObjects, MD::Threshold& cThreshold, int wt);
// wFactor = 0 selects the "fraction" weight function (the reference default);
// any other value n selects the power_n weights.
Cpptraj::Mdance::MD::Counters calculateCounters(const ArrayXd& cTotal, int nObjects, MD::Threshold& cThreshold, int wFactor = 0);
}
}
#endif /* HAS_EIGEN */
#endif
