#ifndef INC_CPPTRAJ_MDANCE_SCORES_H
#define INC_CPPTRAJ_MDANCE_SCORES_H
#include "types.h"
#ifdef HAS_EIGEN
namespace Cpptraj {
namespace Mdance {
double calinskiHarabaszScore(const ArrayXXd& data, const VectorXi& labels);
double daviesBouldinScore(const ArrayXXd& data, const VectorXi& labels);
}
}
#endif
#endif
