#ifndef INC_STATS_REDUCE_H
#define INC_STATS_REDUCE_H
#ifdef MPI
#include <vector>
#include "Parallel.h"
#include "OnlineVarT.h"
namespace Cpptraj {
/// Used to reduce an array of Stats<double> down to master rank
int Stats_Reduce(Parallel::Comm const&, std::vector< Stats<double> >&, unsigned long&);
}
#endif /* MPI */
#endif
