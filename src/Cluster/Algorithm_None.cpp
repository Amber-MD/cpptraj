#include "Algorithm_None.h"
#include "../CpptrajStdio.h"
#include "../CpptrajFile.h"

using namespace Cpptraj::Cluster;

void Algorithm_None::Help() {
  mprintf("\t[nocluster]\n");
}

/** CONSTRUCTOR */
Algorithm_None::Algorithm_None() :
  Algorithm(NOCLUSTER)
{}

int Algorithm_None::Setup(ArgList& analyzeArgs) {
  return 0;
}

void Algorithm_None::Info() const {
  mprintf("\tNo Clustering.\n");
}

int Algorithm_None::DoClustering(List& clusters,
                                 Cframes const& framesToCluster,
                                 MetricArray& pmatrix)
{
  return 0;
}

void Algorithm_None::Results(CpptrajFile& outfile) const {
  outfile.Printf("#Algorithm: No clustering\n");
}
