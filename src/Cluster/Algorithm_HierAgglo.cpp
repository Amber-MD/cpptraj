#include "Algorithm_HierAgglo.h"
#include "../CpptrajStdio.h"

Cpptraj::Cluster::Algorithm_HierAgglo::Algorithm_HierAgglo() :
  nclusters_(-1),
  epsilon_(-1.0),
  linkage_(AVERAGELINK)//,
//  includeSievedFrames_(false)
{}

void Cpptraj::Cluster::Algorithm_HierAgglo::Help() {
  mprintf("\t[hieragglo [epsilon <e>] [clusters <n>] [linkage|averagelinkage|complete]\n"
          "\t  [epsilonplot <file>] [includesieved_cdist]]\n");
}

static const char* LinkageString[] = {
  "single-linkage", "average-linkage", "complete-linkage"
};

int Cpptraj::Cluster::Algorithm_HierAgglo::Setup(ArgList& analyzeArgs) {
  nclusters_ = analyzeArgs.getKeyInt("clusters", -1);
  epsilon_ = analyzeArgs.getKeyDouble("epsilon", -1.0);
  if (analyzeArgs.hasKey("linkage"))             linkage_ = SINGLELINK;
  else if (analyzeArgs.hasKey("averagelinkage")) linkage_ = AVERAGELINK;
  else if (analyzeArgs.hasKey("complete"))       linkage_ = COMPLETELINK;
  else linkage_ = AVERAGELINK; // DEFAULT linkage
//  includeSievedFrames_ = analyzeArgs.hasKey("includesieved_cdist");
  std::string epsilonPlot = analyzeArgs.GetStringKey("epsilonplot");
  if (!epsilonPlot.empty()) {
    if (eps_v_n_.OpenWrite( epsilonPlot )) return 1;
    eps_v_n_.Printf("%-12s %12s\n", "#Epsilon", "Nclusters");
  }
  // Determine finish criteria. If nothing specified default to 10 clusters.
  if (nclusters_==-1 && epsilon_==-1.0) {
    mprintf("Warning: cluster: Neither target # of clusters nor epsilon given.\n");
    nclusters_ = 10;
    mprintf("Warning: cluster: Defaulting to %i clusters.\n", nclusters_);
  }
  return 0;
}

void Cpptraj::Cluster::Algorithm_HierAgglo::Info() const {
    mprintf("\tHierarchical Agglomerative:");
  if (nclusters_ != -1)
    mprintf(" %i clusters,",nclusters_);
  if (epsilon_ != -1.0)
    mprintf(" epsilon %.3f,",epsilon_);
  mprintf(" %s.\n", LinkageString[linkage_]);
  if (eps_v_n_.IsOpen())
    mprintf("\tWriting epsilon vs # clusters to '%s'\n", eps_v_n_.Filename().full());
  /*if (includeSievedFrames_)
    mprintf("\tSieved frames will be included in final cluster distance calculation.\n"
            "Warning: 'includesieved_cdist' may be very slow.\n");
  else
    mprintf("\tSieved frames will not be included in final cluster distance calculation.\n");*/
}

