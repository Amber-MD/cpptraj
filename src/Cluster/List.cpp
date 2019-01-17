#include "List.h"
#include "../CpptrajStdio.h"

void Cpptraj::Cluster::List::PrintClusters() const {
  //mprintf("CLUSTER: %u clusters, %u frames.\n", clusters_.size(),
  //        FrameDistances().OriginalNframes() );
  mprintf("CLUSTER: %u clusters.\n", clusters_.size());
  for (cluster_iterator C = begincluster(); C != endcluster(); C++) {
    mprintf("\t%8i : ",C->Num());
    for (Node::frame_iterator fnum = C->beginframe();
                              fnum != C->endframe(); ++fnum)
      mprintf("%i,",(*fnum)+1);
    mprintf("\n");
  }
}

