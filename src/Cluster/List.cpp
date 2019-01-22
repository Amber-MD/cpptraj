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

int Cpptraj::Cluster::List::CreateCnumVsTime(DataSet_integer* ds, unsigned int maxFrames)
const
{
  if (ds == 0) {
    mprinterr("Internal Error: CreateCnumVsTime() called with null data set\n");
    return 1;
  }
  DataSet_integer& cnum_temp = static_cast<DataSet_integer&>( *ds );
  cnum_temp.Resize( maxFrames );
  // Make all clusters start at -1. This way cluster algorithms that
  // have noise points (i.e. no cluster assigned) will be distinguished.
  std::fill(cnum_temp.begin(), cnum_temp.end(), -1);

  for (cluster_iterator C = begincluster(); C != endcluster(); C++)
  {
    //mprinterr("Cluster %i:\n",CList->CurrentNum());
    int cnum = C->Num();
    // Loop over all frames in the cluster
    for (Node::frame_iterator frame = C->beginframe(); frame != C->endframe(); frame++)
    {
      //mprinterr("%i,",*frame);
      cnum_temp[ *frame ] = cnum;
    }
    //mprinterr("\n");
    //break;
  }
  return 0;
}
