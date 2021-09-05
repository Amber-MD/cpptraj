#include "AssignRefsToClusters.h"
#include "List.h"
#include "Node.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"
#include "../DataSetList.h"
#include "../Frame.h"

int Cpptraj::Cluster::AssignRefsToClusters( DataSetList const& refSets,
                                            std::string const& refmaskexpr,
                                            double refCut,
                                            bool useMass,
                                            DataSet_Coords& coords,
                                            List& CList )
{
  // Pre-center all reference coords at the origin. No need to store trans vectors.
  std::vector<Frame> refFrames;
  refFrames.reserve( refSets.size() );
  for (unsigned int idx = 0; idx != refSets.size(); idx++) {
    AtomMask rMask( refmaskexpr );
    DataSet_Coords_REF* REF_ds = (DataSet_Coords_REF*)refSets[idx];
    if ( REF_ds->Top().SetupIntegerMask( rMask, REF_ds->RefFrame() ) ) {
      mprintf("Warning: Could not set up mask for reference '%s'\n", REF_ds->legend());
      continue;
    }
    refFrames.push_back( Frame(REF_ds->RefFrame(), rMask) );
    refFrames.back().CenterOnOrigin( useMass );
  }
  // For each cluster, assign the reference name with the lowest RMSD
  // to the representative frame that is below the cutoff.
  AtomMask tMask( refmaskexpr );
  if (coords.Top().SetupIntegerMask( tMask )) {
    mprinterr("Error: Could not set up mask for assigning references.\n");
    return 1;
  }
  Frame TGT( coords.AllocateFrame(), tMask );
  unsigned int cidx = 0;
  for (List::cluster_it cluster = CList.begin();
                        cluster != CList.end(); ++cluster, ++cidx)
  {
    coords.GetFrame( cluster->BestRepFrame(), TGT, tMask );
    double minRms = TGT.RMSD_CenteredRef( refFrames[0], useMass );
    unsigned int minIdx = 0;
    for (unsigned int idx = 1; idx < refSets.size(); idx++) {
      double rms = TGT.RMSD_CenteredRef( refFrames[idx], useMass );
      if (rms < minRms) {
        minRms = rms;
        minIdx = idx;
      }
    }
    if (minRms < refCut) {
      //mprintf("DEBUG: Assigned cluster %i to reference \"%s\" (%g)\n", cidx,
      //        refSets[minIdx]->Meta().Name().c_str(), minRms);
      cluster->SetNameAndRms( refSets[minIdx]->Meta().Name(), minRms );
    } else {
      //mprintf("DEBUG: Cluster %i was closest to reference \"(%s)\" (%g)\n", cidx,
      //        refSets[minIdx]->Meta().Name().c_str(), minRms);
      cluster->SetNameAndRms( "(" + refSets[minIdx]->Meta().Name() + ")", minRms );
    }
  }
  return 0;
}
