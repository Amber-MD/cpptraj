#ifndef INC_ENSEMBLENAVIGATOR_H
#define INC_ENSEMBLENAVIGATOR_H
#include "TrajinList.h"
#include "TrajFrameIndex.h"
/// Can be used to randomly access a collection of EnsembleIn trajectories.
/** Currently this is only used in CpptrajState::RunParaEnsemble() as an
  * analogue to DataSet_Coords_TRJ in CpptrajState::RunParallel() for
  * facilitating access to a collection of ensembles.
  */
class EnsembleNavigator {
  public:
    EnsembleNavigator() : currentEns_(0), FirstParm_(0) {}
    int AddEnsembles(TrajinList::ensemble_it const&, TrajinList::ensemble_it const&);
    inline int GetEnsemble(int, FrameArray&, FramePtrArray&);
    EnsembleIn* CurrentEns()                   { return currentEns_;   }
    Topology* FirstParm()                      { return FirstParm_;    }
    CoordinateInfo const& EnsCoordInfo() const { return ensCoordInfo_; }
    TrajFrameIndex const& IDX()          const { return IDX_;          }
  private:
    typedef std::vector<EnsembleIn*> Earray;
    Earray Ensembles_;            ///< Array of input ensembles
    CoordinateInfo ensCoordInfo_; ///< Coordinate info for all ensembles
    TrajFrameIndex IDX_;          ///< Used to convert global index to individual index
    EnsembleIn* currentEns_;      ///< Currently open ensemble
    Topology* FirstParm_;         ///< Topology associated with first ensemble
};
/// Get ensemble set specified by global set index
int EnsembleNavigator::GetEnsemble(int set, FrameArray& FrameEnsemble,
                                   FramePtrArray& SortedFrames)
{
  int internalIdx = IDX_.FindIndex( set );
  // If desired ensemble is different than current, open desired ensemble.
  if (IDX_.TrajHasChanged()) {
    if (currentEns_ != 0) currentEns_->EndEnsemble();
    currentEns_ = Ensembles_[ IDX_.CurrentTrajNum() ];
    // NOTE: Need to check for reallocation? TODO better error check
    if (currentEns_->BeginEnsemble()) return 1;
  }
  // Read the frame
  currentEns_->ReadEnsemble( internalIdx, FrameEnsemble, SortedFrames );
  return 0;
}
#endif
