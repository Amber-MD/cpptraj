#include "GridMover.h"
#include "CpptrajStdio.h"
#include "DataSet_3D.h"
#include "Topology.h"

using namespace Cpptraj;

/** CONSTRUCTOR. */
GridMover::GridMover() :
  gridMoveType_(NO_MOVE),
  firstFrame_(true),
  x_align_(false)
{}

/** Set GridMover options. */
int GridMover::MoverInit(MoveType typeIn, bool x_alignIn) {
  gridMoveType_ = typeIn;
  x_align_ = x_alignIn;
  firstFrame_ = true;
  return 0;
}

/** Setup GridMover */
int GridMover::MoverSetup(Topology const& topIn, AtomMask const& maskIn) {
  // Set up frames if needed
  if (gridMoveType_ == RMS_FIT) {
    tgt_.SetupFrameFromMask(maskIn, topIn.Atoms());
    ref_ = tgt_;
    //firstFrame_ = true;
  }
  return 0;
}

/** Set the coordinates of the first frame. Set original grid unit cell vectors. */
int GridMover::SetTgt(Frame const& frameIn, Matrix_3x3 const& gridUcell, AtomMask const& maskIn)
{
  tgt_.SetFrame( frameIn, maskIn );
  tgtUcell_ = gridUcell;
# ifdef MPI
  // Ensure all processes are using the same reference. Just broadcast the coords.
  trajComm_.MasterBcast( tgt_.xAddress(), tgt_.size(), MPI_DOUBLE );
  // Ensure all processes have the same unit cell vecs
  trajComm_.MasterBcast( tgtUcell_.Dptr(), 9, MPI_DOUBLE );
  //rprintf("DEBUG: Ucell0: %f %f %f %f %f %f %f %f %f\n", tgtUcell_[0], tgtUcell_[1], tgtUcell_[2], tgtUcell_[3], tgtUcell_[4], tgtUcell_[5], tgtUcell_[6], tgtUcell_[7], tgtUcell_[8]);
# endif
  return 0;
}

/** Any final actions to grid. */
void GridMover::MoverFinish(DataSet_3D& grid) const {
  //rprintf("DEBUG: Final Grid origin: %f %f %f\n", grid.Bin().GridOrigin()[0], grid.Bin().GridOrigin()[1], grid.Bin().GridOrigin()[2]);
  if (x_align_) {
    if (!grid.Bin().IsXalignedGrid()) {
      mprintf("\tEnsuring grid '%s' is X-aligned.\n", grid.legend());
      grid.Xalign_3D_Grid();
    }
  }
}
