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

/** Print info to stdout. */
void GridMover::MoverInfo(AtomMask const& maskIn) const {
  if (gridMoveType_ == NO_MOVE)
    mprintf("\tGrid will not move.\n");
  else if (gridMoveType_ == TO_BOX_CTR)
    mprintf("\tGrid will be kept centered at the box center.\n");
  else if (gridMoveType_ == TO_MASK_CTR)
    mprintf("\tGrid will be kept centered on atoms in mask [%s]\n",
            maskIn.MaskString());
  else if (gridMoveType_ == RMS_FIT) {
    mprintf("\tGrid will be RMS-fit using atoms in mask [%s]\n",
            maskIn.MaskString());
    if (x_align_)
      mprintf("\tGrid will be realigned with Cartesian axes after binning is complete.\n");
    else
      mprintf("\tGrid will not be realigned with Cartesian axes after binning is complete.\n");
  }
}

#ifdef MPI
/** Set the trajectory comm for RMS fit */
void GridMover::MoverSetComm(Parallel::Comm const& commIn) {
  trajComm_ = commIn;
}
#endif

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

/** Move/reorient grid if necessary. */ // TODO inline?
void GridMover::MoveGrid(Frame const& currentFrame, AtomMask const& maskIn, DataSet_3D& grid)
{
  if (gridMoveType_ == TO_BOX_CTR)
    grid.SetGridCenter( currentFrame.BoxCrd().Center() );
  else if (gridMoveType_ == TO_MASK_CTR)
    grid.SetGridCenter( currentFrame.VGeometricCenter( maskIn ) );
  else if (gridMoveType_ == RMS_FIT) {
    grid.SetGridCenter( currentFrame.VGeometricCenter( maskIn ) );
#   ifdef MPI
    // Ranks > 0 still need to do the rotation on the first frame.
    bool doRotate = true;
    if (firstFrame_) {
      SetTgt(currentFrame, grid.Bin().Ucell(), maskIn);
      if (trajComm_.Rank() == 0)
        doRotate = false;
      firstFrame_ = false;
    }
    if (doRotate) {
      // Want to rotate to coordinates in current frame. Make them the ref.
      ref_.SetFrame( currentFrame, maskIn );
      // Reset to original grid.
      grid.Assign_Grid_UnitCell( tgtUcell_ );
      // Do not want to modify original coords. Make a copy.
      Frame tmpTgt( tgt_ );
      // Rot will contain rotation from original grid to current frame.
      Matrix_3x3 Rot;
      Vec3 T1, T2;
      tmpTgt.RMSD( ref_, Rot, T1, T2, false );
      grid.Rotate_3D_Grid( Rot );
    }
#   else
    if (firstFrame_) {
      SetTgt(currentFrame, grid.Bin().Ucell(), maskIn);
      //grid.SetGridCenter( tgt_.VGeometricCenter( 0, tgt_.Natom() ) );
      firstFrame_ = false;
    } else {
      //grid.SetGridCenter( currentFrame.VGeometricCenter( maskIn ) );
      // Want to rotate to coordinates in current frame. Make them the ref.
      ref_.SetFrame( currentFrame, maskIn );
      // Reset to original grid.
      grid.Assign_Grid_UnitCell( tgtUcell_ );
      // Do not want to modify original coords. Make a copy.
      Frame tmpTgt( tgt_ );
      // Rot will contain rotation from original grid to current frame.
      Matrix_3x3 Rot;
      Vec3 T1, T2;
      tmpTgt.RMSD( ref_, Rot, T1, T2, false );
      grid.Rotate_3D_Grid( Rot );
      //tgt_.SetFrame( currentFrame, maskIn );
    }
#   endif
  }
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
