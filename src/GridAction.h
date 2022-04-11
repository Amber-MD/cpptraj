#ifndef INC_GRIDACTION_H
#define INC_GRIDACTION_H
#include "AtomMask.h"
#include "Vec3.h"
#include "Frame.h"
#include "DataSet_GridFlt.h"
#include "DataSetList.h"
class Topology;
class ArgList;
class CoordinateInfo;
/// Class for setting up a grid within an action.
class GridAction {
  public:
    /// Indicate whether to apply an offset to coords before gridding.
    enum OffsetType { NO_OFFSET = 0, BOX_CENTER, MASK_CENTER };
    /// Indicate where grid should be located
    enum MoveType { NO_MOVE = 0, TO_BOX_CTR, TO_MASK_CTR, RMS_FIT };
    /// CONSTRUCTOR
    GridAction();
    /// \return List of keywords recognized by GridInit.
    static const char* HelpText;
    /// \return Set-up grid (added to given DataSetList) after processing keywords.
    DataSet_GridFlt* GridInit(const char*, ArgList&, DataSetList&);
#   ifdef MPI
    /// Perform any parallel initialization
    int ParallelGridInit(Parallel::Comm const&, DataSet_GridFlt*);
#   endif
    /// Print information on given grid to STDOUT
    void GridInfo(DataSet_GridFlt const&);
    /// Perform any setup necessary for given Topology/CoordinateInfo
    int GridSetup(Topology const&, CoordinateInfo const&);
    /// Place atoms selected by given mask in given Frame on the given grid.
    inline void GridFrame(Frame const&, AtomMask const&, DataSet_GridFlt&) const;
    /// Move grid if necessary
    inline void MoveGrid(Frame const&, DataSet_GridFlt&);
    /// Anything needed to finalize the grid
    void FinishGrid(DataSet_GridFlt&) const;
    /// \return Type of offset to apply to coords before gridding.
    OffsetType GridOffsetType()  const { return gridOffsetType_;       }
    /// \return Mask to use for centering grid
    AtomMask const& CenterMask() const { return centerMask_; }
    /// \return Amount voxels should be incremented by
    float Increment()            const { return increment_;  }
  private:
    /// Set first frame selected coords (tgt_) and original grid unit cell vectors (tgtUcell_).
    int SetTgt(Frame const&, Matrix_3x3 const&);

    OffsetType gridOffsetType_;
    MoveType gridMoveType_;
    AtomMask centerMask_;
    float increment_;     ///< Set to -1 if negative, 1 if not.
    Frame tgt_;           ///< For MoveType RMS_FIT, first frames selected coordinates
    Matrix_3x3 tgtUcell_; ///< For MoveType RMS_FIT, original grid unit cell vectors
    Frame ref_;           ///< For MoveType RMS_FIT, current frames selected coordinates
    bool firstFrame_;     ///< For MoveType RMS_FIT, true if this is the first frame (no fit needed)
    bool x_align_;        ///< For MoveType RMS_FIT, if true ensure grid is X-aligned in FinishGrid().
#   ifdef MPI
    Parallel::Comm trajComm_;
#   endif
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
void GridAction::GridFrame(Frame const& currentFrame, AtomMask const& mask, 
                           DataSet_GridFlt& grid)
const
{
  if (gridOffsetType_ == BOX_CENTER) {
    Vec3 offset = currentFrame.BoxCrd().Center();
    for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
      grid.Increment( Vec3(currentFrame.XYZ(*atom)) - offset, increment_ );
  } else if (gridOffsetType_ == MASK_CENTER) {
    Vec3 offset = currentFrame.VGeometricCenter( centerMask_ );
    for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
      grid.Increment( Vec3(currentFrame.XYZ(*atom)) - offset, increment_ );
  } else {// mode_==NO_OFFSET
    for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
      grid.Increment( currentFrame.XYZ(*atom), increment_ );
  }
}

/** Move/reorient grid if necessary. */
void GridAction::MoveGrid(Frame const& currentFrame, DataSet_GridFlt& grid)
{
  if (gridMoveType_ == TO_BOX_CTR)
    grid.SetGridCenter( currentFrame.BoxCrd().Center() );
  else if (gridMoveType_ == TO_MASK_CTR)
    grid.SetGridCenter( currentFrame.VGeometricCenter( centerMask_ ) );
  else if (gridMoveType_ == RMS_FIT) {
    grid.SetGridCenter( currentFrame.VGeometricCenter( centerMask_ ) );
#   ifdef MPI
    // Ranks > 0 still need to do the rotation on the first frame.
    bool doRotate = true;
    if (firstFrame_) {
      SetTgt(currentFrame, grid.Bin().Ucell());
      if (trajComm_.Rank() == 0)
        doRotate = false;
      firstFrame_ = false;
    }
    if (doRotate) {
      // Want to rotate to coordinates in current frame. Make them the ref.
      ref_.SetFrame( currentFrame, centerMask_ );
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
      SetTgt(currentFrame, grid.Bin().Ucell());
      //grid.SetGridCenter( tgt_.VGeometricCenter( 0, tgt_.Natom() ) );
      firstFrame_ = false;
    } else {
      //grid.SetGridCenter( currentFrame.VGeometricCenter( centerMask_ ) );
      // Want to rotate to coordinates in current frame. Make them the ref.
      ref_.SetFrame( currentFrame, centerMask_ );
      // Reset to original grid.
      grid.Assign_Grid_UnitCell( tgtUcell_ );
      // Do not want to modify original coords. Make a copy.
      Frame tmpTgt( tgt_ );
      // Rot will contain rotation from original grid to current frame.
      Matrix_3x3 Rot;
      Vec3 T1, T2;
      tmpTgt.RMSD( ref_, Rot, T1, T2, false );
      grid.Rotate_3D_Grid( Rot );
      //tgt_.SetFrame( currentFrame, centerMask_ );
    }
#   endif
  }
}
#endif
