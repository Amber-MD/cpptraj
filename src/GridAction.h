#ifndef INC_GRIDACTION_H
#define INC_GRIDACTION_H
#include "AtomMask.h"
#include "Vec3.h"
#include "Frame.h"
#include "DataSet_GridFlt.h"
#include "DataSetList.h"
#include "GridMover.h"
class Topology;
class ArgList;
class CoordinateInfo;
/// Class for setting up a grid within an action.
class GridAction {
  public:
    /// Indicate whether to apply an offset to coords before gridding.
    enum OffsetType { NO_OFFSET = 0, BOX_CENTER, MASK_CENTER };
    /// CONSTRUCTOR
    GridAction();
    /// \return List of keywords recognized by GridInit.
    static const char* HelpText;
    /// \return Set-up grid (added to given DataSetList) after processing keywords.
    DataSet_GridFlt* GridInit(const char*, ArgList&, DataSetList&);
#   ifdef MPI
    /// Perform any parallel initialization
    int ParallelGridInit(Parallel::Comm const&, DataSet_GridFlt*);
    /// \return the current Comm
    Parallel::Comm const& TrajComm() const { return trajComm_; }
#   endif
    /// Print information on given grid to STDOUT
    void GridInfo(DataSet_GridFlt const&) const;
    /// Perform any setup necessary for given Topology/CoordinateInfo
    int GridSetup(Topology const&, CoordinateInfo const&);

    /// Place atoms selected by given mask in given Frame on the given grid.
    inline void GridFrame(Frame const&, AtomMask const&, DataSet_GridFlt&) const;
    /// Move grid if necessary
    inline void MoveGrid(Frame const&, DataSet_GridFlt&);
    /// Anything needed to finalize the grid
    inline void FinishGrid(DataSet_GridFlt&) const;

    /// \return Type of offset to apply to coords before gridding.
    OffsetType GridOffsetType()  const { return gridOffsetType_;       }
    /// \return Mask to use for centering grid
    AtomMask const& CenterMask() const { return centerMask_; }
    /// \return Amount voxels should be incremented by
    float Increment()            const { return increment_;  }
  private:
    /// Determine move type if any based on arguments. Set centerMask_ for rms fit
    int determineMoveType(ArgList&, Cpptraj::GridMover::MoveType&, bool&);

    OffsetType gridOffsetType_;
    AtomMask centerMask_;
    float increment_;           ///< Set to -1 if negative, 1 if not.
    Cpptraj::GridMover mover_;  ///< Used to translate/rotate grid if needed.
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
  mover_.MoveGrid(currentFrame, centerMask_, grid);
}

/** Finish moving grid if necessary. */
void GridAction::FinishGrid(DataSet_GridFlt& grid) const {
  mover_.MoverFinish(grid);
}
#endif
