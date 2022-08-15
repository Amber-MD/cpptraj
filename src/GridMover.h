#ifndef INC_GRIDMOVER_H
#define INC_GRIDMOVER_H
#include "Frame.h"
// Forward declares
//class AtomMask;
//class Frame;
//class Matrix_3x3;
class Topology;
namespace Cpptraj {
/// Helper class for moving/rotating a GridBin
class GridMover {
  public:
    /// Indicate where grid should be located
    enum MoveType { NO_MOVE = 0, TO_BOX_CTR, TO_MASK_CTR, RMS_FIT };

    /// CONSTRUCTOR
    GridMover();
    /// Initialize (move type, do final x-align)
    int MoverInit(MoveType, bool);
#   ifdef MPI
    /// Set the trajectory comm
    int MoverSetComm(Parallel::Comm const&);
#   endif
    /// Setup Tgt and Ref frames for given topology and mask
    int MoverSetup(Topology const&, AtomMask const&);
  private:
    /// Set first frame selected coords (tgt_) and original grid unit cell vectors (tgtUcell_).
    int SetTgt(Frame const&, Matrix_3x3 const&);

    MoveType gridMoveType_; ///< The move type
    Frame tgt_;             ///< For RMS_FIT, first frames selected coordinates
    Matrix_3x3 tgtUcell_;   ///< For RMS_FIT, original grid unit cell vectors
    Frame ref_;             ///< For RMS_FIT, current frames selected coordinates
    bool firstFrame_;       ///< For RMS_FIT, true if this is the first frame (no fit needed)
    bool x_align_;          ///< For RMS_FIT, if true ensure grid is X-aligned in FinishGrid().
#   ifdef MPI
    Parallel::Comm trajComm_; ///< For RMS_FIT used to ensure all procs use same reference
#   endif
};
} // end namespace Cpptraj
#endif
