#ifndef INC_REFERENCEACTION_H
#define INC_REFERENCEACTION_H
#include "DataSetList.h"
#include "DataSet_Coords.h"
/// Class that can be used by Actions to hold a COORDS DataSet to use as reference.
class ReferenceAction {
  public:
    ReferenceAction();
    ~ReferenceAction();
    /// Process all reference-related arguments, figure out reference mode.
    int InitRef(ArgList&, DataSetList const&, bool, bool);
#   ifdef MPI
    /// Set reference trajectory comm so SetupRef() can properly broadcast reference.
    void SetTrajComm(Parallel::Comm const& c) { trajComm_ = c; }
#   endif
    /// Set reference mask string.
    int SetRefMask(std::string const& m) { return refMask_.SetMaskString( m ); }
    /// \return String describing current reference mode.
    std::string RefModeString() const;
    /// Setup reference mask. Check that # selected reference atoms matches given # target atoms.
    int SetupRef(Topology const&, int);
    /// Peform necessary reference action based on mode
    inline void ActionRef(int, Frame const&);
    /// Store frame for previous if necessary
    inline void PreviousRef(Frame const&);
    /// \return Current entire reference frame.
    Frame const& CurrentReference() const { return refFrame_;    }
    /// \return Current reference frame, selected atoms only.
    Frame const& SelectedRef()      const { return selectedRef_; }
    /// \return Translation vector from origin to original ref center.
    Vec3 const& RefTrans()          const { return refTrans_;    }
    /// \return Help text
    static const char* Help() { return help_.c_str(); }
    /// \return Pointer to reference COORDS topology if possible.
    Topology* RefCrdTopPtr() const {
      if (refCrd_ != 0) return refCrd_->TopPtr();
      return 0;
    }
  private:
    /// Set up refFrame_ and selectedRef_ coordinates base on given frame and current mask.
    void SelectRefAtoms(Frame const&);
    /// Set up ref mask for given topology. Allocate space for selected ref atoms.
    int SetupRefMask(Topology const&);
    /// Modes: FIRST=first frame, FRAME=given frame, TRAJ=reference traj
    enum RefModeType { FIRST = 0, FRAME, TRAJ };

    RefModeType refMode_;    ///< Reference mode.
    DataSet_Coords* refCrd_; ///< Reference COORDS DataSet.
    DataSet_Coords* traj_;   ///< Allocated reference COORDS if not present in DataSetList
    AtomMask refMask_;       ///< Atoms to use from refCrd_.
    Frame refFrame_;         ///< Current reference frame.
    Frame selectedRef_;      ///< Atoms from reference frame selected by refMask_.
    Vec3 refTrans_;          ///< If fitting, translation from origin to original ref center.
    static std::string help_;///< Help text.
    bool previous_;          ///< True if current reference is previous frame (only RMSD now)
    bool needsSetup_;        ///< True if ref from COORDS needs to be set up during SetupRef()
    bool fitRef_;            ///< If true, move reference to origin for RMS fitting
    bool useMass_;           ///< (If fitRef_) If true, move COM, otherwise geometric center.
#   ifdef MPI
    Parallel::Comm trajComm_; ///< Comm across trajectory for synchronizing ref from master.
#   endif
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
// ReferenceAction::ActionRef()
void ReferenceAction::ActionRef(int frameNum, Frame const& frameIn)
{
  if (refMode_ == FIRST) {
    SelectRefAtoms( frameIn );
    refMode_ = FRAME;
  } else if (refMode_ == TRAJ) {
    refCrd_->GetFrame( frameNum, refFrame_ );
    selectedRef_.SetCoordinates(refFrame_, refMask_);
    if (fitRef_)
      refTrans_ = selectedRef_.CenterOnOrigin(useMass_);
  }
}

void ReferenceAction::PreviousRef(Frame const& frameIn) {
  if (previous_) SelectRefAtoms( frameIn );
}
#endif
