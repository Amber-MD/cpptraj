#ifndef INC_REFERENCEACTION_H
#define INC_REFERENCEACTION_H
#include "Trajin_Single.h"
#include "FrameList.h" // ReferenceFrame
/// Class that can be used by Actions to hold reference structure/trajectory.
class ReferenceAction {
  public:
    enum RefModeType { UNKNOWN_REF=0, FIRST, REFFRAME, REFTRAJ };
    ReferenceAction() {}
    /// Set up ref mask for given topology. Allocate space for selected ref atoms.
    int SetRefMask(Topology const&, const char*);
    /// Set up selected ref coordinates base on given frame and ref mask.
    void SetRefStructure(Frame const& frameIn, bool nofitIn, bool useMassIn) {
      RefFrame_ = frameIn;
      SelectedRef_.SetCoordinates( RefFrame_, RefMask_ );
      if (!nofitIn)
        refTrans_ = SelectedRef_.CenterOnOrigin( useMassIn );
    }
    /// Initialize
    int InitRef(bool, bool, bool, bool, std::string const&, ReferenceFrame&,
                Topology*, std::string const&, ArgList&, const char*);
    /// Setup
    int SetupRef(Topology const&, int, const char*);
    /// Action
    void ActionRef(Frame const& frameIn, bool nofitIn, bool useMassIn) {
      if (refmode_ == FIRST) {
        SetRefStructure( frameIn, nofitIn, useMassIn );
        refmode_ = REFFRAME;
      } else if (refmode_ == REFTRAJ) {
        RefTraj_.GetNextFrame( RefFrame_ );
        SelectedRef_.SetCoordinates(RefFrame_, RefMask_);
        if (!nofitIn)
          refTrans_ = SelectedRef_.CenterOnOrigin(useMassIn);
      }
    }

    bool Previous() const { return previous_; }
        
  private:
    RefModeType refmode_;
    Frame RefFrame_;
    Frame SelectedRef_;
    AtomMask RefMask_;
    Vec3 refTrans_;
    Trajin_Single RefTraj_;
    bool previous_;
};
#endif
