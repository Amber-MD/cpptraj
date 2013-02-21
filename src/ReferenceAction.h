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
    void SetRefStructure(Frame const& frameIn, bool fitIn, bool useMassIn) {
      refFrame_ = frameIn;
      selectedRef_.SetCoordinates( refFrame_, refMask_ );
      if (fitIn)
        refTrans_ = selectedRef_.CenterOnOrigin( useMassIn );
    }
    /// Initialize
    int InitRef(bool, bool, bool, bool, std::string const&, ReferenceFrame&,
                Topology*, std::string const&, ArgList&, const char*);
    /// Setup
    int SetupRef(Topology const&, int, const char*);
    /// Action
    void ActionRef(Frame const& frameIn, bool fitIn, bool useMassIn) {
      if (refmode_ == FIRST) {
        SetRefStructure( frameIn, fitIn, useMassIn );
        refmode_ = REFFRAME;
      } else if (refmode_ == REFTRAJ) {
        refTraj_.GetNextFrame( refFrame_ );
        selectedRef_.SetCoordinates(refFrame_, refMask_);
        if (fitIn)
          refTrans_ = selectedRef_.CenterOnOrigin(useMassIn);
      }
    }

    bool Previous()             const { return previous_;           } 
    const char* RefModeString() const { return modeString_.c_str(); }
    Frame const& RefFrame()     const { return refFrame_;           }
    Frame const& SelectedRef()  const { return selectedRef_;        }
    Vec3 const& RefTrans()      const { return refTrans_;           }
  private:
    RefModeType refmode_;
    Frame refFrame_;
    Frame selectedRef_;
    AtomMask refMask_;
    Vec3 refTrans_;
    Trajin_Single refTraj_;
    bool previous_;
    std::string modeString_;
};
#endif
